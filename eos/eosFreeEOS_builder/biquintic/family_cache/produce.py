"""Build an audited, standalone FreeEOS potential-jet family from source tables."""

import argparse
import hashlib
import json
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from freeeos_potential import CRAD, LN10, checksum, construct_twists, source_constraints, valid_cells


ROOT = Path(__file__).resolve().parents[5]
SOURCE_DIR = ROOT / "mesa/data/eosFreeEOS_data"
DEFAULT_OUTPUT = ROOT / "freeeos_biquintic_progress/experiments/062_family_cache/family"
NAME = re.compile(r"mesa-FreeEOS_(\d+)z(\d+)x(_CO)?\.data$")
HEADER_FIELDS = ("version", "X", "Z", "num_logT", "logT_min", "logT_max", "delta_logT",
                 "num_logQ", "logQ_min", "logQ_max", "delta_logQ")
TOLERANCE = 1e-10


class SourceAuditError(ValueError):
    def __init__(self, message, audit):
        super().__init__(message)
        self.audit = audit


def source_identity(source):
    match = NAME.fullmatch(source.name)
    if match is None:
        raise ValueError(f"Unexpected table name: {source.name}")
    return int(match[2]) / 100, int(match[1]) / 100, bool(match[3])


def read_audited(source):
    expected_x, expected_z, co_variant = source_identity(source)
    digest = hashlib.sha256()
    blocks, labels = [], []
    header = None
    expect_label = False
    with source.open("rb") as stream:
        for line_number, raw in enumerate(stream, 1):
            digest.update(raw)
            line = raw.decode("ascii").strip()
            if line_number == 2:
                header = np.fromstring(line, sep=" ")
                if header.shape != (11,) or not np.all(np.isfinite(header)) or header[0] != 51:
                    raise ValueError("Header is not finite FreeEOS format 51")
                if abs(header[1] - expected_x) > TOLERANCE or abs(header[2] - expected_z) > TOLERANCE:
                    raise ValueError("Composition header differs from filename")
                continue
            if line == "logQ = logRho - 2*logT + 12":
                if expect_label:
                    raise ValueError("Missing Q label")
                expect_label = True
                continue
            if expect_label:
                if not line:
                    continue
                parts = line.split()
                if len(parts) != 1:
                    raise ValueError("Q label is not a single number")
                labels.append(float(parts[0]))
                blocks.append([])
                expect_label = False
                continue
            parts = line.split()
            if len(parts) == 22 and parts[0][:1] in "+-0123456789":
                if not blocks:
                    raise ValueError("Source record precedes Q label")
                row = np.fromstring(line, sep=" ")
                if row.shape != (22,):
                    raise ValueError("Malformed source record")
                blocks[-1].append(row)
    if expect_label or len(blocks) < 2 or min(map(len, blocks)) < 2:
        raise ValueError("Incomplete printed Q blocks")
    lengths = [len(block) for block in blocks]
    if any(length != lengths[0] for length in lengths):
        raise ValueError("Printed Q blocks do not form a rectangle")
    rows = np.asarray(blocks)
    logq = np.asarray(labels)
    logt = rows[0, :, 0].copy()
    if not np.all(np.isfinite(logq)) or not np.all(np.isfinite(logt)):
        raise ValueError("Nonfinite printed grid coordinate")
    delta_q, delta_t = np.diff(logq), np.diff(logt)
    if np.any(delta_q <= 0) or np.any(delta_t <= 0):
        raise ValueError("Printed coordinates are not strictly increasing")
    if np.max(abs(delta_q - delta_q[0])) > TOLERANCE or np.max(abs(delta_t - delta_t[0])) > TOLERANCE:
        raise ValueError("Printed coordinates are not a uniform native grid")
    observed = np.asarray([51, expected_x, expected_z, len(logt), logt[0], logt[-1], delta_t[0],
                           len(logq), logq[0], logq[-1], delta_q[0]])
    anomalies = {name: {"declared": float(header[index]), "observed": float(observed[index])}
                 for index, name in enumerate(HEADER_FIELDS)
                 if abs(header[index] - observed[index]) > TOLERANCE}
    audit = {
        "source_sha256": digest.hexdigest(), "source_filename": source.name,
        "X": expected_x, "Z": expected_z, "co_variant": co_variant,
        "declared_header": dict(zip(HEADER_FIELDS, map(float, header))),
        "observed_grid": dict(zip(HEADER_FIELDS, map(float, observed))),
        "header_anomalies": anomalies, "printed_Q_blocks": len(logq),
        "printed_records": int(rows.shape[0] * rows.shape[1]),
    }
    temperature_error = abs(rows[:, :, 0] - logt[None, :])
    audit["max_temperature_coordinate_error"] = float(np.max(temperature_error))
    if np.max(temperature_error) > TOLERANCE:
        raise SourceAuditError("Temperature records differ between Q blocks", audit)
    density_error = rows[:, :, 18] - (logq[:, None] + 2 * logt[None, :] - 12)
    known_exclusion = (abs(expected_x - .7) < TOLERANCE and abs(expected_z - .02) < TOLERANCE
                       and not co_variant and len(logq) > 520)
    if known_exclusion:
        if not np.all(abs(density_error[520] + 1e-4) <= TOLERANCE):
            raise SourceAuditError("Known X=.70/Z=.02 density exclusion changed", audit)
        density_error[520] = 0
    bad_density = np.argwhere(abs(density_error) > TOLERANCE)
    audit["max_density_coordinate_error_excluding_known_block"] = float(np.max(abs(density_error)))
    audit["density_coordinate_mismatch_count"] = len(bad_density)
    audit["density_coordinate_mismatch_examples"] = [
        {"iq_zero_based": int(index_q), "it_zero_based": int(index_t),
         "signed_delta_logrho": float(density_error[index_q, index_t])}
        for index_q, index_t in bad_density[:8]]
    audit["known_density_exclusion"] = (
        {"iq_zero_based": 520, "records": len(logt), "delta_logrho": -1e-4}
        if known_exclusion else None)
    audit["coordinate_exclusion_policy"] = (
        "printed logRho discrepancies remain unchanged and are excluded by source_constraints")
    audit["grid_audit"] = (
        "all printed Q labels, T coordinates, uniform spacing and rectangularity verified; "
        "logRho mismatches audited and masked")
    return observed, logq, logt, rows, audit


def write_payload(output, metadata, logq, logt, nodes, complete):
    temporary = output.with_suffix(output.suffix + ".partial")
    with temporary.open("w") as stream:
        stream.write("FREEEOS_BIQUINTIC_V1\n")
        stream.write(f"{len(logq)} {len(logt)}\n")
        stream.write(f"{metadata['X']:.17e} {metadata['Z']:.17e} {CRAD:.17e} 1\n")
        stream.write(metadata["source_sha256"] + "\n")
        stream.write(f"{metadata['preprocessing_sha256']} 2 mean\n")
        stream.write(" ".join(f"{value:.17e}" for value in
                              (logq[0] * LN10, (logq[1] - logq[0]) * LN10,
                               logt[0] * LN10, (logt[1] - logt[0]) * LN10)) + "\n")
        for temperature_index in range(len(logt)):
            for density_index in range(len(logq)):
                accepted = bool(complete[density_index, temperature_index])
                corner = nodes[:, :, density_index, temperature_index].ravel(order="F") if accepted else np.zeros(9)
                stream.write(f"{density_index + 1} {temperature_index + 1} {int(accepted)} "
                             + " ".join(f"{value:.17e}" for value in corner) + "\n")
    os.replace(temporary, output)


def build_one(source, output):
    try:
        header, logq, logt, rows, audit = read_audited(source)
        low, valid, _ = source_constraints(rows, logq, logt, split_radiation=True, minimum_logT=3.2)
        nodes, complete, _ = construct_twists(low, valid, (logq[1] - logq[0]) * LN10,
                                              (logt[1] - logt[0]) * LN10, 2, "mean")
        audit.update({"format": "FREEEOS_BIQUINTIC_V1", "radius": 2, "route": "mean",
                      "split_radiation": True, "minimum_logT": 3.2,
                      "radiation_constant": CRAD,
                      "preprocessing_sha256": checksum(Path(__file__).resolve().parents[1] / "freeeos_potential.py"),
                      "source_valid_nodes": int(valid.sum()), "complete_jet_nodes": int(complete.sum()),
                      "complete_cells": int(valid_cells(complete).sum())})
        if not np.isfinite(nodes[:, :, complete]).all() or not complete.any():
            raise ValueError("No finite complete canonical potential jets")
        filename = source.stem + ".bq"
        write_payload(output / filename, audit, logq, logt, nodes, complete)
        audit["payload_sha256"] = checksum(output / filename)
        audit["payload_filename"] = filename
        audit["available"] = 1
    except Exception as error:
        audit = (error.audit if isinstance(error, SourceAuditError) else
                 locals().get("audit", {"source_filename": source.name,
                                        "X": source_identity(source)[0], "Z": source_identity(source)[1],
                                        "co_variant": source_identity(source)[2]}))
        audit["source_sha256"] = checksum(source)
        audit["available"] = 0
        audit["error"] = f"{type(error).__name__}: {error}"
    audit_path = output / (source.stem + ".audit.json")
    temporary = audit_path.with_suffix(audit_path.suffix + ".partial")
    temporary.write_text(json.dumps(audit, indent=2) + "\n")
    os.replace(temporary, audit_path)
    return audit


def priority(source):
    x, z, co = source_identity(source)
    plane = {0: 0, .02: 1, .04: 2, .06: 3, .08: 4}.get(z, 5)
    return (int(co), int(abs(x - .7) > TOLERANCE), plane, abs(x - .7), z, x)


def read_audit(path):
    if not path.exists() or not path.stat().st_size:
        return None
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError:
        return None


def publish(output, sources):
    records = []
    for source in sorted(sources, key=priority):
        audit_path = output / (source.stem + ".audit.json")
        audit = read_audit(audit_path)
        if audit is None:
            x, z, co = source_identity(source)
            audit = {"source_filename": source.name, "X": x, "Z": z, "co_variant": co,
                     "source_sha256": checksum(source), "available": 0, "error": "Pending producer"}
        records.append(audit)
    manifest = ["FREEEOS_BIQUINTIC_FAMILY_V1", str(len(records))]
    for audit in records:
        available = int(audit["available"])
        manifest.append(f"{audit['X']:.17e} {audit['Z']:.17e} {available} "
                        f"{audit.get('payload_filename', '-') if available else '-'} "
                        f"{audit['source_sha256']} {audit.get('payload_sha256', '-') if available else '-'}")
    temporary = output / "family.manifest.partial"
    temporary.write_text("\n".join(manifest) + "\n")
    os.replace(temporary, output / "family.manifest")
    summary = {"format": "FREEEOS_BIQUINTIC_FAMILY_V1", "entry_count": len(records),
               "available": sum(item["available"] for item in records),
               "unavailable": sum(not item["available"] for item in records), "tables": records}
    temporary = output / "family.audit.json.partial"
    temporary.write_text(json.dumps(summary, indent=2) + "\n")
    os.replace(temporary, output / "family.audit.json")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=SOURCE_DIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--retry-unavailable", action="store_true")
    arguments = parser.parse_args()
    if not 1 <= arguments.workers <= 8:
        parser.error("Use one through eight workers")
    sources = sorted(arguments.source_dir.glob("mesa-FreeEOS_*.data"), key=priority)
    if len(sources) != 89:
        parser.error(f"Expected 89 source tables, found {len(sources)}")
    arguments.output.mkdir(parents=True, exist_ok=True)
    completed = set()
    for source in sources:
        audit_path = arguments.output / (source.stem + ".audit.json")
        audit = read_audit(audit_path)
        if (audit is not None and audit.get("error") != "Pending producer" and
                not (arguments.retry_unavailable and not audit.get("available"))):
            completed.add(source)
    for source in sources:
        audit_path = arguments.output / (source.stem + ".audit.json")
        if read_audit(audit_path) is not None:
            continue
        x, z, co = source_identity(source)
        audit_path.write_text(json.dumps({"source_filename": source.name, "X": x, "Z": z,
                                          "co_variant": co, "source_sha256": checksum(source),
                                          "available": 0, "error": "Pending producer"}, indent=2) + "\n")
    publish(arguments.output, sources)
    pending = [source for source in sources if not (arguments.resume and source in completed)]
    if arguments.limit is not None:
        pending = pending[:arguments.limit]
    with ProcessPoolExecutor(max_workers=arguments.workers) as pool:
        futures = {pool.submit(build_one, source, arguments.output): source for source in pending}
        for future in as_completed(futures):
            audit = future.result()
            summary = publish(arguments.output, sources)
            print(f"{audit['source_filename']}: available={audit['available']} "
                  f"nodes={audit.get('complete_jet_nodes', 0)} "
                  f"error={audit.get('error', '-')} total={summary['available']}/{len(sources)}", flush=True)


if __name__ == "__main__":
    main()
