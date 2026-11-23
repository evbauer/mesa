"""Generate lossless FreeEOS mixed-derivative supplements from distributed sources.

Only source tables and the MESA source checkout are required. An optional
authenticated V1 family supplies a numerical oracle, never production input.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import h5py
import numpy as np
from family_cache.produce import read_audited
from freeeos_potential import CRAD, LN10, checksum, construct_twists, source_constraints, valid_cells

POLICY_SHA = checksum(Path(__file__).resolve().parents[1] / "freeeos_potential.py")
MIXED = {"f_uu_tau": (2, 1), "f_u_tau_tau": (1, 2), "f_uu_tau_tau": (2, 2)}


def digest_bytes(hex_digest):
    return np.asarray(list(bytes.fromhex(hex_digest)), dtype=np.int32)


def read_reference(path, expected_sha, source_sha):
    if checksum(path) != expected_sha:
        raise ValueError(f"V1 payload SHA mismatch: {path}")
    with path.open() as stream:
        header = [next(stream).strip() for _ in range(6)]
    if header[0] != "FREEEOS_BIQUINTIC_V1" or header[3] != source_sha:
        raise ValueError("V1 source identity mismatch")
    if header[4].split() != [POLICY_SHA, "2", "mean"]:
        raise ValueError("V1 preprocessing policy mismatch")
    nu, nt = map(int, header[1].split())
    records = np.loadtxt(path, skiprows=6).reshape(nt, nu, 12)
    if not np.array_equal(records[:, :, 0], np.broadcast_to(np.arange(1, nu+1), (nt, nu))):
        raise ValueError("V1 u record order mismatch")
    if not np.array_equal(records[:, :, 1], np.broadcast_to(np.arange(1, nt+1)[:, None], (nt, nu))):
        raise ValueError("V1 tau record order mismatch")
    return records


def build_one(source, output, reference):
    source, output = Path(source), Path(output)
    header, logq, logt, rows, audit = read_audited(source)
    low, valid, _ = source_constraints(rows, logq, logt, split_radiation=True, minimum_logT=3.2)
    du, dtau = (logq[1]-logq[0])*LN10, (logt[1]-logt[0])*LN10
    nodes, complete, _ = construct_twists(low, valid, du, dtau, 2, "mean")
    if not complete.any() or not np.isfinite(nodes[:, :, complete]).all():
        raise ValueError("No finite accepted source jet")
    if reference is not None:
        records = read_reference(Path(reference[0]), reference[1], audit["source_sha256"])
        if not np.array_equal(records[:, :, 2].T, complete):
            raise ValueError("V1 accepted mask differs from source generation")
        oracle = records[:, :, 3:].transpose(2, 1, 0).reshape(3, 3, len(logq), len(logt), order="F")
        # Each text float round-trips exactly; this tests the complete construction.
        if not np.array_equal(oracle[:, :, complete], nodes[:, :, complete]):
            raise ValueError("Source-generated nine constraints differ from authenticated V1")
        audit["reference_payload_sha256"] = reference[1]
        audit["reference_all_nine_constraints_exact"] = True
    attrs = {"format_version": 1, "source_version": 51, "X": float(header[1]), "Z": float(header[2]),
             "family_id": int(audit["co_variant"]), "nu": len(logq), "ntau": len(logt),
             "u_min": float(logq[0]*LN10), "du": float(du), "tau_min": float(logt[0]*LN10),
             "dtau": float(dtau), "radiation_split": 1, "radiation_constant": float(CRAD),
             "stencil_radius": 2, "twist_route": 1, "minimum_log10_T": 3.2, "preprocessing_version": 1}
    destination = output / (source.stem + ".bq.h5")
    temporary = destination.with_suffix(".h5.partial")
    with h5py.File(temporary, "w", libver="earliest") as h5:
        for name, value in attrs.items():
            h5.attrs[name] = np.int32(value) if isinstance(value, int) else np.float64(value)
        for name, orders in MIXED.items():
            values = np.where(complete, nodes[orders], 0.0).T
            h5.create_dataset(name, data=values, dtype="<f8", compression="gzip", compression_opts=6,
                              shuffle=True, fletcher32=True, track_times=False)
        h5.create_dataset("valid", data=complete.T.astype(np.int8), dtype="i1", compression="gzip",
                          compression_opts=6, track_times=False)
        h5.create_dataset("source_sha256", data=digest_bytes(audit["source_sha256"]), track_times=False)
        h5.create_dataset("preprocessing_sha256", data=digest_bytes(POLICY_SHA), track_times=False)
    temporary.chmod(0o644)
    os.replace(temporary, destination)
    audit.update({"format": "FREEEOS_BIQUINTIC_SUPPLEMENTS_V1", "generation": "distributed source preprocessing",
                  "supplement_filename": destination.name, "supplement_sha256": checksum(destination),
                  "supplement_bytes": destination.stat().st_size, "preprocessing_sha256": POLICY_SHA,
                  "source_valid_nodes": int(valid.sum()), "complete_jet_nodes": int(complete.sum()),
                  "complete_cells": int(valid_cells(complete).sum()), "metadata": attrs})
    (output / (source.stem + ".bq.audit.json")).write_text(json.dumps(audit, indent=2, sort_keys=True)+"\n")
    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--reference-dir", type=Path, help="Optional authenticated V1 family")
    parser.add_argument("--workers", type=int, default=4, choices=range(1,9))
    parser.add_argument("--limit", type=int, help="Partial developer fixture; never a distribution bundle")
    args = parser.parse_args()
    sources = sorted(args.source_dir.glob("mesa-FreeEOS_*.data"))
    if len(sources) != 89:
        parser.error(f"Expected all 89 distributed sources, found {len(sources)}")
    if args.limit:
        sources = sources[:args.limit]
    references = {}
    if args.reference_dir:
        lines = (args.reference_dir / "family.manifest").read_text().splitlines()
        if lines[0] != "FREEEOS_BIQUINTIC_FAMILY_V1" or int(lines[1]) != len(lines)-2:
            parser.error("Malformed V1 reference manifest")
        for line in lines[2:]:
            x, z, available, filename, source_sha, payload_sha = line.split()
            if available != "1":
                parser.error("Reference family is incomplete")
            references[filename] = (str(args.reference_dir / filename), payload_sha, source_sha)
    args.output.mkdir(parents=True, exist_ok=True)
    records = []
    with ProcessPoolExecutor(args.workers) as pool:
        jobs = {}
        for source in sources:
            reference = references.get(source.stem+".bq") if args.reference_dir else None
            if args.reference_dir and reference is None:
                parser.error(f"Missing reference {source.name}")
            jobs[pool.submit(build_one, str(source), str(args.output), reference)] = source
        for job in as_completed(jobs):
            record = job.result()
            records.append(record)
            print(f"{len(records)}/{len(sources)} {record['supplement_filename']} {record['supplement_bytes']} bytes", flush=True)
    records.sort(key=lambda item: item["source_filename"])
    manifest = ["FREEEOS_BIQUINTIC_SUPPLEMENTS_V1", str(len(records))]
    for item in records:
        manifest.append(f"{item['X']:.17e} {item['Z']:.17e} 1 {item['supplement_filename']} "
                        f"{item['source_sha256']} {item['supplement_sha256']}")
    (args.output / "supplements.manifest").write_text("\n".join(manifest)+"\n")
    summary = {"format": manifest[0], "entry_count": len(records), "source_count": len(records),
               "co_count": sum(item["co_variant"] for item in records),
               "supplement_bytes": sum(item["supplement_bytes"] for item in records), "tables": records}
    (args.output / "supplements.audit.json").write_text(json.dumps(summary, indent=2, sort_keys=True)+"\n")


if __name__ == "__main__":
    main()
