"""Classify identical edge-control holdouts/faces with one actual-MESA batch."""

import os

os.environ["OPENBLAS_NUM_THREADS"] = "8"
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

import io
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

from edge_control import (
    JACOBIANS, VALUES, build_jets, comparison, distribution, face_comparison,
    positivity, state_at,
)
from active_faces import classifications
from freeeos_potential import LN10, checksum, construct_twists, query_cells, read_table, source_constraints, valid_cells
from mesa_probe import run_probe
from run_audit import own_use
from validate_table import material_nodes, table_state


CORE = Path(__file__).resolve().parents[1]
WORKSPACE = CORE.parents[3]
OUTPUT = WORKSPACE / "freeeos_biquintic_progress/experiments/020_edge_control"
CATEGORIES = ("all_successful", "active", "pure_FreeEOS", "blend", "inactive", "failed", "own_use", "active_own_use")
TABLE_FIELDS = ("Pgas", "Cv", "Cp", "chiRho", "chiT", "gamma1", "gamma3", "grad_ad")


def probe_child():
    run_probe.REPORT = OUTPUT / "active_quality_probe_runtime"
    sys.argv = [str(CORE / "mesa_probe/run_probe.py")]
    run_probe.main()


def holdout_population(rows, valid, logq, logt, expected):
    coarse_q, coarse_t, coarse_rows = logq[::2], logt[::2], rows[::2, ::2]
    coarse_low, coarse_valid, _ = source_constraints(coarse_rows, coarse_q, coarse_t, minimum_logT=3.2)
    du, dtau = (coarse_q[1] - coarse_q[0]) * LN10, (coarse_t[1] - coarse_t[0]) * LN10
    nodes, complete, _ = build_jets(coarse_low, coarse_valid, du, dtau)
    _, core_complete, _ = construct_twists(coarse_low, coarse_valid, du, dtau)
    density, temperature = np.nonzero(valid & material_nodes(rows, logq, logt)
        & (logt[None, :] >= 3.2 - 1e-12) & (logt[None, :] <= 3.6 + 1e-12))
    withheld = (density % 2 == 1) & (temperature % 2 == 1)
    density, temperature = density[withheld], temperature[withheld]
    sample_q, sample_t = logq[density], logt[temperature]
    cell_q, cell_t, inside = query_cells(coarse_q, coarse_t, sample_q, sample_t)
    accepted = inside & valid_cells(complete)[cell_q, cell_t]
    overlap = inside & valid_cells(core_complete)[cell_q, cell_t]
    assert len(sample_q) == expected["source_native_material_withheld"]
    assert int(accepted.sum()) == expected["accepted_control"]
    assert int(overlap.sum()) == expected["accepted_core"]
    sample_q, sample_t = sample_q[accepted], sample_t[accepted]
    density, temperature, overlap = density[accepted], temperature[accepted], overlap[accepted]
    state = state_at(nodes, complete, coarse_q, coarse_t, sample_q, sample_t, False)
    reference = table_state(rows[density, temperature])
    for name, selected in (("new_coverage", ~overlap), ("overlap", overlap)):
        candidate = {field: values[selected] for field, values in state.items()}
        target = {field: values[selected] for field, values in reference.items()}
        actual = comparison(candidate, target, sample_q[selected], sample_t[selected], tuple(reference))
        if actual != expected["groups"][name]["archived_relative_departure"]:
            raise AssertionError("Holdout sampling or values differ from experiment020")
    return {"q": sample_q, "t": sample_t, "overlap": overlap, "state": state, "reference": reference,
            "native_density_indices": density, "native_temperature_indices": temperature}


def face_populations(nodes, complete, core_complete, logq, logt, rows, expected):
    cells = valid_cells(complete) & valid_cells(material_nodes(rows, logq, logt))
    cells &= logt[None, :-1] < 3.36 + 1e-12
    generator = np.random.default_rng(20261028020)
    populations = {}
    for axis, direction in ((0, "u"), (1, "tau")):
        adjacent = cells[:-1, :] & cells[1:, :] if axis == 0 else cells[:, :-1] & cells[:, 1:]
        density, temperature = np.nonzero(adjacent)
        count = expected[direction]["count"]
        chosen = generator.choice(len(density), min(count, len(density)), replace=False)
        density, temperature = density[chosen], temperature[chosen]
        fraction = generator.uniform(.1, .9, len(chosen))
        if axis == 0:
            sample_q = logq[density + 1]
            sample_t = logt[temperature] + fraction * (logt[1] - logt[0])
            right_cells = density + 1, temperature
        else:
            sample_q = logq[density] + fraction * (logq[1] - logq[0])
            sample_t = logt[temperature + 1]
            right_cells = density, temperature + 1
        left_cells = density, temperature
        left = state_at(nodes, complete, logq, logt, sample_q, sample_t, False, left_cells)
        right = state_at(nodes, complete, logq, logt, sample_q, sample_t, False, right_cells)
        overlap = valid_cells(core_complete)[left_cells] & valid_cells(core_complete)[right_cells]
        if face_comparison(left, right, sample_q, sample_t) != expected[direction]["scaled_jumps"]:
            raise AssertionError("Exact face sampling or states differ from experiment020")
        assert int(overlap.sum()) == expected[direction]["overlap_core_faces"]
        populations[direction] = {"q": sample_q, "t": sample_t, "overlap": overlap, "state": left, "reference": right,
            "left_cells": np.column_stack(left_cells), "right_cells": np.column_stack(right_cells)}
    return populations


def signed_comparison(candidate, reference, sample_q, sample_t, fractions, faces=False):
    fields = VALUES + JACOBIANS if faces else TABLE_FIELDS
    metrics = {}
    for field in fields:
        scale = abs(reference[field])
        if faces:
            scale = np.maximum(scale, abs(candidate[field]))
            if field in JACOBIANS:
                scale = np.maximum.reduce([scale, abs(candidate[field[:-2]]), abs(reference[field[:-2]])])
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            signed = (candidate[field] - reference[field]) / np.maximum(scale, np.finfo(float).tiny)
        metric = distribution(signed)
        finite = np.isfinite(signed)
        if finite.any():
            values = signed[finite]
            metric["signed"] = dict(zip(("p01", "median", "p99", "minimum", "maximum"),
                map(float, np.quantile(values, (.01, .5, .99, 0, 1)))))
            index = int(np.argmax(np.where(finite, abs(signed), -1)))
            metric["worst"] = {"logQ": float(sample_q[index]), "logT": float(sample_t[index]),
                "signed_scaled_difference": float(signed[index]),
                "fraction_FreeEOS": float(fractions[index]) if np.isfinite(fractions[index]) else None,
                "candidate_or_left": float(candidate[field][index]), "printed_or_right": float(reference[field][index])}
        metrics[field] = metric
    return metrics


def classified_groups(population, probe, faces=False):
    own = own_use(population["q"], population["t"])
    masks = classifications(probe[:, 2], probe[:, 5])
    masks["own_use"] = own & masks["all_successful"]
    masks["active_own_use"] = own & masks["active"]
    groups = {}
    for name, selected in (("new_coverage", ~population["overlap"]), ("overlap", population["overlap"])):
        categories = {}
        for category in CATEGORIES:
            chosen = selected & masks[category]
            state = {field: values[chosen] for field, values in population["state"].items()}
            target = {field: values[chosen] for field, values in population["reference"].items()}
            sample_q, sample_t = population["q"][chosen], population["t"][chosen]
            data = {"count": int(chosen.sum()), "own_use_count": int((chosen & own).sum()),
                "metrics": signed_comparison(state, target, sample_q, sample_t, probe[chosen, 5], faces),
                "candidate_positivity": positivity(state, sample_q, sample_t)}
            if chosen.any():
                data["coordinates"] = {"logQ": [float(sample_q.min()), float(sample_q.max())],
                    "logT": [float(sample_t.min()), float(sample_t.max())]}
                if np.isfinite(probe[chosen, 5]).any():
                    data["fraction_FreeEOS_range"] = [float(np.nanmin(probe[chosen, 5])), float(np.nanmax(probe[chosen, 5]))]
            if faces:
                data["candidate_right_positivity"] = positivity(target, sample_q, sample_t)
            categories[category] = data
        groups[name] = {"count": int(selected.sum()), "categories": categories}
    return groups


def append_summary(result):
    lines = ["", "## Actual-MESA active-quality appendix", "", "One batch classifies the IDENTICAL 1,399 accepted stride2 material holdouts and 2,000 exact material faces from the main experiment. All original unclassified errors/jumps are reproduced exactly before probing. These material populations are not prefiltered by `own_use`; that predicate is evaluated separately here.",
        "", "Default-MESA scaled-solar X=.70/Y=.28/Z=.02 fractions classify the coordinates; the candidate is not installed in MESA. Printed rows are fidelity references, not independent accuracy. New/overlap populations still differ; no stencil-causation or runtime-backend claim.",
        "", "### Holdout classification", "", "| Population / actual MESA class | Rows | Also own_use | Cv p99 / max | Cp p99 / max | chiRho p99 / max | chiT p99 / max |", "|---|---:|---:|---|---|---|---|"]
    for name, group in result["holdouts"].items():
        for category in ("all_successful", "active", "pure_FreeEOS", "blend", "inactive"):
            data = group["categories"][category]
            values = []
            for field in ("Cv", "Cp", "chiRho", "chiT"):
                metric = data["metrics"][field]
                values.append(f"{metric['p99']:.6g} / {metric['max']:.6g}" if data["count"] else "none")
            lines.append(f"| {name} / {category} | {data['count']} | {data['own_use_count']} | " + " | ".join(values) + " |")
    worst = result["holdouts"]["new_coverage"]["categories"]["all_successful"]["metrics"]["Cv"]["worst"]
    lines.extend(["", f"The original Cv maximum at logQ={worst['logQ']:.12g}, logT={worst['logT']:.12g} has actual fraction_FreeEOS={worst['fraction_FreeEOS']:.12g}; signed relative departure={worst['signed_scaled_difference']:.6g}. Do not attribute an inactive-state tail to active FreeEOS.",
        "", "### Exact response jumps by actual MESA fraction", "", "| Orientation / population / class | Faces | Also own_use | Cv_t median / p99 / max | chiT_t median / p99 / max | Cp_t median / p99 / max |", "|---|---:|---:|---|---|---|"])
    for direction, populations in result["faces"].items():
        for name, group in populations.items():
            for category in ("active", "pure_FreeEOS", "blend", "inactive"):
                data = group["categories"][category]
                values = []
                for field in ("Cv_t", "chiT_t", "Cp_t"):
                    metric = data["metrics"][field]
                    values.append(" / ".join(f"{metric[key]:.6g}" for key in ("median", "p99", "max")) if data["count"] else "none")
                lines.append(f"| {direction} / {name} / {category} | {data['count']} | {data['own_use_count']} | " + " | ".join(values) + " |")
    failed = sum(group["categories"]["failed"]["count"] for group in result["holdouts"].values())
    failed += sum(group["categories"]["failed"]["count"] for populations in result["faces"].values() for group in populations.values())
    violations = sum(metric["count"] for group in result["holdouts"].values() for metric in group["categories"]["active"]["candidate_positivity"].values())
    violations += sum(metric["count"] for populations in result["faces"].values() for group in populations.values()
        for side in ("candidate_positivity", "candidate_right_positivity") for metric in group["categories"]["active"][side].values())
    lines.extend(["", f"Probe failures={failed}; sampled active candidate positivity violations={violations} (summed across fields/sides, not unique states). E/S sign checks depend on the printed convention, not gauge-invariant stability.",
        "", "Table departures use abs(printed reference) denominators. Exact Jacobian jumps use max(abs(left/right derivative), abs(left/right response)); x/t are physical natural-log rho/T. Tables give magnitude quantiles; active_quality.json also preserves signed quantiles, worst coordinates and actual fractions. Gamma3 uses the printed Gamma3, not Gamma3-minus1.",
        "", "Source constraints/donors remain clean logT>=3.2; stride2 starts at3.22, with no extrapolation to3.2 and no bug repair below3.2. C2, not C3, remains the continuity class. Fractions are measured for the baseline MESA blend, not recomputed for an installed candidate; no source-accuracy, final-blend-error or solver-safety assertion."])
    marker = "\n## Actual-MESA active-quality appendix"
    path = OUTPUT / "SUMMARY.md"
    original = path.read_text().split(marker)[0].rstrip()
    path.write_text(original + "\n" + "\n".join(lines) + "\n")


def main():
    if sys.argv[1:] == ["--probe-child"]:
        probe_child()
        return
    reuse_probe = sys.argv[1:] == ["--reuse-probe"]
    if sys.argv[1:] and not reuse_probe:
        raise SystemExit("Usage: active_quality.py [--reuse-probe]")
    started = time.perf_counter()
    baseline_path = OUTPUT / "results.json"
    baseline = json.loads(baseline_path.read_text())
    if baseline["method"]["split_radiation"]:
        raise ValueError("This appendix reproduces the primary total-potential run only")
    for path, expected in baseline["provenance"]["input_sha256"].items():
        if checksum(path) != expected and Path(path).resolve() != (CORE / "freeeos_potential.py").resolve():
            raise ValueError(f"Main experiment snapshot changed: {path}")
    inputs = (baseline_path, CORE / "edge_control/edge_control.py", CORE / "freeeos_potential.py", CORE / "validate_table.py",
        CORE / "active_faces.py", CORE / "run_audit.py", CORE / "mesa_probe/run_probe.py", CORE / "mesa_probe/probe.f90",
        CORE / "mesa_probe/build/mesa_probe", CORE / "mesa_probe/inlist_probe", CORE.parent / "mass_frac_solar.txt", Path(__file__))
    fingerprints = {str(path): checksum(path) for path in inputs}
    table_path = WORKSPACE / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
    fingerprints[str(table_path)] = checksum(table_path)
    _, logq, logt, rows = read_table(table_path)
    low, valid, _ = source_constraints(rows, logq, logt, minimum_logT=3.2)
    du, dtau = (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10
    nodes, complete, _ = build_jets(low, valid, du, dtau)
    _, core_complete, _ = construct_twists(low, valid, du, dtau)
    with np.load(OUTPUT / "stencils.npz") as original:
        np.testing.assert_array_equal(valid, original["source_valid"])
        np.testing.assert_array_equal(complete, original["control_complete"])
        np.testing.assert_array_equal(core_complete, original["core_complete"])
    holdouts = holdout_population(rows, valid, logq, logt, baseline["stride2_lowedge_holdouts"])
    faces = face_populations(nodes, complete, core_complete, logq, logt, rows, baseline["lowedge_material_faces"])
    populations = {"holdouts": holdouts, **faces}
    query_q = np.concatenate([data["q"] for data in populations.values()])
    query_t = np.concatenate([data["t"] for data in populations.values()])
    if (query_t < 3.2).any():
        raise AssertionError("Query batch crossed the clean source floor")
    query_text = "".join(f"{density:.17g} {temperature:.17g}\n" for density, temperature in zip(query_q + 2 * query_t - 12, query_t))
    query_path = OUTPUT / "active_quality_queries.txt"
    command = [sys.executable, "-B", str(Path(__file__).resolve()), "--probe-child"]
    probe_path = OUTPUT / "active_quality_probe.tsv"
    execution_driver_sha256 = fingerprints[str(Path(__file__))]
    if reuse_probe:
        previous = json.loads((OUTPUT / "active_quality.json").read_text())
        if query_path.read_text() != query_text:
            raise ValueError("Existing probe coordinates differ from this reconstruction")
        for path, expected in previous["provenance"]["input_sha256"].items():
            if Path(path).resolve() not in (Path(__file__).resolve(), (CORE / "freeeos_potential.py").resolve()) and checksum(path) != expected:
                raise ValueError(f"Original probe snapshot changed: {path}")
        if "probe_output_sha256" in previous["provenance"] and checksum(probe_path) != previous["provenance"]["probe_output_sha256"]:
            raise ValueError("Recorded probe output changed")
        execution_driver_sha256 = previous["provenance"].get("probe_driver_sha256_at_execution",
            previous["provenance"]["input_sha256"][str(Path(__file__))])
    else:
        query_path.write_text(query_text)
        completed = subprocess.run(command, input=query_text, capture_output=True, text=True, timeout=600)
        probe_path.write_text(completed.stdout)
        (OUTPUT / "active_quality_probe.stderr").write_text(completed.stderr)
        if completed.returncode:
            raise RuntimeError(f"Probe child failed: {completed.returncode}; see experiment020 stderr")
    probe = np.loadtxt(io.StringIO(probe_path.read_text()), skiprows=1, ndmin=2)
    if probe.shape != (len(query_q), 14):
        raise ValueError(f"Unexpected probe shape: {probe.shape}")
    np.testing.assert_allclose(probe[:, :2], np.column_stack((query_q + 2 * query_t - 12, query_t)), rtol=0, atol=1e-12)
    metadata_path = OUTPUT / "active_quality_probe_runtime/probe_metadata.json"
    metadata = json.loads(metadata_path.read_text())
    if not metadata["installed_cache_inventory_unchanged"] or metadata["private_replaced_cache_links"]:
        raise RuntimeError("Probe cache safety failed")
    result = {"purpose": "One actual-MESA fraction batch for identical experiment020 accepted material holdouts and exact faces",
        "provenance": {"input_sha256": fingerprints, "query_sha256": checksum(query_path), "command": command,
            "original_main_input_sha256": baseline["provenance"]["input_sha256"],
            "current_explicit_3.2_masks_and_unclassified_metrics_match_original_exactly": True,
            "probe_output_sha256": checksum(probe_path), "reanalysis_without_new_probe": reuse_probe,
            "probe_driver_sha256_at_execution": execution_driver_sha256,
            "python": sys.version, "numpy": np.__version__, "probe_metadata": metadata},
        "scope": {"queries": len(query_q), "same_unclassified_statistics_exactly_reproduced": True,
            "material_not_prefiltered_own_use": True, "own_use_definition": "existing run_audit.own_use evaluated separately; not an actual-MESA fraction proxy",
            "source_floor": 3.2, "minimum_query_logT": float(query_t.min()), "coarse_donor_floor": 3.22,
            "classification": "active fraction>1e-12; pure>=1-1e-12; blend strictly between; ierr0 and finite fraction required",
            "not_measured": "independent source accuracy, final blended candidate errors, stencil causality across unlike populations, installed candidate fractions, solver/backend safety",
            "normalization": "holdout candidate-minus-printed / abs(printed); faces left-minus-right / max(abs(left/right)), additionally abs(left/right response) for Jacobians; physical natural-log rho/T partials"},
        "holdouts": {}, "faces": {}}
    offset = 0
    archives = {}
    for name, population in populations.items():
        count = len(population["q"])
        local_probe = probe[offset:offset + count]
        offset += count
        groups = classified_groups(population, local_probe, name != "holdouts")
        if name == "holdouts":
            result["holdouts"] = groups
        else:
            result["faces"][name] = groups
        for field in ("q", "t", "overlap", "native_density_indices", "native_temperature_indices", "left_cells", "right_cells"):
            if field in population:
                archives[name + "_" + field] = population[field]
        archives[name + "_probe"] = local_probe
        archives[name + "_own_use"] = own_use(population["q"], population["t"])
    result["elapsed_seconds"] = time.perf_counter() - started
    for path, expected in fingerprints.items():
        if checksum(path) != expected:
            raise RuntimeError(f"Input snapshot changed during run: {path}")
    (OUTPUT / "active_quality.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    np.savez_compressed(OUTPUT / "active_quality_queries.npz", **archives)
    append_summary(result)
    print(json.dumps({"queries": len(query_q), "holdouts": {name: {category: data["count"] for category, data in group["categories"].items()}
        for name, group in result["holdouts"].items()}, "faces": {direction: {name: {category: data["count"] for category, data in group["categories"].items()}
        for name, group in groups.items()} for direction, groups in result["faces"].items()}, "elapsed_seconds": result["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
