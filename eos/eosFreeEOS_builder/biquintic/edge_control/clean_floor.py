"""Explicit clean-floor native shifted-stencil CONTROL, not runtime integration."""

import os

os.environ["OPENBLAS_NUM_THREADS"] = "8"
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

import argparse
import io
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

from edge_control import JACOBIANS, VALUES, build_jets, comparison, face_comparison, state_at
from active_quality import classified_groups
from freeeos_potential import LN10, checksum, construct_twists, query_cells, read_table, source_constraints, valid_cells
from mesa_probe import run_probe
from validate_table import material_nodes, table_state


CORE = Path(__file__).resolve().parents[1]
WORKSPACE = CORE.parents[3]


def output_directory(source_floor, split):
    output = WORKSPACE / "freeeos_biquintic_progress/experiments/027_clean_floor"
    if source_floor != 3.1:
        output /= f"floor_{source_floor:g}"
    if split:
        output /= "split_radiation"
    return output


def support_report(complete, metadata, logq, logt, source_floor):
    density, temperature = np.indices(complete.shape)
    bounds = metadata["support"]
    minimum = float(logt[bounds[..., 2][complete]].min())
    reach_q = np.maximum(density - bounds[..., 0], bounds[..., 1] - density)[complete]
    reach_t = np.maximum(temperature - bounds[..., 2], bounds[..., 3] - temperature)[complete]
    if minimum < source_floor - 1e-12 or source_floor < 3.1 or not np.isfinite(source_floor):
        raise AssertionError("Donor support crossed the tested clean source floor")
    if np.max(reach_q) > 8 or np.max(reach_t) > 8:
        raise AssertionError("Two-stage native support exceeded eight grid spacings")
    return {"tested_source_floor": source_floor, "minimum_donor_logT": minimum,
        "minimum_complete_node_logT": float(logt[np.nonzero(complete)[1]].min()),
        "maximum_density_grid_reach": int(reach_q.max()), "maximum_temperature_grid_reach": int(reach_t.max()),
        "maximum_logQ_reach": float(reach_q.max() * (logq[1] - logq[0])),
        "maximum_logT_reach": float(reach_t.max() * (logt[1] - logt[0])),
        "support_interpretation": "conservative unions of numerical and eligibility donors; not every enclosed point is used",
        "shifted_window_counts": {name: int(np.sum((data["position"] >= 0) & (data["position"] != 2)))
            for name, data in metadata["stencils"].items()}}


def accepted_queries(complete, logq, logt, sample_q, sample_t):
    density, temperature, inside = query_cells(logq, logt, sample_q, sample_t)
    return inside & valid_cells(complete)[density, temperature]


def build_control(rows, logq, logt, source_floor, split):
    low, valid, _ = source_constraints(rows, logq, logt, split_radiation=split, minimum_logT=source_floor)
    du, dtau = (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10
    nodes, complete, metadata = build_jets(low, valid, du, dtau)
    core_nodes, core_complete, _ = construct_twists(low, valid, du, dtau)
    report = support_report(complete, metadata, logq, logt, source_floor)
    return {"low": low, "valid": valid, "nodes": nodes, "complete": complete, "metadata": metadata,
        "centered_nodes": core_nodes, "centered_complete": core_complete, "support": report}


def nodal_population(rows, logq, logt, candidate, reference_control, source_floor, split):
    density, temperature = np.nonzero(candidate["valid"] & material_nodes(rows, logq, logt)
        & (logt[None, :] >= source_floor - 1e-12) & (logt[None, :] <= 3.2 + 1e-12))
    sample_q, sample_t = logq[density], logt[temperature]
    selected = accepted_queries(candidate["complete"], logq, logt, sample_q, sample_t)
    excluded = int((~selected).sum())
    density, temperature = density[selected], temperature[selected]
    sample_q, sample_t = sample_q[selected], sample_t[selected]
    return {"q": sample_q, "t": sample_t,
        "overlap": accepted_queries(reference_control["complete"], logq, logt, sample_q, sample_t),
        "state": state_at(candidate["nodes"], candidate["complete"], logq, logt, sample_q, sample_t, split),
        "reference": table_state(rows[density, temperature]), "density_index": density, "temperature_index": temperature,
        "excluded_native_material_nodes_without_complete_query_cell": excluded, "tested_source_floor": source_floor}


def holdout_population(rows, logq, logt, source_floor, split):
    coarse_q, coarse_t, coarse_rows = logq[::2], logt[::2], rows[::2, ::2]
    candidate = build_control(coarse_rows, coarse_q, coarse_t, source_floor, split)
    reference_control = build_control(coarse_rows, coarse_q, coarse_t, 3.2, split)
    _, full_valid, _ = source_constraints(rows, logq, logt, split_radiation=split, minimum_logT=source_floor)
    density, temperature = np.nonzero(full_valid & material_nodes(rows, logq, logt)
        & (logt[None, :] >= source_floor - 1e-12) & (logt[None, :] <= 3.4 + 1e-12))
    withheld = (density % 2 == 1) & (temperature % 2 == 1)
    density, temperature = density[withheld], temperature[withheld]
    sample_q, sample_t = logq[density], logt[temperature]
    selected = accepted_queries(candidate["complete"], coarse_q, coarse_t, sample_q, sample_t)
    total = len(sample_q)
    density, temperature, sample_q, sample_t = density[selected], temperature[selected], sample_q[selected], sample_t[selected]
    overlap = accepted_queries(reference_control["complete"], coarse_q, coarse_t, sample_q, sample_t)
    centered = accepted_queries(candidate["centered_complete"], coarse_q, coarse_t, sample_q, sample_t)
    state = state_at(candidate["nodes"], candidate["complete"], coarse_q, coarse_t, sample_q, sample_t, split)
    reference = table_state(rows[density, temperature])
    comparison_control = {}
    for label, chosen, control_nodes, control_complete in (
        ("same_floor_centered", centered, candidate["centered_nodes"], candidate["centered_complete"]),
        ("floor3.2_shifted", overlap, reference_control["nodes"], reference_control["complete"])):
        control_state = state_at(control_nodes, control_complete, coarse_q, coarse_t, sample_q[chosen], sample_t[chosen], split)
        comparison_control[label] = {"count": int(chosen.sum()), "tested_source_floor": source_floor if label == "same_floor_centered" else 3.2,
            "printed_relative_departures": comparison(control_state, {field: values[chosen] for field, values in reference.items()},
                sample_q[chosen], sample_t[chosen], tuple(reference))}
    return {"q": sample_q, "t": sample_t, "overlap": overlap, "state": state, "reference": reference,
        "density_index": density, "temperature_index": temperature, "same_floor_centered_accepted": centered,
        "total_withheld_material_nodes": total, "accepted_withheld_material_nodes": len(sample_q),
        "coarse_support": candidate["support"], "tested_source_floor": source_floor, "comparison_controls": comparison_control}


def face_populations(rows, logq, logt, candidate, reference_control, source_floor, split, maximum):
    cells = valid_cells(candidate["complete"]) & valid_cells(material_nodes(rows, logq, logt))
    cells &= (logt[None, :-1] >= source_floor - 1e-12) & (logt[None, :-1] <= 3.22 + 1e-12)
    generator = np.random.default_rng(20261028027)
    populations = {}
    for axis, direction in ((0, "u"), (1, "tau")):
        adjacent = cells[:-1, :] & cells[1:, :] if axis == 0 else cells[:, :-1] & cells[:, 1:]
        density, temperature = np.nonzero(adjacent)
        chosen = generator.choice(len(density), min(maximum, len(density)), replace=False)
        density, temperature = density[chosen], temperature[chosen]
        fraction = generator.uniform(.1, .9, len(chosen))
        if axis == 0:
            sample_q, sample_t = logq[density + 1], logt[temperature] + fraction * (logt[1] - logt[0])
            right_cells = density + 1, temperature
        else:
            sample_q, sample_t = logq[density] + fraction * (logq[1] - logq[0]), logt[temperature + 1]
            right_cells = density, temperature + 1
        left_cells = density, temperature
        overlap = valid_cells(reference_control["complete"])[left_cells] & valid_cells(reference_control["complete"])[right_cells]
        populations[direction] = {"q": sample_q, "t": sample_t, "overlap": overlap,
            "state": state_at(candidate["nodes"], candidate["complete"], logq, logt, sample_q, sample_t, split, left_cells),
            "reference": state_at(candidate["nodes"], candidate["complete"], logq, logt, sample_q, sample_t, split, right_cells),
            "left_cells": np.column_stack(left_cells), "right_cells": np.column_stack(right_cells), "tested_source_floor": source_floor}
    return populations


def summary(result):
    lines = ["# Explicit clean-floor native CONTROL", "", f"Tested source floor={result['method']['tested_source_floor']}; total/split reference={result['method']['split_radiation']}. No donors below the tested floor, no repair below3.1, no default/runtime/backend change.",
        "", "## Coverage", "", "| Archived population | Core3.2 centered | Control3.2 shifted | Tested-floor centered | Tested-floor shifted |", "|---|---:|---:|---:|---:|"]
    for band, counts in result["coverage"]["bands"].items():
        lines.append(f"| {band} | {counts['core3.2_centered']} | {counts['control3.2_shifted']} | {counts['tested_floor_centered']} | {counts['tested_floor_shifted']} |")
    lines.extend(["", f"Additional accepted archived queries over shifted3.2={result['coverage']['additional_over_shifted3.2']}; default centered complete jets match bitwise={result['coverage']['default_centered_jets_unchanged']}. Existing shifted3.2 edge jets can change because earlier clean donors now permit centered windows; do not assume shifted-overlap invariance.",
        "", f"Minimum complete node and donor logT={result['support']['minimum_complete_node_logT']} / {result['support']['minimum_donor_logT']}; maximum two-stage donor reach q/T={result['support']['maximum_density_grid_reach']} / {result['support']['maximum_temperature_grid_reach']} native spacings.",
        "", "## Actual-MESA classification and quality", "", "One batch evaluates baseline-MESA fractions at all additional accepted archived queries, native material closure nodes, coarse withheld material nodes and exact shared material faces. `own_use` is a separate predicate, not an actual-fraction proxy. New/overlap refer to tested-floor shifted versus shifted3.2; their populations differ, so no causality claim.",
        "", "| Population / class | Count | Also own_use | Cv p99 / max | Cp p99 / max | chiRho p99 / max | chiT p99 / max |", "|---|---:|---:|---|---|---|---|"])
    for label in ("native_closure", "coarse_holdouts"):
        for name, group in result["classified"][label].items():
            for category in ("all_successful", "active", "pure_FreeEOS", "blend", "inactive"):
                data = group["categories"][category]
                metrics = [f"{data['metrics'][field]['p99']:.6g} / {data['metrics'][field]['max']:.6g}" if data["count"] else "none"
                    for field in ("Cv", "Cp", "chiRho", "chiT")]
                lines.append(f"| {label}/{name}/{category} | {data['count']} | {data['own_use_count']} | " + " | ".join(metrics) + " |")
    lines.extend(["", "### Recovered archived queries", ""])
    for category, data in result["classified"]["recovered_archive"]["new_coverage"]["categories"].items():
        if category not in ("all_successful", "active", "pure_FreeEOS", "blend", "inactive", "failed"):
            continue
        failed = ", ".join(f"{field}={metric['count']}" for field, metric in data["candidate_positivity"].items() if metric["count"])
        lines.append(f"- {category}: {data['count']} queries, nonpositive/nonfinite counts {failed or 'all tested zero'}.")
    lines.extend(["", "### Exact shared-face response jumps", "", "| Orientation / class (new coverage) | Faces | Cv_t median / p99 / max | chiT_t median / p99 / max | Cp_t median / p99 / max |", "|---|---:|---|---|---|"])
    for direction in ("u", "tau"):
        for category in ("active", "pure_FreeEOS", "blend", "inactive"):
            data = result["classified"][direction]["new_coverage"]["categories"][category]
            metrics = [" / ".join(f"{data['metrics'][field][key]:.6g}" for key in ("median", "p99", "max")) if data["count"] else "none"
                for field in ("Cv_t", "chiT_t", "Cp_t")]
            lines.append(f"| {direction}/{category} | {data['count']} | " + " | ".join(metrics) + " |")
    lines.extend(["", "## Interpretation and limits", "", "Printed-node closure constrains six selected derivatives; Cp/chis/Gamma redundancy is checked but missing third/fourth twists are estimates, not supplied exact derivatives. Off-node recovered states have no source oracle. Coarse-only holdouts measure printed-table fidelity, not independent accuracy. Global even stride2 has its own explicitly recorded donor floor; no extrapolated heldouts.",
        "", "Exact face values retain C2 to numerical roundoff; response Jacobians can jump because this is not C3. Holdout departures are signed candidate-minus-printed divided by abs(printed). Faces use signed left-minus-right divided by max(abs(left/right derivative), abs(left/right response)); physical x/t are natural-log rho/T. Gamma3 is the printed Gamma3, not Gamma3-minus1. Magnitude quantiles and signed quantiles/extrema/coordinates/fractions are preserved in results.json.",
        "", "E/S signs depend on printed convention, not gauge-invariant stability. Sampled positivity is not a guarantee. Measured fractions refer to the baseline MESA blend, not an installed candidate. No source solver accuracy, gauge verification, mixture/backend safety, or recommendation to lower the runtime/default source floor follows."])
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-floor", type=float, default=3.1)
    parser.add_argument("--split-radiation", action="store_true")
    parser.add_argument("--faces", type=int, default=1000)
    parser.add_argument("--probe-child", action="store_true", help=argparse.SUPPRESS)
    arguments = parser.parse_args()
    source_floor, split = arguments.source_floor, arguments.split_radiation
    if not np.isfinite(source_floor) or source_floor < 3.1 or source_floor > 3.2 or arguments.faces < 1:
        parser.error("Bounded control requires finite source floor3.1..3.2 and a positive face count")
    output = output_directory(source_floor, split)
    if arguments.probe_child:
        run_probe.REPORT = output / "probe_runtime"
        sys.argv = [str(CORE / "mesa_probe/run_probe.py")]
        run_probe.main()
        return
    started = time.perf_counter()
    table_path = WORKSPACE / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
    archive_path = WORKSPACE / "freeeos_report/plots/offgrid_sample_data.npz"
    inputs = (table_path, archive_path, Path(__file__), CORE / "edge_control/edge_control.py", CORE / "edge_control/active_quality.py",
        CORE / "freeeos_potential.py", CORE / "validate_table.py", CORE / "active_faces.py", CORE / "run_audit.py",
        CORE / "mesa_probe/run_probe.py", CORE / "mesa_probe/build/mesa_probe", CORE / "mesa_probe/inlist_probe", CORE.parent / "mass_frac_solar.txt")
    fingerprints = {str(path): checksum(path) for path in inputs}
    _, logq, logt, rows = read_table(table_path)
    candidate = build_control(rows, logq, logt, source_floor, split)
    reference_control = build_control(rows, logq, logt, 3.2, split)
    np.testing.assert_array_equal(candidate["valid"] & (logt[None, :] >= 3.2 - 1e-12), reference_control["valid"])
    np.testing.assert_array_equal(candidate["low"][:, :, reference_control["valid"]], reference_control["low"][:, :, reference_control["valid"]])
    np.testing.assert_array_equal(candidate["nodes"][:, :, reference_control["centered_complete"]], reference_control["centered_nodes"][:, :, reference_control["centered_complete"]])
    with np.load(archive_path) as archive:
        sample_q = np.asarray(archive["logQ"]).ravel()
        sample_t = np.asarray(archive["logT"]).ravel()
    masks = {"core3.2_centered": accepted_queries(reference_control["centered_complete"], logq, logt, sample_q, sample_t),
        "control3.2_shifted": accepted_queries(reference_control["complete"], logq, logt, sample_q, sample_t),
        "tested_floor_centered": accepted_queries(candidate["centered_complete"], logq, logt, sample_q, sample_t),
        "tested_floor_shifted": accepted_queries(candidate["complete"], logq, logt, sample_q, sample_t)}
    clean = accepted_queries(candidate["valid"], logq, logt, sample_q, sample_t)
    accepted = masks["tested_floor_shifted"]
    if (accepted & ~clean).any() or accepted[sample_t < source_floor - 1e-12].any() or candidate["valid"][:, logt < source_floor - 1e-12].any():
        raise AssertionError("Coverage crossed native donor/cell or clean temperature masks")
    if (masks["control3.2_shifted"] & ~accepted).any():
        raise AssertionError("Lower clean floor lost previously accepted shifted coverage")
    recovered_indices = np.flatnonzero(accepted & ~masks["control3.2_shifted"])
    recovered_q, recovered_t = sample_q[recovered_indices], sample_t[recovered_indices]
    recovered_state = state_at(candidate["nodes"], candidate["complete"], logq, logt, recovered_q, recovered_t, split)
    populations = {"recovered_archive": {"q": recovered_q, "t": recovered_t, "overlap": np.zeros(len(recovered_q), bool),
        "state": recovered_state, "reference": recovered_state, "archive_flat_indices": recovered_indices, "tested_source_floor": source_floor},
        "native_closure": nodal_population(rows, logq, logt, candidate, reference_control, source_floor, split),
        "coarse_holdouts": holdout_population(rows, logq, logt, source_floor, split),
        **face_populations(rows, logq, logt, candidate, reference_control, source_floor, split, arguments.faces)}
    query_q = np.concatenate([data["q"] for data in populations.values()])
    query_t = np.concatenate([data["t"] for data in populations.values()])
    if (query_t < source_floor - 1e-12).any():
        raise AssertionError("Probe query batch crossed source floor")
    output.mkdir(parents=True, exist_ok=True)
    query_text = "".join(f"{density:.17g} {temperature:.17g}\n" for density, temperature in zip(query_q + 2 * query_t - 12, query_t))
    (output / "queries.txt").write_text(query_text)
    command = [sys.executable, "-B", str(Path(__file__).resolve()), "--probe-child", "--source-floor", str(source_floor)]
    if split:
        command.append("--split-radiation")
    completed = subprocess.run(command, input=query_text, capture_output=True, text=True, timeout=600)
    (output / "probe.tsv").write_text(completed.stdout)
    (output / "probe.stderr").write_text(completed.stderr)
    if completed.returncode:
        raise RuntimeError(f"MESA probe failed: {completed.returncode}; see experiment027 stderr")
    probe = np.loadtxt(io.StringIO(completed.stdout), skiprows=1, ndmin=2)
    if probe.shape != (len(query_q), 14):
        raise ValueError(f"Unexpected probe output shape: {probe.shape}")
    np.testing.assert_allclose(probe[:, :2], np.column_stack((query_q + 2 * query_t - 12, query_t)), rtol=0, atol=1e-12)
    bands = {"all": np.ones(len(sample_q), bool), "below_tested_floor": sample_t < source_floor,
        "tested_floor_to3.2": (sample_t >= source_floor) & (sample_t < 3.2), "atleast3.2": sample_t >= 3.2}
    metadata = json.loads((output / "probe_runtime/probe_metadata.json").read_text())
    result = {"method": {"tested_source_floor": source_floor, "reference_control_source_floor": 3.2, "split_radiation": split,
            "degree": 4, "contiguous_native_donors": 5, "radius": 2, "window_preference": [2, 1, 3, 0, 4],
            "source_printed_constraints_unchanged": True, "flags_density_finite_positivity_exclusions_unchanged": True,
            "meaning": "explicit clean-floor missing-twist coverage CONTROL; not independent accuracy or runtime/backend/default recommendation"},
        "provenance": {"input_sha256": fingerprints, "query_sha256": checksum(output / "queries.txt"), "probe_output_sha256": checksum(output / "probe.tsv"),
            "command": command, "python": sys.version, "numpy": np.__version__, "probe_metadata": metadata},
        "support": candidate["support"], "coarse_support": populations["coarse_holdouts"]["coarse_support"],
        "coverage": {"bands": {band: {name: int((selected & mask).sum()) for name, mask in masks.items()} for band, selected in bands.items()},
            "archived_queries": len(sample_q), "additional_over_shifted3.2": len(recovered_q), "default_centered_jets_unchanged": True,
            "source_clean_still_rejected": int((clean & ~accepted).sum())},
        "classified": {}, "source_floor_metadata_all_populations": {name: source_floor for name in populations},
        "native_closure": {"source_floor": source_floor, "excluded_nodes": populations["native_closure"]["excluded_native_material_nodes_without_complete_query_cell"],
            "printed_relative_departures": comparison(populations["native_closure"]["state"], populations["native_closure"]["reference"],
                populations["native_closure"]["q"], populations["native_closure"]["t"], tuple(populations["native_closure"]["reference"]))},
        "coarse_holdout_controls": populations["coarse_holdouts"]["comparison_controls"]}
    arrays = {"tested_source_floor": np.array(source_floor), "source_valid": candidate["valid"], "complete": candidate["complete"],
        "source_support_bounds": candidate["metadata"]["support"], **masks}
    offset = 0
    for name, population in populations.items():
        count = len(population["q"])
        local_probe = probe[offset:offset + count]
        offset += count
        groups = classified_groups(population, local_probe, name in ("u", "tau"))
        if name == "recovered_archive":
            for group in groups.values():
                for data in group["categories"].values():
                    data.pop("metrics")
                    data["interpretation"] = "positivity only; no accepted reference off-node prediction or source oracle"
        result["classified"][name] = groups
        if name in ("u", "tau"):
            result.setdefault("unclassified_face_jumps", {})[name] = face_comparison(population["state"], population["reference"], population["q"], population["t"])
        for field in ("q", "t", "overlap", "density_index", "temperature_index", "archive_flat_indices", "left_cells", "right_cells", "same_floor_centered_accepted"):
            if field in population:
                arrays[name + "_" + field] = population[field]
        arrays[name + "_probe"] = local_probe
    for name, stencil in candidate["metadata"]["stencils"].items():
        for field, values in stencil.items():
            arrays["stencil_" + name + "_" + field] = values
    result["elapsed_seconds"] = time.perf_counter() - started
    for path, expected in fingerprints.items():
        if checksum(path) != expected:
            raise RuntimeError(f"Read-only input changed during run: {path}")
    if not metadata["installed_cache_inventory_unchanged"] or metadata["private_replaced_cache_links"]:
        raise RuntimeError("MESA cache safety failed")
    (output / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    np.savez_compressed(output / "controls.npz", **arrays)
    (output / "SUMMARY.md").write_text(summary(result))
    print(json.dumps({"source_floor": source_floor, "coverage": result["coverage"], "recovered_active": result["classified"]["recovered_archive"]["new_coverage"]["categories"]["active"]["count"],
        "elapsed_seconds": result["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
