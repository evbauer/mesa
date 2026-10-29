"""Measure source compatibility at the frozen solver-face locations."""

import argparse
import json
import os
from pathlib import Path
import sys

import numpy as np


CONTROL = Path(__file__).resolve().parent
BIQUINTIC = CONTROL.parent
WORKSPACE = CONTROL.parents[4]
sys.path.insert(0, str(BIQUINTIC))

from compatibility import face_source_nodes, face_state, knot_line, nodal_closure
from freeeos_potential import (
    LN10, LOW_ORDERS, checksum, construct_twists, evaluate, read_table,
    source_constraints, thermodynamics, valid_cells,
)
from refinement_control.refinement import build_grids, material_mask, native_mask, own_use
from review.rounding_sensitivity import perturb, read_intervals, stencil_support


MANIFEST = WORKSPACE / "freeeos_biquintic_progress/experiments/038_solver_manifest/queries_v2.json"
SUPERSEDED_MANIFEST = MANIFEST.with_name("queries.json")
SOURCE = WORKSPACE / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
OUTPUT = WORKSPACE / "freeeos_biquintic_progress/experiments/040_matched_compatibility"
FROZEN_MANIFEST_SHA = "afc6fe676ba13bf858fc0a611268402b3663f4ad1e0ff991049713fd78788f19"
VARIANTS = (("split_r2_mean", True, 2, "mean"),
            ("total_r2_mean", False, 2, "mean"),
            ("split_r1_mean", True, 1, "mean"),
            ("split_r3_mean", True, 3, "mean"),
            ("split_r2_first", True, 2, "first"),
            ("split_r2_second", True, 2, "second"))
ROUNDING_CASES = ("u_ordinary_06", "u_tail_01", "tau_ordinary_06", "tau_tail_01")
ROUNDING_SEEDS = (2026102901, 2026102902, 2026102903, 2026102904)


def source_grids(rows, logq, logt):
    return {split: source_constraints(rows, logq, logt, split_radiation=split, minimum_logT=3.2)
            for split in (True, False)}


def candidate(low, valid, logq, logt, radius, route):
    spacings = (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10
    return construct_twists(low, valid, *spacings, radius=radius, route=route)[:2]


def validate_composition(header, manifest):
    if not np.allclose(header[[1, 2]], [0.70, 0.02], rtol=0, atol=1e-12):
        raise ValueError("Unexpected source-table composition")
    if not np.allclose([manifest["candidate"]["X"], manifest["candidate"]["Z"]],
                       header[[1, 2]], rtol=0, atol=1e-12):
        raise ValueError("Canonical manifest composition disagrees with source header")


def validate_manifest(manifest, header, logq, logt, split_complete):
    if manifest["schema"] != "FREEEOS_CONSUMER_MANIFEST_V2" or len(manifest["cases"]) != 24:
        raise ValueError("Unexpected frozen manifest schema or population")
    if checksum(MANIFEST) != FROZEN_MANIFEST_SHA:
        raise ValueError("Canonical frozen manifest SHA mismatch")
    old_manifest = json.loads(SUPERSEDED_MANIFEST.read_text())
    if (checksum(SUPERSEDED_MANIFEST) != manifest["supersedes"]["sha256"]
            or old_manifest["cases"] != manifest["cases"]):
        raise ValueError("Metadata correction changed a selected case")
    if manifest["provenance"]["source_sha256"] != checksum(SOURCE):
        raise ValueError("Manifest/source identity mismatch")
    if manifest["provenance"]["core_sha256"] != checksum(BIQUINTIC / "freeeos_potential.py"):
        raise ValueError("Manifest/kernel identity mismatch")
    validate_composition(header, manifest)
    seen = set()
    counts = {(direction, category): 0 for direction in ("u", "tau") for category in ("ordinary", "tail")}
    cells = valid_cells(split_complete)
    for case in manifest["cases"]:
        case_id = case["case_id"]
        if case_id in seen:
            raise ValueError("Repeated case id")
        seen.add(case_id)
        direction, category = case["direction"], case["category"]
        counts[direction, category] += 1
        left, right = np.asarray(case["left_cell"]), np.asarray(case["right_cell"])
        if not cells[tuple(left)] or not cells[tuple(right)]:
            raise ValueError(f"Incomplete frozen case {case_id}")
        expected_step = np.array([1, 0] if direction == "u" else [0, 1])
        if not np.array_equal(right - left, expected_step):
            raise ValueError(f"Nonadjacent frozen case {case_id}")
        face = logq[right[0]] if direction == "u" else logt[right[1]]
        coordinate = case["logQ"] if direction == "u" else case["logT"]
        if abs(face - coordinate) > 1e-12:
            raise ValueError(f"Face coordinate mismatch {case_id}")
        if case["fraction_FreeEOS_at_seed"] < 1 - 1e-12:
            raise ValueError(f"Case is not pure at seed {case_id}")
    if any(count != 6 for count in counts.values()):
        raise ValueError("Frozen manifest has unbalanced strata")
    return {f"{direction}_{category}": count for (direction, category), count in counts.items()}


def summarize(values):
    values = np.abs(np.asarray(values, dtype=float))
    if values.size == 0:
        return {"count": 0}
    return {"count": int(values.size), "median": float(np.median(values)),
            "p90": float(np.quantile(values, .9)), "max": float(np.max(values))}


def measure_faces(manifest, sources, rows, logq, logt):
    cases = {case["case_id"]: {"case_id": case["case_id"], "direction": case["direction"],
                               "category": case["category"], "logQ": case["logQ"],
                               "logT": case["logT"], "logRho": case["logRho"],
                               "archive_index": case["archive_index"],
                               "archived_scaled_Cv_t_jump": case["archived_scaled_Cv_t_jump"],
                               "fraction_FreeEOS_at_seed": case["fraction_FreeEOS_at_seed"]}
             for case in manifest["cases"]}
    for label, split, radius, route in VARIANTS:
        low, valid, source = sources[split]
        nodes, complete = candidate(low, valid, logq, logt, radius, route)
        cells = valid_cells(complete)
        for case in manifest["cases"]:
            record = cases[case["case_id"]]
            if not cells[tuple(case["left_cell"])] or not cells[tuple(case["right_cell"])]:
                record[label] = {"status": "incomplete_support"}
                continue
            state = face_state(nodes, complete, logq, logt, case, split)
            record[label] = {"status": "complete", **state}
            knots = knot_line(case, low, logq, logt, state["scale"])
            tangent_key = "logT" if case["direction"] == "u" else "logQ"
            knot_jumps = [face_state(nodes, complete, logq, logt,
                                     {**case, tangent_key: knot["tangent_coordinate"]}, split)["Cv_t_jump"]
                          for knot in knots]
            record[label]["knot_evaluator_Cv_t_jumps"] = knot_jumps
            record[label]["knot_formula_max_discrepancy_over_seed_scale"] = max(
                abs(actual - knot["Cv_t_left_minus_right"]) / state["scale"]
                for actual, knot in zip(knot_jumps, knots))
            if label == "split_r2_mean":
                record["source_knot_lines"] = knots
                record["source_nodal_closure_max_relative"] = nodal_closure(rows, source, face_source_nodes(case))
                record["manifest_scaled_difference"] = state["scaled_Cv_t_jump"] - case["split_scaled_Cv_t_jump"]
                record["manifest_side_max_relative_difference"] = max(
                    abs(state[side][field] - case[f"split_{side}"][field]) /
                    max(1, abs(case[f"split_{side}"][field]))
                    for side in ("left", "right") for field in ("Cv", "Cv_t", "Cv_x"))
            elif label == "total_r2_mean":
                record["total_source_knot_lines"] = knots
        del nodes, complete
    low, valid, _ = sources[True]
    nodes, complete = candidate(low, valid, logq, logt, 2, "mean")
    nodes[2, 1] = 0
    nodes[1, 2] = 0
    nodes[2, 2] = 0
    for case in manifest["cases"]:
        record = cases[case["case_id"]]
        zero_twist = face_state(nodes, complete, logq, logt, case, True)
        base = record["split_r2_mean"]
        record["zero_twist_algebraic_control"] = {
            "normal_third_jump": zero_twist["normal_third_jump"],
            "Cv_t_jump": zero_twist["Cv_t_jump"],
            "twist_contribution_Cv_t_jump": base["Cv_t_jump"] - zero_twist["Cv_t_jump"],
            "twist_contribution_over_seed_scale": (base["Cv_t_jump"] - zero_twist["Cv_t_jump"]) / base["scale"],
        }
    return cases


def local_holdouts(rows, logq, logt, manifest):
    grids = build_grids(rows, logq, logt)
    density_indices, temperature_indices = np.indices(rows.shape[:2])
    eligible = (density_indices % 2 == 1) & (temperature_indices % 2 == 1)
    eligible &= native_mask(rows, logq, logt) & material_mask(rows, logq, logt)
    eligible &= own_use(logq[:, None], logt[None, :])
    coordinates = np.nonzero(eligible)
    common = np.ones(len(coordinates[0]), dtype=bool)
    for grid in grids.values():
        from freeeos_potential import query_cells
        density, temperature, inside = query_cells(grid["q"], grid["t"],
                                                   logq[coordinates[0]], logt[coordinates[1]])
        common &= inside & grid["material_cells"][density, temperature]
        common &= grid["candidates"]["BQ_split"]["cells"][density, temperature]
    available = np.column_stack((coordinates[0][common], coordinates[1][common]))
    results = {}
    for case in manifest["cases"]:
        center_q = (case["logQ"] - logq[0]) / (logq[1] - logq[0])
        center_t = (case["logT"] - logt[0]) / (logt[1] - logt[0])
        distance = np.hypot(available[:, 0] - center_q, available[:, 1] - center_t)
        choice = np.flatnonzero((abs(available[:, 0] - center_q) <= 8)
                                    & (abs(available[:, 1] - center_t) <= 8))
        choice = choice[np.lexsort((available[choice, 1], available[choice, 0], distance[choice]))[:4]]
        records = []
        for density_index, temperature_index in available[choice]:
            sample_q = np.array([logq[density_index]])
            sample_t = np.array([logt[temperature_index]])
            reference = float(rows[density_index, temperature_index, 7])
            errors = {}
            for stride, grid in grids.items():
                candidate_grid = grid["candidates"]["BQ_split"]
                state = thermodynamics(evaluate(candidate_grid["nodes"], candidate_grid["complete"],
                                                grid["q"], grid["t"], sample_q, sample_t,
                                                orders=LOW_ORDERS, add_radiation=True), sample_q, sample_t)
                prediction = float(state["Cv"][0])
                errors[f"stride{stride}"] = abs(prediction - reference) / abs(reference)
            records.append({"source_node": [int(density_index), int(temperature_index)],
                            "logQ": float(sample_q[0]), "logT": float(sample_t[0]),
                            "source_Cv": reference, "relative_archived_Cv_errors": errors})
        results[case["case_id"]] = records
    del grids
    return results


def rounding_control(rows, logq, logt, manifest, baseline_cases, sources):
    selected_cases = [case for case in manifest["cases"] if case["case_id"] in ROUNDING_CASES]
    selected = np.zeros(rows.shape[:2], dtype=bool)
    for case in selected_cases:
        selected[face_source_nodes(case)] = True
    support = stencil_support(selected, sources[True][1])
    indices, _, lower, upper, examples = read_intervals(SOURCE, rows, support)
    changes = {case["case_id"]: {"seed_changes_over_baseline_scale": [],
                                 "knot_changes_over_baseline_scale": []} for case in selected_cases}
    for seed in ROUNDING_SEEDS:
        changed = perturb(rows, indices, lower, upper, (1, 2, 3, 4, 5, 7), seed)
        low, valid, _ = source_constraints(changed, logq, logt, split_radiation=True, minimum_logT=3.2)
        if not np.array_equal(valid, sources[True][1]):
            raise ValueError("Printed rounding changed validity mask")
        nodes, complete = candidate(low, valid, logq, logt, 2, "mean")
        for case in selected_cases:
            case_id = case["case_id"]
            base = baseline_cases[case_id]
            actual = face_state(nodes, complete, logq, logt, case, True)
            change = (actual["Cv_t_jump"] - base["split_r2_mean"]["Cv_t_jump"]) / base["split_r2_mean"]["scale"]
            knots = knot_line(case, low, logq, logt, base["split_r2_mean"]["scale"])
            knot_change = max(abs(new["Cv_t_jump_over_seed_scale"] - old["Cv_t_jump_over_seed_scale"])
                              for new, old in zip(knots, base["source_knot_lines"]))
            changes[case_id]["seed_changes_over_baseline_scale"].append(float(change))
            changes[case_id]["knot_changes_over_baseline_scale"].append(float(knot_change))
        del changed, low, nodes, complete
    for case_id, record in changes.items():
        max_change = max(abs(value) for value in record["seed_changes_over_baseline_scale"])
        base = baseline_cases[case_id]["split_r2_mean"]["scaled_Cv_t_jump"]
        record["max_abs_seed_change_over_scale"] = max_change
        record["max_seed_change_over_abs_baseline_jump"] = max_change / abs(base)
        record["max_abs_knot_change_over_scale"] = max(record["knot_changes_over_baseline_scale"])
    return {"cases": changes, "seeds": list(ROUNDING_SEEDS), "support_valid_nodes": int(support.sum()),
            "printed_interval_examples": examples,
            "model": "Independent uniform perturbation of six consumed printed fields within their displayed half-decimal intervals; four deterministic seeds; limited formatting sensitivity, not source uncertainty or an interval bound"}


def make_figure(cases, output):
    os.environ.setdefault("MPLCONFIGDIR", str(output / ".mplconfig"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), constrained_layout=True)
    for direction, marker in (("u", "o"), ("tau", "s")):
        for category, color in (("ordinary", "#277da1"), ("tail", "#d1495b")):
            group = [case for case in cases.values() if case["direction"] == direction and case["category"] == category]
            x = [max(abs(knot["Cv_t_jump_over_seed_scale"]) for knot in case["source_knot_lines"]) for case in group]
            y = [abs(case["split_r2_mean"]["scaled_Cv_t_jump"]) for case in group]
            axes[0].scatter(x, y, marker=marker, color=color, label=f"{direction} {category}", s=39)
            sensitivity = []
            for case in group:
                baseline = case["split_r2_mean"]["scaled_Cv_t_jump"]
                variants = [case[name]["Cv_t_jump"] / case["split_r2_mean"]["scale"] for name, _, _, _ in VARIANTS[2:]
                            if case[name]["status"] == "complete"]
                sensitivity.append(max(abs(value - baseline) for value in variants) if variants else np.nan)
            axes[1].scatter(y, sensitivity, marker=marker, color=color, s=39)
    for axis in axes:
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.grid(True, which="both", alpha=.2)
    axes[0].plot([1e-9, 1], [1e-9, 1], color=".6", linestyle="--", linewidth=1)
    axes[0].set_xlabel("max knot-line |source-implied Cv,t jump| / seed scale")
    axes[0].set_ylabel("off-node |Cv,t jump| / seed scale")
    axes[1].set_xlabel("off-node |Cv,t jump| / seed scale")
    axes[1].set_ylabel("max radius/route change / seed scale")
    axes[0].legend(fontsize=8)
    figure.suptitle("Frozen pure-at-seed FreeEOS faces: source-line compatibility and twist sensitivity")
    figure.savefig(output / "matched_compatibility.pdf")
    figure.savefig(output / "matched_compatibility.png", dpi=180)
    plt.close(figure)


def run(output):
    manifest = json.loads(MANIFEST.read_text())
    header, logq, logt, rows = read_table(SOURCE)
    sources = source_grids(rows, logq, logt)
    baseline_nodes, baseline_complete = candidate(*sources[True][:2], logq, logt, 2, "mean")
    strata = validate_manifest(manifest, header, logq, logt, baseline_complete)
    del baseline_nodes, baseline_complete
    cases = measure_faces(manifest, sources, rows, logq, logt)
    holdouts = local_holdouts(rows, logq, logt, manifest)
    rounding = rounding_control(rows, logq, logt, manifest, cases, sources)
    holdout_pairs = [target for targets in holdouts.values() for target in targets]
    unique_holdout_nodes = {tuple(target["source_node"]) for target in holdout_pairs}
    result = {
        "status": "completed frozen-case compatibility diagnostic; not source accuracy or solver result",
        "provenance": {"manifest_sha256": checksum(MANIFEST), "expected_frozen_manifest_sha256": FROZEN_MANIFEST_SHA,
                       "superseded_manifest_sha256": checksum(SUPERSEDED_MANIFEST),
                       "cases_identical_to_superseded_manifest": True,
                       "source_sha256": checksum(SOURCE), "core_sha256": checksum(BIQUINTIC / "freeeos_potential.py"),
                       "script_sha256": checksum(__file__), "python": sys.version, "numpy": np.__version__,
                       "command": sys.argv, "node": os.uname().nodename, "job_id": os.environ.get("SLURM_JOB_ID")},
        "composition": {"source_header_X": float(header[1]), "source_header_Z": float(header[2]),
                        "manifest_candidate_X": manifest["candidate"]["X"], "manifest_candidate_Z": manifest["candidate"]["Z"],
                        "manifest_candidate_metadata_matches_source": bool(np.allclose(
                            [manifest["candidate"]["X"], manifest["candidate"]["Z"]], header[[1, 2]]))},
        "definitions": {"source_knot_jump": "right-minus-left native normal third derivative from Appendix A using native natural-log spacing and low-order supplied nodal constraints",
                        "physical_jump": "left-minus-right Cv partial at fixed physical ln rho, ln T; u normal has Cv_t=8*g_uuu/T and Cv_x=-4*g_uuu/T; tau normal has Cv_t=-g_tautautau/T and Cv_x=0",
                        "seed_scale": "max(1,abs(left/right Cv),abs(left/right Cv_t)) in cgs units",
                        "source_closure": "max relative discrepancy at six neighboring source nodes, redundant printed mixed-energy/mixed-entropy/Cv-entropy routes",
                        "zero_twist": "algebraic low-order-only control; missing mixed nodal derivatives set to zero, not a recommended interpolant",
                        "holdouts": "up to four nearest odd/odd archived source nodes per seed within 8 native index steps in each direction, common accepted material/own-use and split-complete target across strides 2/4/8; source Cv column reference; stride1 retained and not a holdout"},
        "strata": strata, "variants": [name for name, _, _, _ in VARIANTS],
        "cases": cases, "local_holdouts": holdouts,
        "local_holdout_population": {"case_target_associations": len(holdout_pairs),
                                     "unique_withheld_source_nodes": len(unique_holdout_nodes),
                                     "covered_cases": sum(bool(targets) for targets in holdouts.values()),
                                     "not_off_node_seed_accuracy": True},
        "rounding_control": rounding,
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    make_figure(cases, output)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    if arguments.output.resolve() != OUTPUT.resolve():
        parser.error("Outputs are restricted to experiment 040")
    result = run(arguments.output)
    print(json.dumps({"case_count": len(result["cases"]),
                      "manifest_sha256": result["provenance"]["manifest_sha256"],
                      "holdout_targets": sum(len(targets) for targets in result["local_holdouts"].values()),
                      "output": str(arguments.output)}, indent=2))


if __name__ == "__main__":
    main()
