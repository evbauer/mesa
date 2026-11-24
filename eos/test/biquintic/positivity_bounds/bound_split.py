"""Reference-split polynomial plus monotone analytic radiation bounds; sidecar only."""

import os
import sys

sys.dont_write_bytecode = True
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import argparse
import datetime
import json
from pathlib import Path
import time

import numpy as np

from bounds import (
    DTYPE, FIELDS, SAFETY, STATUS, bernstein_bounds, bernstein_value, cell_power,
    normalize_closed_cell, polynomial_value, sign_powers, subdivide,
)
from run_bounds import BIQUINTIC, OUTPUT as TOTAL_OUTPUT, PHASES, TABLE, TABLE_SHA, WORKSPACE, collect_jets, location
from freeeos_potential import (
    CRAD, LN10, checksum, construct_twists, evaluate, read_table,
    source_constraints, thermodynamics, valid_cells,
)

OUTPUT = WORKSPACE / "freeeos_biquintic_progress/experiments/028_split_bounds"


def radiation_cv_numerator(logq, logt):
    return 4 * DTYPE(CRAD) * np.power(DTYPE(10), 2*np.asarray(logt, dtype=DTYPE)-np.asarray(logq, dtype=DTYPE)+12)


def radiation_bounds(logbox):
    density_min, density_max, temperature_min, temperature_max = logbox
    return (radiation_cv_numerator(density_max, temperature_min),
            radiation_cv_numerator(density_min, temperature_max))


def radiation_guard(field, maximum_radiation_cv, restored_core):
    analytic_margin = SAFETY * DTYPE(np.finfo(DTYPE).eps) * maximum_radiation_cv if field == "Cv_numerator" else DTYPE(0)
    restoration_factor = DTYPE(20)/12 if field == "Cv_numerator" else DTYPE(2)/12
    restoration_margin = (SAFETY * DTYPE(np.finfo(float).eps) * maximum_radiation_cv * restoration_factor
                          if restored_core else DTYPE(0))
    return analytic_margin + restoration_margin


def adaptive_split(coefficients, polynomial_margin, logbox, field, witness=None,
                   restored_core=True, max_depth=4, max_nodes=341):
    if field not in FIELDS or not 0 <= max_depth <= 4 or not 1 <= max_nodes <= 341:
        raise ValueError("Require known field, depth 0..4, nodes 1..341")
    radiation_min, radiation_max = radiation_bounds(logbox)
    initial_lower = coefficients.min() + (radiation_min if field == "Cv_numerator" else 0) - polynomial_margin - radiation_guard(field, radiation_max, restored_core)
    stack = [(coefficients, DTYPE(polynomial_margin), 0, (0., 1., 0., 1.), tuple(logbox))]
    visited, positive_leaves, unresolved_leaves, deepest, witness_queries = 0, 0, 0, 0, 0
    minimum_lower = DTYPE(np.inf)
    while stack:
        current, current_poly_margin, depth, phasebox, current_logbox = stack.pop()
        visited += 1
        deepest = max(deepest, depth)
        radiation_min, radiation_max = radiation_bounds(current_logbox)
        total_margin = current_poly_margin + radiation_guard(field, radiation_max, restored_core)
        lower = current.min() + (radiation_min if field == "Cv_numerator" else 0) - total_margin
        if lower > 0:
            minimum_lower = min(minimum_lower, lower)
            positive_leaves += 1
            continue
        density_min, density_max, temperature_min, temperature_max = current_logbox
        phase_density_min, phase_density_max, phase_temperature_min, phase_temperature_max = phasebox
        if witness is not None:
            control_density, control_temperature = np.unravel_index(np.argmin(current), (6, 6))
            candidates = [(control_density/5, control_temperature/5), (.5, .5),
                          (0., 0.), (0., 1.), (1., 0.), (1., 1.)]
            for density_phase, temperature_phase in dict.fromkeys(candidates):
                sample_q = density_min + (density_max-density_min)*density_phase
                sample_t = temperature_min + (temperature_max-temperature_min)*temperature_phase
                value = bernstein_value(current, DTYPE(density_phase), DTYPE(temperature_phase))
                if field == "Cv_numerator":
                    value += radiation_cv_numerator(sample_q, sample_t)
                if value < -total_margin:
                    evidence = witness(phase_density_min+(phase_density_max-phase_density_min)*density_phase,
                                       phase_temperature_min+(phase_temperature_max-phase_temperature_min)*temperature_phase)
                    witness_queries += 1
                    if evidence is not None:
                        return {"status": "negative_witness", "visited_nodes": visited,
                                "maximum_depth": deepest, "positive_leaves": positive_leaves,
                                "unresolved_leaves": 0, "lower_bound": float(initial_lower),
                                "witness_core_queries": witness_queries, "witness": evidence}
        if depth == max_depth or visited+len(stack)+4 > max_nodes:
            unresolved_leaves += 1
            minimum_lower = min(minimum_lower, lower)
            continue
        density_middle, temperature_middle = (density_min+density_max)/2, (temperature_min+temperature_max)/2
        phase_density_middle = (phase_density_min+phase_density_max)/2
        phase_temperature_middle = (phase_temperature_min+phase_temperature_max)/2
        child_margin = current_poly_margin + SAFETY * DTYPE(np.finfo(DTYPE).eps)*np.max(abs(current))
        logboxes = [(density_min, density_middle, temperature_min, temperature_middle),
                    (density_min, density_middle, temperature_middle, temperature_max),
                    (density_middle, density_max, temperature_min, temperature_middle),
                    (density_middle, density_max, temperature_middle, temperature_max)]
        phaseboxes = [(phase_density_min, phase_density_middle, phase_temperature_min, phase_temperature_middle),
                      (phase_density_min, phase_density_middle, phase_temperature_middle, phase_temperature_max),
                      (phase_density_middle, phase_density_max, phase_temperature_min, phase_temperature_middle),
                      (phase_density_middle, phase_density_max, phase_temperature_middle, phase_temperature_max)]
        stack.extend((child, child_margin, depth+1, child_phasebox, child_logbox)
                     for child, child_phasebox, child_logbox in zip(subdivide(current), phaseboxes, logboxes))
    return {"status": "positive" if unresolved_leaves == 0 else "unresolved",
            "visited_nodes": visited, "maximum_depth": deepest, "positive_leaves": positive_leaves,
            "unresolved_leaves": unresolved_leaves, "lower_bound": float(minimum_lower),
            "witness_core_queries": witness_queries, "witness": None}


def split_witness(field, density_index, temperature_index, nodes, complete, logq, logt):
    def witness(density_phase, temperature_phase):
        sample_q = np.array([logq[density_index]+(logq[density_index+1]-logq[density_index])*density_phase])
        sample_t = np.array([logt[temperature_index]+(logt[temperature_index+1]-logt[temperature_index])*temperature_phase])
        sample_q = np.clip(sample_q, logq[density_index], logq[density_index+1])
        sample_t = np.clip(sample_t, logt[temperature_index], logt[temperature_index+1])
        native = evaluate(nodes, complete, logq, logt, sample_q, sample_t, add_radiation=True,
                          cell_indices=(np.array([density_index]), np.array([temperature_index])))
        state = thermodynamics(native, sample_q, sample_t)
        value = state["Cv" if field == "Cv_numerator" else "P_rho"][0]
        if np.isfinite(value) and value < 0:
            return {"logQ": float(sample_q[0]), "logT": float(sample_t[0]),
                    "logRho": float(sample_q[0]+2*sample_t[0]-12),
                    "core_field": "Cv" if field == "Cv_numerator" else "P_rho",
                    "core_value": float(value), "core_Cv": float(state["Cv"][0]),
                    "core_P_rho": float(state["P_rho"][0]), "core_add_radiation": True}
        return None
    return witness


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-depth", type=int, choices=range(5), default=4)
    parser.add_argument("--max-nodes", type=int, default=341)
    parser.add_argument("--batch-size", type=int, default=2048)
    arguments = parser.parse_args()
    if not 1 <= arguments.max_nodes <= 341 or not 1 <= arguments.batch_size <= 4096:
        parser.error("Bounded node budget 1..341 and batch size 1..4096 required")
    started, started_at = time.monotonic(), datetime.datetime.now().astimezone().isoformat()
    baseline = json.loads((TOTAL_OUTPUT / "results.json").read_text())
    baseline_archive = np.load(TOTAL_OUTPUT / "cell_bounds.npz", allow_pickle=False)
    baseline_before = {path.name: checksum(path) for path in (TOTAL_OUTPUT / "results.json", TOTAL_OUTPUT / "cell_bounds.npz")}
    if baseline_before["cell_bounds.npz"] != baseline["artifacts_sha256"]["cell_bounds.npz"]:
        raise ValueError("Experiment024 archive does not match its calibrated metadata")
    if not baseline["calibration_passed"] or checksum(TABLE) != TABLE_SHA:
        raise ValueError("Require calibrated experiment024 and the declared unchanged source")
    before = {str(path.relative_to(WORKSPACE)): checksum(path) for path in (TABLE, BIQUINTIC / "freeeos_potential.py")}
    _, logq, logt, rows = read_table(TABLE)
    low, source_valid, state = source_constraints(rows, logq, logt, split_radiation=True)
    du, dtau = (logq[1]-logq[0])*LN10, (logt[1]-logt[0])*LN10
    nodes, complete, _ = construct_twists(low, source_valid, du, dtau, radius=2, route="mean")
    density_indices, temperature_indices = np.nonzero(valid_cells(complete))
    np.testing.assert_array_equal(density_indices, baseline_archive["density_cell_index"])
    np.testing.assert_array_equal(temperature_indices, baseline_archive["temperature_cell_index"])
    count = len(density_indices)
    if count > 150000:
        raise ValueError("150000-cell budget exceeded")
    full_own, material = baseline_archive["full_one_branch_own_use"], baseline_archive["material_corner_screen"]
    subsets = {"all_complete": np.ones(count, dtype=bool), "full_one_branch_own_use": full_own,
               "material_full_one_branch_own_use": material & full_own,
               "material_not_full_one_branch_own_use": material & ~full_own,
               "not_full_one_branch_own_use": ~full_own}
    archive = {"density_cell_index": density_indices, "temperature_cell_index": temperature_indices,
               "full_one_branch_own_use": full_own, "material_corner_screen": material,
               "maximum_corner_Prad_over_Pgas": baseline_archive["maximum_corner_Prad_over_Pgas"]}
    records = {field: {name: np.zeros(count, dtype=dtype) for name, dtype in (
        ("status", np.int8), ("representation_status", np.int8), ("root_polynomial_minimum", DTYPE),
        ("root_polynomial_maximum", DTYPE), ("polynomial_margin", DTYPE), ("restoration_margin", DTYPE),
        ("analytic_margin", DTYPE), ("root_lower_bound", DTYPE), ("final_lower_bound", DTYPE),
        ("representation_lower_bound", DTYPE), ("visited_nodes", np.int16),
        ("representation_visited_nodes", np.int16), ("maximum_depth", np.int8),
        ("positive_leaves", np.int16), ("unresolved_leaves", np.int16), ("witness_core_queries", np.int16))}
               for field in FIELDS}
    witnesses, unresolved = {field: [] for field in FIELDS}, {field: [] for field in FIELDS}
    agreement = {field: {"sample_count": 0, "outside_effective_margin_count": 0,
                         "outside_polynomial_margin_count": 0, "negative_core_samples": 0,
                         "zero_core_samples": 0, "nonfinite_core_samples": 0,
                         "maximum_difference_over_effective_margin": 0., "maximum_difference_over_polynomial_margin": 0.,
                         "representation_positive_but_core_nonpositive_samples": 0, "worst": None} for field in FIELDS}
    print(f"{count} identical cells; declared reference split, dual representation/restored-core guards", flush=True)
    for start in range(0, count, arguments.batch_size):
        stop = min(count, start+arguments.batch_size)
        density_batch, temperature_batch = density_indices[start:stop], temperature_indices[start:stop]
        power, magnitude = cell_power(collect_jets(nodes, density_batch, temperature_batch), du, dtau)
        powers, scales = sign_powers(power, magnitude, du, dtau)
        density_width = logq[density_batch+1].astype(DTYPE)-logq[density_batch].astype(DTYPE)
        temperature_width = logt[temperature_batch+1].astype(DTYPE)-logt[temperature_batch].astype(DTYPE)
        density_stretch, temperature_stretch = density_width/DTYPE(logq[1]-logq[0]), temperature_width/DTYPE(logt[1]-logt[0])
        powers = {field: normalize_closed_cell(values, density_stretch, temperature_stretch) for field, values in powers.items()}
        scales = {field: normalize_closed_cell(values, density_stretch, temperature_stretch) for field, values in scales.items()}
        radiation_min = radiation_cv_numerator(logq[density_batch+1], logt[temperature_batch])
        radiation_max = radiation_cv_numerator(logq[density_batch], logt[temperature_batch+1])
        archive.setdefault("radiation_Cv_numerator_minimum", np.zeros(count, dtype=DTYPE))[start:stop] = radiation_min
        archive.setdefault("radiation_Cv_numerator_maximum", np.zeros(count, dtype=DTYPE))[start:stop] = radiation_max
        margins, effective_margins = {}, {}
        for field in FIELDS:
            coefficients, polynomial_margin, *_ = bernstein_bounds(powers[field], scales[field])
            if not np.all(np.isfinite(coefficients)) or not np.all(np.isfinite(polynomial_margin)):
                raise AssertionError("Nonfinite split polynomial")
            analytic_margin = radiation_guard(field, radiation_max, False)
            restoration_margin = radiation_guard(field, radiation_max, True)-analytic_margin
            analytic_lower = radiation_min if field == "Cv_numerator" else DTYPE(0)
            representation_lower = coefficients.min(axis=(1, 2))+analytic_lower-polynomial_margin-analytic_margin
            root_lower = representation_lower-restoration_margin
            record = records[field]
            for name, values in (("root_polynomial_minimum", coefficients.min(axis=(1, 2))),
                                 ("root_polynomial_maximum", coefficients.max(axis=(1, 2))),
                                 ("polynomial_margin", polynomial_margin), ("analytic_margin", analytic_margin),
                                 ("restoration_margin", restoration_margin), ("root_lower_bound", root_lower),
                                 ("final_lower_bound", root_lower), ("representation_lower_bound", representation_lower)):
                record[name][start:stop] = values
            record["status"][start:stop] = (root_lower > 0).astype(np.int8)
            record["representation_status"][start:stop] = (representation_lower > 0).astype(np.int8)
            record["visited_nodes"][start:stop] = 1
            record["representation_visited_nodes"][start:stop] = 1
            record["positive_leaves"][start:stop] = (root_lower > 0).astype(np.int16)
            margins[field], effective_margins[field] = polynomial_margin+analytic_margin, polynomial_margin+analytic_margin+restoration_margin
            for local_index in np.flatnonzero(root_lower <= 0):
                index = start+local_index
                density_index, temperature_index = int(density_batch[local_index]), int(temperature_batch[local_index])
                logbox = tuple(map(DTYPE, (logq[density_index], logq[density_index+1], logt[temperature_index], logt[temperature_index+1])))
                callback = split_witness(field, density_index, temperature_index, nodes, complete, logq, logt)
                result = adaptive_split(coefficients[local_index], polynomial_margin[local_index], logbox, field, callback,
                                        True, arguments.max_depth, arguments.max_nodes)
                record["status"][index] = STATUS[result["status"]]
                record["final_lower_bound"][index] = result["lower_bound"]
                for name in ("visited_nodes", "maximum_depth", "positive_leaves", "unresolved_leaves", "witness_core_queries"):
                    record[name][index] = result[name]
                if representation_lower[local_index] <= 0:
                    representation = adaptive_split(coefficients[local_index], polynomial_margin[local_index], logbox, field,
                                                    callback, False, arguments.max_depth, arguments.max_nodes)
                    record["representation_status"][index] = STATUS[representation["status"]]
                    record["representation_lower_bound"][index] = representation["lower_bound"]
                    record["representation_visited_nodes"][index] = representation["visited_nodes"]
                if result["status"] != "positive":
                    entry = {**location(index, density_indices, temperature_indices, logq, logt),
                             "full_one_branch_own_use": bool(full_own[index]), "material_corner_screen": bool(material[index]),
                             "representation_status": int(record["representation_status"][index]),
                             "polynomial_margin": float(polynomial_margin[local_index]),
                             "root_restoration_margin": float(restoration_margin[local_index]), **result}
                    (witnesses if result["status"] == "negative_witness" else unresolved)[field].append(entry)
        for density_phase, temperature_phase in PHASES:
            sample_q = np.asarray(logq[density_batch]+density_width*density_phase, dtype=float)
            sample_t = np.asarray(logt[temperature_batch]+temperature_width*temperature_phase, dtype=float)
            actual_density = (sample_q.astype(DTYPE)-logq[density_batch].astype(DTYPE))/density_width
            actual_temperature = (sample_t.astype(DTYPE)-logt[temperature_batch].astype(DTYPE))/temperature_width
            native = evaluate(nodes, complete, logq, logt, sample_q, sample_t, add_radiation=True,
                              cell_indices=(density_batch, temperature_batch))
            state = thermodynamics(native, sample_q, sample_t)
            for field in FIELDS:
                prediction = polynomial_value(powers[field], actual_density, actual_temperature)
                if field == "Cv_numerator":
                    prediction += radiation_cv_numerator(sample_q, sample_t)
                core_value = state["Cv"]*state["T"] if field == "Cv_numerator" else state["P_rho"]
                difference = abs(prediction-core_value)
                normalized = difference/effective_margins[field]
                polynomial_normalized = difference/margins[field]
                summary = agreement[field]
                summary["sample_count"] += len(sample_q)
                summary["outside_effective_margin_count"] += int(np.count_nonzero(normalized > 1))
                summary["outside_polynomial_margin_count"] += int(np.count_nonzero(polynomial_normalized > 1))
                summary["negative_core_samples"] += int(np.count_nonzero(core_value < 0))
                summary["zero_core_samples"] += int(np.count_nonzero(core_value == 0))
                summary["nonfinite_core_samples"] += int(np.count_nonzero(~np.isfinite(core_value)))
                summary["maximum_difference_over_polynomial_margin"] = max(summary["maximum_difference_over_polynomial_margin"], float(polynomial_normalized.max()))
                record = records[field]
                summary["representation_positive_but_core_nonpositive_samples"] += int(np.count_nonzero((record["representation_status"][start:stop] == 1) & (core_value <= 0)))
                if np.any((record["status"][start:stop] == 1) & (core_value <= 0)):
                    raise AssertionError("Restored-core guarded positive disagrees in sampled sign")
                worst = int(np.argmax(normalized))
                if normalized[worst] > summary["maximum_difference_over_effective_margin"]:
                    summary["maximum_difference_over_effective_margin"] = float(normalized[worst])
                    summary["worst"] = {**location(start+worst, density_indices, temperature_indices, logq, logt),
                                        "logQ": float(sample_q[worst]), "logT": float(sample_t[worst]),
                                        "representation_value": float(prediction[worst]), "restored_core_value": float(core_value[worst]),
                                        "absolute_difference": float(difference[worst]), "effective_margin": float(effective_margins[field][worst])}
        if start % (arguments.batch_size*10) == 0:
            print(f"{stop}/{count} bounded and restored-core checked", flush=True)
    summaries = {}
    for name, mask in subsets.items():
        summaries[name] = {"cell_count": int(mask.sum()), "fields": {}}
        both_positive, both_representation = mask.copy(), mask.copy()
        for field, record in records.items():
            summary = {status: int(np.count_nonzero(mask & (record["status"] == status_value))) for status, status_value in STATUS.items()}
            summary["representation"] = {status: int(np.count_nonzero(mask & (record["representation_status"] == status_value))) for status, status_value in STATUS.items()}
            summary["root_positive"] = int(np.count_nonzero(mask & (record["root_lower_bound"] > 0)))
            summary["adaptive_positive"] = summary["positive"]-summary["root_positive"]
            summary["restoration_guard_only_unresolved"] = int(np.count_nonzero(mask & (record["status"] == 0) & (record["representation_status"] == 1)))
            summary["visited_nodes"] = int(record["visited_nodes"][mask].sum())
            summary["representation_visited_nodes"] = int(record["representation_visited_nodes"][mask].sum())
            summary["unresolved_leaves"] = int(record["unresolved_leaves"][mask].sum())
            summary["maximum_depth"] = int(record["maximum_depth"][mask].max()) if mask.any() else 0
            summary["extrema"] = {}
            for metric, maximize in (("final_lower_bound", False), ("polynomial_margin", True), ("restoration_margin", True)):
                index = int(np.argmax(np.where(mask, record[metric], -np.inf)) if maximize else np.argmin(np.where(mask, record[metric], np.inf)))
                summary["extrema"][metric] = {"value": float(record[metric][index]), **location(index, density_indices, temperature_indices, logq, logt)}
            baseline_status = baseline_archive[f"{field}__status"]
            summary["total024_unresolved_now_positive"] = int(np.count_nonzero(mask & (baseline_status == 0) & (record["status"] == 1)))
            summary["total024_unresolved_now_representation_positive"] = int(np.count_nonzero(mask & (baseline_status == 0) & (record["representation_status"] == 1)))
            summaries[name]["fields"][field] = summary
            both_positive &= record["status"] == 1
            both_representation &= record["representation_status"] == 1
        summaries[name]["both_restored_core_guarded_positive"] = int(both_positive.sum())
        summaries[name]["both_representation_positive"] = int(both_representation.sum())
    after = {name: checksum(WORKSPACE / name) for name in before}
    if before != after:
        raise AssertionError("Core/source changed during reference-split audit")
    baseline_after = {name: checksum(TOTAL_OUTPUT / name) for name in baseline_before}
    if baseline_before != baseline_after:
        raise AssertionError("Experiment024 input changed during reference-split audit")
    prior_ambiguous = baseline_archive["P_rho__status"] == STATUS["unresolved"]
    reduction = baseline_archive["P_rho__roundoff_margin"][prior_ambiguous] / records["P_rho"]["polynomial_margin"][prior_ambiguous]
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for field, record in records.items():
        archive.update({f"{field}__{name}": values for name, values in record.items()})
    np.savez_compressed(OUTPUT / "cell_bounds.npz", **archive)
    (OUTPUT / "negative_witnesses.json").write_text(json.dumps(witnesses, indent=2, allow_nan=False)+"\n")
    (OUTPUT / "unresolved_cells.json").write_text(json.dumps(unresolved, indent=2, allow_nan=False)+"\n")
    result = {"started_at": started_at, "elapsed_seconds": time.monotonic()-started, "command": sys.argv,
              "candidate": "radius2_mean_declared_reference_split", "source_sha256": before, "source_unchanged": before == after,
              "baseline024_sha256": baseline_before, "baseline024_unchanged": baseline_before == baseline_after,
              "sidecar_sha256": {path.name: checksum(path) for path in (Path(__file__), Path(__file__).with_name("bounds.py"), Path(__file__).with_name("run_bounds.py"))},
              "CRAD_declared": float(CRAD), "floating_proof_claim": False,
              "numpy": np.__version__, "python": sys.executable,
              "prior_ambiguous_P_rho_polynomial_guard_reduction": {"cell_count": int(prior_ambiguous.sum()),
                                                                   "minimum_factor": float(reduction.min()),
                                                                   "median_factor": float(np.median(reduction)),
                                                                   "maximum_factor": float(reduction.max())},
              "radiation": {"potential": "-CRAD/3 * exp(2*tau-u+12*ln(10))",
                            "P_rho_analytic_contribution": "Exactly zero: -Frad + Frad; certify remainder polynomial",
                            "Cv_numerator_analytic_contribution": "4*CRAD*exp(2*tau-u+12*ln(10)); lower at tau_min,u_max; upper at tau_max,u_min",
                            "historical_remainder": "Existing declared-reference subtraction retained, including known ~4.4 ppm radiation remainder; not verified material EOS or source-constant determination"},
              "margin_policy": {"safety": SAFETY, "longdouble_epsilon": float(np.finfo(DTYPE).eps),
                                "float64_epsilon": float(np.finfo(float).eps),
                                "polynomial": "Same absolute-contribution 4096*(eps_longdouble+eps_float64)*M as024, now on remainder jets",
                                "analytic_exp": "4096*eps_longdouble*maximum radiation Cv numerator",
                                "restored_core": "4096*eps_float64*max_rad_Cv*(20/12 for Cv numerator, 2/12 for P_rho), per leaf",
                                "scope": "Heuristic conservative numerical guards, not directed-rounded exp intervals or rigorous IEEE proof"},
              "domain": baseline["domain"], "same_cells_and_frozen_domain_as024": True,
              "budgets": {"max_depth": arguments.max_depth, "max_nodes_per_mode_cell_field": arguments.max_nodes,
                          "max_modes": 2, "batch_size": arguments.batch_size, "max_cells": 150000, "max_live_dfs_nodes": 1+3*arguments.max_depth},
              "counts": summaries, "sample_agreement_add_radiation_true": agreement,
              "calibration_passed": all(summary["outside_effective_margin_count"] == 0 and summary["nonfinite_core_samples"] == 0 for summary in agreement.values()),
              "status_codes": STATUS,
              "limitations": ["Reference split is not a verified material EOS, source truth, global stability, actual radiation constant or runtime claim",
                              "Representation positivity is separate from restored-core guarded positivity and floating cancellation",
                              "No Pgas/E/S certificate; material mask is still source four-corner screening only; partial own-use cells reported separately",
                              "Zero/nonpositive bound minima are ambiguous; only actual finite negative restored-core points are witnesses"],
              "artifacts_sha256": {path.name: checksum(path) for path in (OUTPUT / "cell_bounds.npz", OUTPUT / "negative_witnesses.json", OUTPUT / "unresolved_cells.json")}}
    (OUTPUT / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Experiment028: declared-reference split positivity bounds", "",
             f"Completed in {result['elapsed_seconds']:.1f} seconds, on the identical {count} experiment024 cells.",
             "Numerical polynomial-plus-analytic bounds, NOT rigorous IEEE intervals or verified material EOS.", "",
             "| Subset | Cells | Cv positive / witness / unresolved | P_rho positive / witness / unresolved | Both guarded positive | Both representation positive |",
             "| --- | ---: | --- | --- | ---: | ---: |"]
    for name, summary in summaries.items():
        entries = [" / ".join(str(summary["fields"][field][key]) for key in ("positive", "negative_witness", "unresolved")) for field in FIELDS]
        lines.append(f"| {name} | {summary['cell_count']} | {entries[0]} | {entries[1]} | {summary['both_restored_core_guarded_positive']} | {summary['both_representation_positive']} |")
    all_prho = summaries["all_complete"]["fields"]["P_rho"]
    lines += ["", "P_rho analytically cancels radiation exactly; Cv adds monotone positive 4*CRAD*exp(2*tau-u+12*ln(10)) on every leaf.",
              f"Of024's {baseline['counts']['all_complete']['fields']['P_rho']['unresolved']} unresolved P_rho cells: {all_prho['total024_unresolved_now_positive']} now guarded-positive, {all_prho['total024_unresolved_now_representation_positive']} now representation-positive.",
              f"Restoration-guard-only P_rho ambiguities: {all_prho['restoration_guard_only_unresolved']}.",
              f"Core add_radiation=True sample calibration passes: {result['calibration_passed']}; four phases in every cell.",
              "Declared modern CRAD and the existing ~4.4 ppm historical radiation remainder remain unchanged; this is not a verified material EOS.",
              "Read results.json for sample zero/negative/cancellation counts, guards, exact worst coordinates and source/baseline SHA.",
              "Lead integration: use the material/full-one-branch primary row; keep representation-only results and restored-core guards distinct.",
              "Experiment024 is read-only to this extension. No source truth, global/blended/runtime stability or Pgas/E/S certificate."]
    (OUTPUT / "SUMMARY.md").write_text("\n".join(lines)+"\n")
    print("\n".join(lines), flush=True)
    if not result["calibration_passed"]:
        raise SystemExit("Restored-core calibration failed; do not use guarded certificates")


if __name__ == "__main__":
    main()
