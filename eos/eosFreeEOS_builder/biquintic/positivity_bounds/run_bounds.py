"""Bounded total-only cell positivity sidecar; no source/runtime stability claim."""

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
    DTYPE, FIELDS, SAFETY, STATUS, adaptive_bound, bernstein_bounds, cell_power,
    derivative, domain_masks, normalize_closed_cell, polynomial_value, sign_powers,
)

SIDECAR = Path(__file__).resolve().parent
BIQUINTIC = SIDECAR.parent
WORKSPACE = BIQUINTIC.parents[3]
sys.path.insert(0, str(BIQUINTIC))

from freeeos_potential import (
    CRAD, LN10, checksum, construct_twists, evaluate, read_table,
    source_constraints, thermodynamics, valid_cells,
)

TABLE = WORKSPACE / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
OUTPUT = WORKSPACE / "freeeos_biquintic_progress/experiments/024_positivity_bounds"
TABLE_SHA = "e6ca7ce27bee8197bb257a8c4886a06366e6951e58958fbed3735c877ee68502"
PHASES = ((.5, .5), (.211, .337), (.789, .663), (.01, .99))


def collect_jets(nodes, density_indices, temperature_indices):
    jets = np.empty((len(density_indices), 2, 2, 3, 3), dtype=DTYPE)
    for density_end in range(2):
        for temperature_end in range(2):
            jets[:, density_end, temperature_end] = nodes[
                :, :, density_indices+density_end, temperature_indices+temperature_end].transpose(2, 0, 1)
    return jets


def location(index, density_indices, temperature_indices, logq, logt):
    density_index, temperature_index = int(density_indices[index]), int(temperature_indices[index])
    return {"cell": [density_index, temperature_index],
            "logQ_bounds": [float(logq[density_index]), float(logq[density_index+1])],
            "logT_bounds": [float(logt[temperature_index]), float(logt[temperature_index+1])],
            "logRho_bounds": [float(logq[density_index]+2*logt[temperature_index]-12),
                              float(logq[density_index+1]+2*logt[temperature_index+1]-12)]}


def make_witness(field, density_index, temperature_index, nodes, complete, logq, logt):
    def witness(density_phase, temperature_phase):
        sample_q = np.array([logq[density_index] + (logq[density_index+1]-logq[density_index])*density_phase])
        sample_t = np.array([logt[temperature_index] + (logt[temperature_index+1]-logt[temperature_index])*temperature_phase])
        sample_q = np.clip(sample_q, logq[density_index], logq[density_index+1])
        sample_t = np.clip(sample_t, logt[temperature_index], logt[temperature_index+1])
        native = evaluate(nodes, complete, logq, logt, sample_q, sample_t,
                          cell_indices=(np.array([density_index]), np.array([temperature_index])))
        state = thermodynamics(native, sample_q, sample_t)
        value = state["Cv" if field == "Cv_numerator" else "P_rho"][0]
        if np.isfinite(value) and value < 0:
            return {"logQ": float(sample_q[0]), "logT": float(sample_t[0]),
                    "logRho": float(sample_q[0]+2*sample_t[0]-12),
                    "actual_density_phase": float((sample_q[0]-logq[density_index])/(logq[density_index+1]-logq[density_index])),
                    "actual_temperature_phase": float((sample_t[0]-logt[temperature_index])/(logt[temperature_index+1]-logt[temperature_index])),
                    "core_field": "Cv" if field == "Cv_numerator" else "P_rho",
                    "core_value": float(value), "core_Cv": float(state["Cv"][0]),
                    "core_P_rho": float(state["P_rho"][0]),
                    "core_Cv_numerator": float(state["Cv"][0]*state["T"][0])}
        return None
    return witness


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--max-depth", type=int, choices=range(5), default=4)
    parser.add_argument("--max-nodes", type=int, default=341)
    parser.add_argument("--batch-size", type=int, default=2048)
    parser.add_argument("--max-cells", type=int, default=150000)
    arguments = parser.parse_args()
    if not 1 <= arguments.batch_size <= 4096 or not 1 <= arguments.max_nodes <= 341:
        parser.error("Require batch size 1..4096 and node budget 1..341")
    if not 1 <= arguments.max_cells <= 150000:
        parser.error("Require cell budget 1..150000")
    output = arguments.output.resolve()
    if output != OUTPUT.resolve() and OUTPUT.resolve() not in output.parents:
        parser.error("Outputs must stay inside the experiment024 directory")
    started = time.monotonic()
    started_at = datetime.datetime.now().astimezone().isoformat()
    before = {str(path.relative_to(WORKSPACE)): checksum(path)
              for path in (TABLE, BIQUINTIC / "freeeos_potential.py")}
    if checksum(TABLE) != TABLE_SHA:
        raise ValueError("Unexpected source SHA: refuse silently changed scientific input")
    header, logq, logt, rows = read_table(TABLE)
    low, source_valid, state = source_constraints(rows, logq, logt, split_radiation=False)
    du, dtau = (logq[1]-logq[0])*LN10, (logt[1]-logt[0])*LN10
    nodes, complete, _ = construct_twists(low, source_valid, du, dtau, radius=2, route="mean")
    density_indices, temperature_indices = np.nonzero(valid_cells(complete))
    count = len(density_indices)
    if count > arguments.max_cells:
        raise ValueError(f"{count} cells exceed budget {arguments.max_cells}; no partial run mislabeled complete")
    print(f"{count} complete total-potential cells; bounded depth {arguments.max_depth}", flush=True)
    ratio = CRAD * state["T"]**4 / (3*state["Pgas"])
    full_own_grid, material_grid, branches, ratio_max = domain_masks(logq, logt, ratio)
    full_own = full_own_grid[density_indices, temperature_indices]
    material = material_grid[density_indices, temperature_indices]
    subsets = {"all_complete": np.ones(count, dtype=bool), "full_one_branch_own_use": full_own,
               "material_full_one_branch_own_use": material & full_own,
               "material_not_full_one_branch_own_use": material & ~full_own,
               "not_full_one_branch_own_use": ~full_own}
    corner_own = ((logt[None, :] >= 3.2) & ((logt[None, :] >= 6.5) | (logq[:, None] <= -.5)
                  | (logq[:, None] + 2*logt[None, :] - 12 <= -6.5)))
    any_corner_own = np.logical_or.reduce([corner_own[:-1, :-1], corner_own[1:, :-1],
                                          corner_own[:-1, 1:], corner_own[1:, 1:]])
    all_corner_own = valid_cells(corner_own)
    archive = {"density_cell_index": density_indices, "temperature_cell_index": temperature_indices,
               "full_one_branch_own_use": full_own, "material_corner_screen": material,
               "maximum_corner_Prad_over_Pgas": ratio_max[density_indices, temperature_indices]}
    records, witnesses, unresolved = {}, {field: [] for field in FIELDS}, {field: [] for field in FIELDS}
    reproduction_indices = np.linspace(0, count-1, min(128, count), dtype=int)
    agreement = {field: {"sample_count": 0, "maximum_absolute_difference": 0.,
                         "maximum_difference_over_margin": 0., "outside_margin_count": 0,
                         "negative_core_samples": 0, "zero_core_samples": 0,
                         "nonfinite_core_samples": 0, "worst": None} for field in FIELDS}
    reproduction = {"cell_count": 0, "jet_count": 0, "maximum_scaled_error": 0.,
                    "maximum_core_scaled_error": 0.}
    for field in FIELDS:
        records[field] = {name: np.zeros(count, dtype=dtype) for name, dtype in (
            ("status", np.int8), ("root_minimum", DTYPE), ("root_maximum", DTYPE),
            ("roundoff_margin", DTYPE), ("assembly_margin", DTYPE), ("core_margin", DTYPE),
            ("coefficient_absolute_scale", DTYPE), ("final_lower_bound", DTYPE),
            ("minimum_leaf_upper_bound", DTYPE), ("visited_nodes", np.int16),
            ("maximum_depth", np.int8), ("positive_leaves", np.int16),
            ("unresolved_leaves", np.int16), ("witness_core_queries", np.int16))}
    for start in range(0, count, arguments.batch_size):
        stop = min(start+arguments.batch_size, count)
        density_batch, temperature_batch = density_indices[start:stop], temperature_indices[start:stop]
        jets = collect_jets(nodes, density_batch, temperature_batch)
        power, magnitude = cell_power(jets, du, dtau)
        powers, scales = sign_powers(power, magnitude, du, dtau)
        density_width = logq[density_batch+1].astype(DTYPE)-logq[density_batch].astype(DTYPE)
        temperature_width = logt[temperature_batch+1].astype(DTYPE)-logt[temperature_batch].astype(DTYPE)
        density_stretch = density_width / DTYPE(logq[1]-logq[0])
        temperature_stretch = temperature_width / DTYPE(logt[1]-logt[0])
        powers = {field: normalize_closed_cell(values, density_stretch, temperature_stretch) for field, values in powers.items()}
        scales = {field: normalize_closed_cell(values, density_stretch, temperature_stretch) for field, values in scales.items()}
        bernsteins, margins = {}, {}
        for field in FIELDS:
            coefficients, margin, assembly_margin, core_margin, scale = bernstein_bounds(powers[field], scales[field])
            if not np.all(np.isfinite(coefficients)) or not np.all(np.isfinite(margin)):
                raise AssertionError("Nonfinite polynomial bounds; refuse certificate")
            bernsteins[field], margins[field] = coefficients, margin
            record = records[field]
            for name, values in (("root_minimum", coefficients.min(axis=(1, 2))),
                                 ("root_maximum", coefficients.max(axis=(1, 2))),
                                 ("roundoff_margin", margin), ("assembly_margin", assembly_margin),
                                 ("core_margin", core_margin), ("coefficient_absolute_scale", scale)):
                record[name][start:stop] = values
            root_positive = coefficients.min(axis=(1, 2)) > margin
            record["status"][start:stop] = root_positive.astype(np.int8)
            record["final_lower_bound"][start:stop] = coefficients.min(axis=(1, 2))-margin
            record["minimum_leaf_upper_bound"][start:stop] = coefficients.max(axis=(1, 2))+margin
            record["visited_nodes"][start:stop] = 1
            record["positive_leaves"][start:stop] = root_positive.astype(np.int16)
            for local_index in np.flatnonzero(~root_positive):
                index = start+local_index
                callback = make_witness(field, int(density_batch[local_index]), int(temperature_batch[local_index]),
                                        nodes, complete, logq, logt)
                result = adaptive_bound(coefficients[local_index], margin[local_index], callback,
                                        arguments.max_depth, arguments.max_nodes)
                record["status"][index] = STATUS[result["status"]]
                for name in ("visited_nodes", "maximum_depth", "positive_leaves", "unresolved_leaves",
                             "witness_core_queries", "minimum_leaf_upper_bound"):
                    record[name][index] = result[name]
                record["final_lower_bound"][index] = result["lower_bound"]
                evidence = {**location(index, density_indices, temperature_indices, logq, logt),
                            "material_corner_screen": bool(material[index]), "full_one_branch_own_use": bool(full_own[index]),
                            "root_minimum": float(record["root_minimum"][index]),
                            "roundoff_margin": float(margin[local_index]), **result}
                if result["status"] == "negative_witness":
                    witnesses[field].append(evidence)
                elif result["status"] == "unresolved":
                    unresolved[field].append(evidence)
        for density_phase, temperature_phase in PHASES:
            sample_q = np.asarray(logq[density_batch] + density_width*density_phase, dtype=float)
            sample_t = np.asarray(logt[temperature_batch] + temperature_width*temperature_phase, dtype=float)
            actual_density_phase = (sample_q.astype(DTYPE)-logq[density_batch].astype(DTYPE))/density_width
            actual_temperature_phase = (sample_t.astype(DTYPE)-logt[temperature_batch].astype(DTYPE))/temperature_width
            native = evaluate(nodes, complete, logq, logt, sample_q, sample_t,
                              cell_indices=(density_batch, temperature_batch))
            core_state = thermodynamics(native, sample_q, sample_t)
            core_values = {"Cv_numerator": core_state["Cv"]*core_state["T"], "P_rho": core_state["P_rho"]}
            for field in FIELDS:
                values = polynomial_value(powers[field], actual_density_phase, actual_temperature_phase)
                difference = abs(values-core_values[field])
                normalized = difference / margins[field]
                summary = agreement[field]
                summary["sample_count"] += len(sample_q)
                summary["outside_margin_count"] += int(np.count_nonzero(normalized > 1))
                summary["negative_core_samples"] += int(np.count_nonzero(core_values[field] < 0))
                summary["zero_core_samples"] += int(np.count_nonzero(core_values[field] == 0))
                summary["nonfinite_core_samples"] += int(np.count_nonzero(~np.isfinite(core_values[field])))
                summary["maximum_absolute_difference"] = max(summary["maximum_absolute_difference"], float(difference.max()))
                worst = int(np.argmax(normalized))
                if normalized[worst] > summary["maximum_difference_over_margin"]:
                    summary["maximum_difference_over_margin"] = float(normalized[worst])
                    summary["worst"] = {**location(start+worst, density_indices, temperature_indices, logq, logt),
                                        "logQ": float(sample_q[worst]), "logT": float(sample_t[worst]),
                                        "polynomial_value": float(values[worst]), "core_value": float(core_values[field][worst]),
                                        "absolute_difference": float(difference[worst]), "margin": float(margins[field][worst])}
                certified = records[field]["status"][start:stop] == STATUS["positive"]
                if np.any(certified & (core_values[field] <= 0)):
                    raise AssertionError("Certified polynomial disagrees in sign with sampled core")
        selected = reproduction_indices[(reproduction_indices >= start) & (reproduction_indices < stop)]-start
        if selected.size:
            reproduction["cell_count"] += len(selected)
            for density_end in range(2):
                for temperature_end in range(2):
                    sample_q = logq[density_batch[selected]+density_end]
                    sample_t = logt[temperature_batch[selected]+temperature_end]
                    core = evaluate(nodes, complete, logq, logt, sample_q, sample_t,
                                    orders=tuple((density_order, temperature_order) for density_order in range(3)
                                                 for temperature_order in range(3)),
                                    cell_indices=(density_batch[selected], temperature_batch[selected]))
                    density_phase = (sample_q-logq[density_batch[selected]])/(logq[1]-logq[0])
                    temperature_phase = (sample_t-logt[temperature_batch[selected]])/(logt[1]-logt[0])
                    for density_order in range(3):
                        for temperature_order in range(3):
                            derived = derivative(power[selected], density_order, temperature_order, du, dtau)
                            values = polynomial_value(derived, density_phase, temperature_phase)
                            expected = jets[selected, density_end, temperature_end, density_order, temperature_order]
                            scale = np.maximum(1, derivative(magnitude[selected], density_order, temperature_order, du, dtau).sum(axis=(1, 2)))
                            reproduction["jet_count"] += len(selected)
                            reproduction["maximum_scaled_error"] = max(reproduction["maximum_scaled_error"], float(np.max(abs(values-expected)/scale)))
                            reproduction["maximum_core_scaled_error"] = max(reproduction["maximum_core_scaled_error"], float(np.max(abs(core[density_order, temperature_order]-expected)/scale)))
        if start % (arguments.batch_size*10) == 0:
            print(f"{stop}/{count} cells bounded and sample-checked", flush=True)
    summaries = {}
    for name, mask in subsets.items():
        summaries[name] = {"cell_count": int(mask.sum()), "fields": {}}
        both_positive = mask.copy()
        for field, record in records.items():
            summary = {status: int(np.count_nonzero(mask & (record["status"] == status_value)))
                       for status, status_value in STATUS.items()}
            summary["root_positive"] = int(np.count_nonzero(mask & (record["root_minimum"] > record["roundoff_margin"])))
            summary["adaptive_positive"] = summary["positive"]-summary["root_positive"]
            summary["visited_nodes"] = int(record["visited_nodes"][mask].sum())
            summary["unresolved_leaves"] = int(record["unresolved_leaves"][mask].sum())
            summary["maximum_depth"] = int(record["maximum_depth"][mask].max()) if mask.any() else 0
            summary["extrema"] = {}
            for metric, maximize in (("root_minimum", False), ("final_lower_bound", False),
                                     ("roundoff_margin", True), ("assembly_margin", True), ("core_margin", True)):
                if mask.any():
                    index = int(np.argmax(np.where(mask, record[metric], -np.inf)) if maximize
                                else np.argmin(np.where(mask, record[metric], np.inf)))
                    summary["extrema"][metric] = {"value": float(record[metric][index]),
                                                  **location(index, density_indices, temperature_indices, logq, logt)}
            summaries[name]["fields"][field] = summary
            both_positive &= record["status"] == STATUS["positive"]
        summaries[name]["both_numerically_certified_positive"] = int(both_positive.sum())
    after = {path: checksum(WORKSPACE / path) for path in before}
    if before != after:
        raise AssertionError("Source/core changed during audit")
    reproduction["passed"] = bool(reproduction["maximum_scaled_error"] < SAFETY*np.finfo(np.float64).eps
                                 and reproduction["maximum_core_scaled_error"] < SAFETY*np.finfo(np.float64).eps)
    for field, record in records.items():
        archive.update({f"{field}__{name}": values for name, values in record.items()})
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output / "cell_bounds.npz", **archive)
    (output / "negative_witnesses.json").write_text(json.dumps(witnesses, indent=2, allow_nan=False)+"\n")
    (output / "unresolved_cells.json").write_text(json.dumps(unresolved, indent=2, allow_nan=False)+"\n")
    result = {
        "started_at": started_at, "elapsed_seconds": time.monotonic()-started,
        "command": sys.argv, "python": sys.executable, "numpy": np.__version__,
        "candidate": "radius2_mean_total", "table_header": header.tolist(),
        "source_sha256": before, "source_unchanged": before == after,
        "sidecar_sha256": {str(path.relative_to(WORKSPACE)): checksum(path) for path in SIDECAR.glob("*.py")},
        "units": "Natural logarithms u=ln(Q), tau=ln(T); normalized x,y; jets/spacings exactly as core float64 inputs",
        "closed_cell_normalization": "Core uses one global dq/dt for H basis; rescale each sign polynomial by actual endpoint width/global spacing so Bernstein [0,1]^2 covers the entire stored closed cell, including binary grid-width variations",
        "formulas": {"Cv_numerator": "f_tau - 2 f_u - f_tautau + 4 f_utau - 4 f_uu = T Cv",
                     "P_rho": "f_u + f_uu"},
        "method": "6x6 tensor power coefficients from all nine endpoint jets; degree-5 Bernstein elevation; strict min>margin; dyadic de Casteljau covering",
        "floating_proof_claim": False,
        "margin_policy": {"safety_factor": SAFETY, "longdouble_epsilon": float(np.finfo(DTYPE).eps),
                          "float64_epsilon": float(np.finfo(np.float64).eps), "longdouble_mantissa_bits": np.finfo(DTYPE).nmant,
                          "formula": "4096*(eps_longdouble + eps_float64)*max Bernstein of absolute assembly/derivative contributions",
                          "subdivision_addition": "Inherited margin plus 4096*eps_longdouble*max(abs(parent Bernstein)) per level",
                          "scope": "Conservative numerical guards, NOT directed-rounded intervals or formal floating-point error bounds; no source jet uncertainty model"},
        "budgets": {"max_depth": arguments.max_depth, "max_nodes_per_cell_field": arguments.max_nodes,
                    "batch_size": arguments.batch_size, "max_cells": arguments.max_cells,
                    "max_live_dfs_nodes": 1+3*arguments.max_depth,
                    "sample_phases": PHASES},
        "domain": {"complete_cells": count, "source_grid_shape": list(rows.shape[:2]),
                   "source_valid_nodes": int(source_valid.sum()), "complete_nodes": int(complete.sum()),
                   "one_branch_definition": "Entire closed cell: tmin>=3.2 and (tmin>=6.5 OR qmax<=-0.5 OR maxlogrho<=-6.5), no boundary tolerance",
                   "material_definition": "Four source corners Prad/Pgas<1e-3, entire q in [-8,5], t in [3.2,7.8]; source-corner screen, NOT a continuous material-fraction bound",
                   "full_branch_counts_overlapping": {name: int(mask[density_indices, temperature_indices].sum()) for name, mask in branches.items()},
                   "valid_partial_or_union_only_own_use_cells": int(np.count_nonzero(any_corner_own[density_indices, temperature_indices] & ~full_own)),
                   "all_corner_own_but_no_single_branch_cells": int(np.count_nonzero(all_corner_own[density_indices, temperature_indices] & ~full_own)),
                   "not_full_own_use_policy": "Reported separately, never counted in primary material/full-own-use claim"},
        "counts": summaries, "sample_agreement_against_core": agreement, "Hermite_reproduction": reproduction,
        "status_codes": STATUS,
        "calibration_passed": reproduction["passed"] and all(item["outside_margin_count"] == 0 and item["nonfinite_core_samples"] == 0 for item in agreement.values()),
        "limitations": ["Only the implemented total-potential cell polynomial, not actual source truth, global stability, runtime or blended EOS",
                        "Analytic radiation split excluded: exponential radiation is not a biquintic polynomial",
                        "No Pgas/E/S full-cell certificate, no intercell smoothness or composition claim",
                        "Nonpositive Bernstein minima are ambiguous, not physical failures; only evaluated core-negative witnesses are failures",
                        "Finite depth/budget unresolved cells remain explicitly unresolved; no sample-only pass"],
        "artifacts_sha256": {path.name: checksum(path) for path in (output / "cell_bounds.npz", output / "negative_witnesses.json", output / "unresolved_cells.json")},
    }
    (output / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    primary = summaries["material_full_one_branch_own_use"]
    lines = ["# Experiment024: total-potential Bernstein positivity bounds", "",
             f"Completed in {result['elapsed_seconds']:.1f} seconds; {count} complete cells.",
             "Numerical certificates use conservative margins, not rigorous floating-point intervals.", "",
             "| Subset | Cells | Cv positive / witness / unresolved | P_rho positive / witness / unresolved | Both positive |",
             "| --- | ---: | --- | --- | ---: |"]
    for name, summary in summaries.items():
        entries = [" / ".join(str(summary["fields"][field][key]) for key in ("positive", "negative_witness", "unresolved")) for field in FIELDS]
        lines.append(f"| {name} | {summary['cell_count']} | {entries[0]} | {entries[1]} | {summary['both_numerically_certified_positive']} |")
    lines += ["", f"Primary: {primary['both_numerically_certified_positive']}/{primary['cell_count']} material-screened cells entirely in one own-use branch have both numerical sign certificates.",
              "Cv uses T Cv = f_tau-2f_u-f_tautau+4f_utau-4f_uu; P_rho uses f_u+f_uu.",
              "Material means only the archived four-corner radiation fraction screen, not a cellwide radiation-fraction theorem.",
              f"Sample calibration passes: {result['calibration_passed']}; four actual-coordinate phases per cell, not the certificate mechanism.",
              "Read results.json for exact worst bounds, margins, source SHA and provenance; cell_bounds.npz retains every cell/field status.",
              "negative_witnesses.json contains actual core-evaluated coordinates; unresolved_cells.json retains all budget/depth ambiguities.",
              "No source accuracy, global/blended stability, runtime/cache change, analytic radiation-split certificate or Pgas/E/S certificate.",
              "Lead integration: use the primary row with numerical scope qualifiers; keep all-complete and partial-domain rows separate."]
    (output / "SUMMARY.md").write_text("\n".join(lines)+"\n")
    print("\n".join(lines), flush=True)
    if not result["calibration_passed"]:
        raise SystemExit("Calibration failed; do not use numerical certificates")


if __name__ == "__main__":
    main()
