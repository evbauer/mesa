"""Matched saved-query departures and coarse-only heldout C3 feasibility audit."""

import argparse
import json
import platform
import sys
import time
import unittest
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

from c3_potential import (
    ESTIMATED_ORDERS, JET_ORDERS, ROUTE_GRAPH, construct_jets, evaluate as evaluate_c3,
)
from freeeos_potential import (
    CRAD, EVAL_ORDERS, LN10, LOW_ORDERS, checksum, construct_twists,
    evaluate as evaluate_biquintic, query_cells, read_table, source_constraints,
    thermodynamic_partials, valid_cells,
)
from run_audit import own_use
from validate_table import material_nodes, table_state


ROOT = Path(__file__).resolve().parents[5]
HERE = Path(__file__).resolve().parent
FIELDS = ("P", "Pgas", "E", "S", "Cv", "Cp", "chiRho", "chiT", "gamma1", "gamma3", "grad_ad")
PARTIAL_FIELDS = tuple(field + suffix for field in FIELDS for suffix in ("_x", "_t"))
POSITIVE_FIELDS = ("P", "Pgas", "Cv", "P_rho", "chiRho", "E", "S")


def summary(values, sample_q=None, sample_t=None):
    values = np.asarray(values).ravel()
    finite = np.isfinite(values)
    result = {"count": int(values.size), "nonfinite": int((~finite).sum())}
    if finite.any():
        result.update(zip(("median", "p90", "p99", "max"),
                          map(float, np.quantile(abs(values[finite]), (.5, .9, .99, 1)))))
        if sample_q is not None:
            index = np.flatnonzero(finite)[np.argmax(abs(values[finite]))]
            result["worst"] = {"logQ": float(sample_q[index]), "logT": float(sample_t[index]),
                               "signed_scaled_difference": float(values[index])}
    return result


def positivity(state, sample_q, sample_t):
    result = {}
    for field in POSITIVE_FIELDS:
        values = state[field]
        bad = (~np.isfinite(values)) | (values <= 0)
        result[field] = {"count": int(bad.sum()), "nonfinite": int((~np.isfinite(values)).sum()),
                         "examples": [{"logQ": float(sample_q[index]), "logT": float(sample_t[index]),
                                       "value": float(values[index]) if np.isfinite(values[index]) else None}
                                      for index in np.flatnonzero(bad)[:8]]}
    return result


def build_models(logq, logt, rows):
    models, constraints = {}, {}
    du, dtau = (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10
    for mode, split in (("total", False), ("analytic_split", True)):
        low, valid, _ = source_constraints(rows, logq, logt, split_radiation=split)
        constraints[mode] = (low, valid)
        for radius in (2, 3):
            nodes, complete, diagnostics = construct_jets(low, valid, du, dtau, radius)
            models[f"c3_r{radius}_{mode}"] = {"nodes": nodes, "complete": complete,
                "evaluate": evaluate_c3, "split": split, "diagnostics": diagnostics}
            nodes, complete, routes = construct_twists(low, valid, du, dtau, radius, "mean")
            models[f"biquintic_r{radius}_{mode}"] = {"nodes": nodes, "complete": complete,
                "evaluate": evaluate_biquintic, "split": split, "diagnostics": routes}
    common_cells = np.ones((len(logq) - 1, len(logt) - 1), bool)
    for model in models.values():
        common_cells &= valid_cells(model["complete"])
    return models, constraints, common_cells


def native_at(model, logq, logt, sample_q, sample_t, orders=EVAL_ORDERS, cell_indices=None):
    return model["evaluate"](model["nodes"], model["complete"], logq, logt,
        sample_q, sample_t, orders=orders, cell_indices=cell_indices, add_radiation=model["split"])


def state_at(model, logq, logt, sample_q, sample_t):
    return thermodynamic_partials(native_at(model, logq, logt, sample_q, sample_t), sample_q, sample_t)


def pairs(models):
    return [(name, name.replace("c3_", "biquintic_")) for name in models if name.startswith("c3_")]


def mask_products(models, constraints, prefix=""):
    products = {prefix + mode + "__source_valid": valid for mode, (_, valid) in constraints.items()}
    for name, model in models.items():
        products[prefix + name + "__complete"] = model["complete"]
        if name.startswith("c3_"):
            for orders, mask in model["diagnostics"]["jet_masks"].items():
                products[prefix + name + "__jet" + "".join(map(str, orders))] = mask
    return products


def estimator_diagnostics(models):
    result = {}
    for name, model in models.items():
        if not name.startswith("c3_"):
            continue
        current = {"complete_nodes": int(model["complete"].sum()),
                   "complete_cells": int(valid_cells(model["complete"]).sum()), "jets": {}}
        for orders in JET_ORDERS:
            current["jets"][str(orders)] = {"accepted_nodes": int(model["diagnostics"]["jet_masks"][orders].sum()),
                                          "supplied": orders in LOW_ORDERS}
        for orders, route in model["diagnostics"]["routes"].items():
            item = current["jets"][str(orders)]
            item["routes"] = [{"parent": list(parent), "axis": axis} for parent, axis in route["parents"]]
            item["combination"] = "arithmetic mean requiring ALL route masks" if len(route["values"]) == 2 else "single estimated pure-third route"
            if len(route["values"]) == 2:
                first, second = (value[route["valid"]] for value in route["values"])
                scale = np.maximum.reduce([abs(first), abs(second), np.ones(first.size)])
                item["route_disagreement_relative_or_floor1"] = summary((first - second) / scale)
        result[name] = current
    return result


def coverage(models, logq, logt, sample_q=None, sample_t=None, eligible=None, material=None):
    result = {}
    for name, model in models.items():
        cells = valid_cells(model["complete"])
        current = {"complete_nodes": int(model["complete"].sum()), "complete_cells": int(cells.sum())}
        if sample_q is not None:
            first, second, inside = query_cells(logq, logt, sample_q, sample_t)
            accepted = inside & eligible & own_use(sample_q, sample_t) & cells[first, second]
            current["own_use_queries_before_common_intersection"] = int(accepted.sum())
            current["material_own_use_queries_before_common_intersection"] = int((accepted & material).sum())
        result[name] = current
    return result


def matched_departures(models, logq, logt, sample_q, sample_t, subsets, chunk_size):
    result = {}
    fields = (*FIELDS, *PARTIAL_FIELDS)
    for candidate_name, reference_name in pairs(models):
        print(f"saved-query comparison {candidate_name}: {len(sample_q):,} queries", flush=True)
        differences = np.empty((len(fields), len(sample_q)))
        failure_masks = {name: {field: np.zeros(len(sample_q), bool) for field in POSITIVE_FIELDS}
                         for name in (candidate_name, reference_name)}
        nonfinite_masks = {name: {field: np.zeros(len(sample_q), bool) for field in POSITIVE_FIELDS}
                           for name in (candidate_name, reference_name)}
        failure_examples = {name: {subset: {field: [] for field in POSITIVE_FIELDS} for subset in subsets}
                            for name in (candidate_name, reference_name)}
        for start in range(0, len(sample_q), chunk_size):
            stop = min(start + chunk_size, len(sample_q))
            query_q, query_t = sample_q[start:stop], sample_t[start:stop]
            candidate = state_at(models[candidate_name], logq, logt, query_q, query_t)
            reference = state_at(models[reference_name], logq, logt, query_q, query_t)
            with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                for index, field in enumerate(fields):
                    scale = np.maximum(abs(reference[field]), 1e-300)
                    if field in PARTIAL_FIELDS:
                        scale = np.maximum.reduce([abs(reference[field]), abs(reference[field[:-2]]), np.ones(stop - start)])
                    differences[index, start:stop] = (candidate[field] - reference[field]) / scale
            for name, state in ((candidate_name, candidate), (reference_name, reference)):
                for field in POSITIVE_FIELDS:
                    bad = (~np.isfinite(state[field])) | (state[field] <= 0)
                    failure_masks[name][field][start:stop] = bad
                    nonfinite_masks[name][field][start:stop] = ~np.isfinite(state[field])
                    for subset, selected in subsets.items():
                        examples = failure_examples[name][subset][field]
                        for index in np.flatnonzero(bad & selected[start:stop])[:max(0, 8 - len(examples))]:
                            value = state[field][index]
                            examples.append({"logQ": float(query_q[index]), "logT": float(query_t[index]),
                                             "value": float(value) if np.isfinite(value) else None})
        current = {"reference": reference_name, "subsets": {}}
        for subset, selected in subsets.items():
            current["subsets"][subset] = {"count": int(selected.sum()),
                "departures_not_source_accuracy": {field: summary(differences[index, selected], sample_q[selected], sample_t[selected])
                                                  for index, field in enumerate(fields)},
                "positivity": {name: {field: {"count": int(failure_masks[name][field][selected].sum()),
                    "nonfinite": int(nonfinite_masks[name][field][selected].sum()),
                    "examples": failure_examples[name][subset][field]} for field in POSITIVE_FIELDS}
                    for name in (candidate_name, reference_name)}}
        result[candidate_name] = current
    return result


def face_audit(models, logq, logt, rows, common_cells, max_faces):
    generator = np.random.default_rng(20261031)
    result, products = {}, {}
    material_cells = valid_cells(material_nodes(rows, logq, logt))
    for axis, direction in ((0, "u"), (1, "tau")):
        adjacent = common_cells[:-1] & common_cells[1:] if axis == 0 else common_cells[:, :-1] & common_cells[:, 1:]
        first, second = np.nonzero(adjacent)
        all_q = logq[first + (axis == 0)]
        all_t = logt[second + (axis == 1)]
        selected = own_use(all_q, all_t)
        first, second = first[selected], second[selected]
        choice = generator.choice(len(first), min(max_faces, len(first)), replace=False)
        first, second = first[choice], second[choice]
        right_first, right_second = first + (axis == 0), second + (axis == 1)
        material = material_cells[first, second] & material_cells[right_first, right_second]
        for mode, fraction in (("nodal_tangent", np.zeros(len(first))),
                               ("offnode_tangent", generator.uniform(.1, .9, len(first)))):
            sample_q = logq[first + 1] if axis == 0 else logq[first] + fraction * (logq[1] - logq[0])
            sample_t = logt[second] + fraction * (logt[1] - logt[0]) if axis == 0 else logt[second + 1]
            accepted = own_use(sample_q, sample_t)
            query_q, query_t = sample_q[accepted], sample_t[accepted]
            left_indices = first[accepted], second[accepted]
            right_indices = right_first[accepted], right_second[accepted]
            key = direction + "__" + mode
            products[key + "__logQ"], products[key + "__logT"] = query_q, query_t
            products[key + "__left_i"], products[key + "__left_j"] = left_indices
            products[key + "__right_i"], products[key + "__right_j"] = right_indices
            products[key + "__material"] = material[accepted]
            result[key] = {}
            for name, model in models.items():
                first_native = native_at(model, logq, logt, query_q, query_t, JET_ORDERS, left_indices)
                second_native = native_at(model, logq, logt, query_q, query_t, JET_ORDERS, right_indices)
                first_state = thermodynamic_partials(first_native, query_q, query_t)
                second_state = thermodynamic_partials(second_native, query_q, query_t)
                current = {}
                for subset, selected in (("own_use", np.ones(query_q.size, bool)), ("material_own_use", material[accepted])):
                    native_metrics = {}
                    for orders in JET_ORDERS:
                        scale = np.maximum.reduce([abs(first_native[orders]), abs(second_native[orders]), np.ones(query_q.size)])
                        native_metrics[str(orders)] = summary(((first_native[orders] - second_native[orders]) / scale)[selected],
                                                            query_q[selected], query_t[selected])
                    state_metrics = {}
                    for field in (*FIELDS, *PARTIAL_FIELDS):
                        scale = np.maximum.reduce([abs(first_state[field]), abs(second_state[field]), np.ones(query_q.size)])
                        if field in PARTIAL_FIELDS:
                            scale = np.maximum.reduce([scale, abs(first_state[field[:-2]]), abs(second_state[field[:-2]])])
                        state_metrics[field] = summary(((first_state[field] - second_state[field]) / scale)[selected],
                                                       query_q[selected], query_t[selected])
                    current[subset] = {"count": int(selected.sum()), "tensor_partial_face_jumps": native_metrics,
                                       "value_and_response_jacobian_face_jumps": state_metrics}
                result[key][name] = current
    return result, products


def nodal_audit(models, logq, logt, common_cells):
    first, second = np.nonzero(common_cells)
    choice = np.random.default_rng(11).choice(len(first), min(2000, len(first)), replace=False)
    first, second = first[choice], second[choice]
    sample_q, sample_t = logq[first], logt[second]
    result = {}
    for name, model in models.items():
        if not name.startswith("c3_"):
            continue
        evaluated = model["evaluate"](model["nodes"], model["complete"], logq, logt,
            sample_q, sample_t, JET_ORDERS, cell_indices=(first, second))
        result[name] = {str(orders): summary((evaluated[orders] - model["nodes"][orders][first, second])
            / np.maximum(abs(model["nodes"][orders][first, second]), 1)) for orders in JET_ORDERS}
    return result


def heldout(logq, logt, rows, full_valid):
    coarse_q, coarse_t, coarse_rows = logq[::2], logt[::2], rows[::2, ::2].copy()
    models, constraints, common_cells = build_models(coarse_q, coarse_t, coarse_rows)
    first, second = np.nonzero(full_valid)
    sample_q, sample_t = logq[first], logt[second]
    coarse_first, coarse_second, inside = query_cells(coarse_q, coarse_t, sample_q, sample_t)
    selected = inside & (first % 2 == 1) & (second % 2 == 1) & own_use(sample_q, sample_t)
    selected &= common_cells[coarse_first, coarse_second]
    first, second = first[selected], second[selected]
    sample_q, sample_t = sample_q[selected], sample_t[selected]
    coarse_first, coarse_second = coarse_first[selected], coarse_second[selected]
    material = material_nodes(rows, logq, logt)[first, second]
    material &= valid_cells(material_nodes(coarse_rows, coarse_q, coarse_t))[coarse_first, coarse_second]
    reference = table_state(rows[first, second])
    reference["P"] = reference["Pgas"] + CRAD * 10.0**(4 * sample_t) / 3
    result = {"stride": 2, "coarse_only_source": True, "withheld_in_both_axes": True,
              "coverage": coverage(models, coarse_q, coarse_t),
              "estimator": estimator_diagnostics(models), "models": {}}
    for name, model in models.items():
        state = state_at(model, coarse_q, coarse_t, sample_q, sample_t)
        current = {}
        for subset, accepted in (("own_use", np.ones(len(sample_q), bool)), ("material_own_use", material)):
            current[subset] = {"count": int(accepted.sum()),
                "relative_printed_node_fidelity": {field: summary(((state[field] - reference[field])
                    / np.maximum(abs(reference[field]), 1e-300))[accepted], sample_q[accepted], sample_t[accepted])
                    for field in FIELDS},
                "positivity": positivity({field: state[field][accepted] for field in POSITIVE_FIELDS}, sample_q[accepted], sample_t[accepted])}
        result["models"][name] = current
    products = mask_products(models, constraints, "coarse__")
    products.update(heldout_i=first, heldout_j=second, heldout_logQ=sample_q, heldout_logT=sample_t,
                    heldout_material=material, coarse_common_cells=common_cells)
    return result, products


def analytic_convergence():
    generator = np.random.default_rng(11)
    sample_q, sample_t = generator.uniform(-.25, .25, (2, 101)) / LN10
    result = {"function": "exp(1.3 u + 0.9 tau)", "queries": 101, "grids": []}
    for count in (33, 49, 65):
        logq = logt = np.linspace(-2, 2, count) / LN10
        first, second = np.meshgrid(logq * LN10, logt * LN10, indexing="ij")
        exact = np.empty((4, 4, count, count))
        for orders in JET_ORDERS:
            exact[orders] = 1.3**orders[0] * .9**orders[1] * np.exp(1.3 * first + .9 * second)
        low = np.full((3, 3, count, count), np.nan)
        for orders in LOW_ORDERS:
            low[orders] = exact[orders]
        valid = np.ones((count, count), bool)
        spacing = (logq[1] - logq[0]) * LN10
        central = (abs(first) <= .25 + 1e-12) & (abs(second) <= .25 + 1e-12)
        entry = {"nodes_per_axis": count, "natural_spacing": float(spacing), "models": {}}
        for radius in (2, 3):
            nodes, complete, _ = construct_jets(low, valid, spacing, spacing, radius)
            nodal_errors = {str(orders): summary((nodes[orders][central & complete] - exact[orders][central & complete])
                / abs(exact[orders][central & complete])) for orders in JET_ORDERS}
            for name, current_nodes, current_complete, evaluator in (
                    (f"estimated_c3_r{radius}", nodes, complete, evaluate_c3),
                    (f"biquintic_r{radius}", *construct_twists(low, valid, spacing, spacing, radius, "mean")[:2], evaluate_biquintic)):
                actual = evaluator(current_nodes, current_complete, logq, logt, sample_q, sample_t)
                metrics = {str(orders): summary((actual[orders] - 1.3**orders[0] * .9**orders[1]
                    * np.exp(1.3 * sample_q * LN10 + .9 * sample_t * LN10))
                    / (1.3**orders[0] * .9**orders[1] * np.exp(1.3 * sample_q * LN10 + .9 * sample_t * LN10)))
                    for orders in EVAL_ORDERS}
                entry["models"][name] = {"ten_interpolated_partials": metrics}
                if name.startswith("estimated_c3"):
                    entry["models"][name]["sixteen_estimated_nodal_jets"] = nodal_errors
        actual = evaluate_c3(exact, valid, logq, logt, sample_q, sample_t)
        entry["exact_jet_c3"] = {str(orders): summary((actual[orders] - 1.3**orders[0] * .9**orders[1]
            * np.exp(1.3 * sample_q * LN10 + .9 * sample_t * LN10))
            / (1.3**orders[0] * .9**orders[1] * np.exp(1.3 * sample_q * LN10 + .9 * sample_t * LN10)))
            for orders in EVAL_ORDERS}
        result["grids"].append(entry)
    for name in ("estimated_c3_r2", "estimated_c3_r3", "biquintic_r2", "biquintic_r3"):
        for previous, current in zip(result["grids"], result["grids"][1:]):
            order = {}
            for derivative in EVAL_ORDERS:
                before = previous["models"][name]["ten_interpolated_partials"][str(derivative)]["max"]
                after = current["models"][name]["ten_interpolated_partials"][str(derivative)]["max"]
                order[str(derivative)] = float(np.log(before / after) / np.log(previous["natural_spacing"] / current["natural_spacing"]))
            current["models"][name]["observed_order_from_previous"] = order
    return result


def write_summary(output, result):
    lines = ["# Offline C3 biseptic control", "", "Not a runtime or primary replacement; no new EOS tables.", "",
        "Same archived queries, all eight models' complete-cell intersection, native-only donors.",
        "Departures are against biquintic, not direct FreeEOS accuracy. Stride-2 errors use printed withheld nodes.", "",
        "| Candidate | Material saved queries | P p99 departure | Cv p99 departure | Cv_t p99 scaled departure |",
        "| --- | ---: | ---: | ---: | ---: |"]
    for name, comparison in result["saved_query_departures"].items():
        subset = comparison["subsets"]["material_own_use"]
        metrics = subset["departures_not_source_accuracy"]
        lines.append(f"| {name} | {subset['count']:,} | {metrics['P']['p99']:.3g} | {metrics['Cv']['p99']:.3g} | {metrics['Cv_t']['p99']:.3g} |")
    lines.extend(["", "| Model | Complete cells | Saved own-use coverage before common intersection |",
                  "| --- | ---: | ---: |"])
    for name, current in result["domain"]["coverage_by_model"].items():
        lines.append(f"| {name} | {current['complete_cells']:,} | {current['own_use_queries_before_common_intersection']:,} |")
    lines.extend(["", "| Model | Material stride-2 nodes | Pgas p99 fidelity | Cv p99 fidelity | Cv maximum fidelity |",
                  "| --- | ---: | ---: | ---: | ---: |"])
    for name, model in result["stride2_heldout"]["models"].items():
        subset = model["material_own_use"]
        metrics = subset["relative_printed_node_fidelity"]
        lines.append(f"| {name} | {subset['count']:,} | {metrics['Pgas']['p99']:.3g} | {metrics['Cv']['p99']:.3g} | {metrics['Cv']['max']:.3g} |")
    lines.extend(["", "| Face sample | Model | Maximum material response-Jacobian jump |", "| --- | --- | ---: |"])
    for face, models in result["exact_faces"].items():
        for name, model in models.items():
            metrics = model["material_own_use"]["value_and_response_jacobian_face_jumps"]
            maximum = max(metrics[field]["max"] for field in PARTIAL_FIELDS)
            lines.append(f"| {face} | {name} | {maximum:.3g} |")
    failed_saved = sum(data["count"] for comparison in result["saved_query_departures"].values()
        for subset in comparison["subsets"].values() for model in subset["positivity"].values() for data in model.values())
    failed_heldout = sum(data["count"] for model in result["stride2_heldout"]["models"].values()
        for subset in model.values() for data in subset["positivity"].values())
    positivity_message = ("No sampled nonfinite or nonpositive P/Pgas/Cv/P_rho/chiRho/E/S in either saved-query subset or either stride-2 subset, for all eight models."
        if failed_saved + failed_heldout == 0 else "Sampled positivity failures occur; consult per-model/subset counts and coordinates in results.json.")
    lines.extend(["", "## Positivity And Analytic Checks", "", positivity_message,
        "This is sampled-domain evidence, not a global positivity guarantee.", "",
        f"{result['tests']['run']} independent tests pass: exact rational/cardinal degree-7 tensor reproduction, degree-4/6 complete estimator, all-stage mask propagation/overflow and poisoning, coarse-only donor isolation, coordinate-exchange symmetry, exponential convergence, bounds, exact faces, and physical-coordinate response finite differences.",
        "", "| Estimated C3 | h=0.125 maximum nodal-jet relative error | h=0.0625 maximum nodal-jet relative error |",
        "| --- | ---: | ---: |"])
    for radius in (2, 3):
        name = f"estimated_c3_r{radius}"
        before = max(value["max"] for value in result["analytic_convergence"]["grids"][0]["models"][name]["sixteen_estimated_nodal_jets"].values())
        after = max(value["max"] for value in result["analytic_convergence"]["grids"][-1]["models"][name]["sixteen_estimated_nodal_jets"].values())
        lines.append(f"| radius {radius} | {before:.3g} | {after:.3g} |")
    lines.extend(["", "C3 continuity is a construction check, not evidence of better thermodynamics or source accuracy.",
        "Ten of sixteen nodal jets are estimated, including both pure thirds and orders four through six.",
        "More derivative stages cost donor coverage and amplify printed-data/FD sensitivity.",
        "Total and declared analytic-radiation-split modes are distinct controls; the split remainder is not verified pure material.",
        "In this material stride-2 test C3 Cv p99 error is roughly twice biquintic and its maximum is over twice biquintic; the smoothness gain does not justify a fidelity-superiority claim.",
        "Full distributions, positivity counts, worst coordinates, masks, fingerprints, analytic convergence, and test logs accompany this summary.", ""])
    (output / "SUMMARY.md").write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, default=ROOT / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data")
    parser.add_argument("--archive", type=Path, default=ROOT / "freeeos_report/plots/offgrid_sample_data.npz")
    parser.add_argument("--output", type=Path, default=ROOT / "freeeos_biquintic_progress/experiments/011_c3_control")
    parser.add_argument("--chunk-size", type=int, default=16384)
    parser.add_argument("--max-faces", type=int, default=2000)
    args = parser.parse_args()
    if args.chunk_size <= 0 or args.max_faces <= 0:
        parser.error("Chunk size and face count must be positive")
    started = time.monotonic()
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "tests.txt").open("w") as stream:
        suite = unittest.defaultTestLoader.discover(str(HERE), pattern="test_*.py")
        tested = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    if not tested.wasSuccessful():
        raise RuntimeError(f"C3 tests failed; see {args.output / 'tests.txt'}")
    _, logq, logt, rows = read_table(args.table)
    models, constraints, common_cells = build_models(logq, logt, rows)
    with np.load(args.archive) as archive:
        sample_q, sample_t = archive["logQ"].ravel(), archive["logT"].ravel()
        eligible, material_saved = archive["eligible"].ravel(), archive["material"].ravel()
        shape = archive["logQ"].shape
    first, second, inside = query_cells(logq, logt, sample_q, sample_t)
    selected = inside & eligible & own_use(sample_q, sample_t) & common_cells[first, second]
    material_all = material_saved & valid_cells(material_nodes(rows, logq, logt))[first, second]
    coverage_by_model = coverage(models, logq, logt, sample_q, sample_t, eligible, material_all)
    ids = np.flatnonzero(selected)
    material = material_all[selected]
    query_q, query_t = sample_q[selected], sample_t[selected]
    if not material.any():
        raise RuntimeError("No common saved material/own-use query domain")
    products = mask_products(models, constraints)
    products.update(logQ_grid=logq, logT_grid=logt, common_cells=common_cells,
        saved_flat_indices=ids, saved_logQ=query_q, saved_logT=query_t, saved_material=material,
        source_native=rows[..., 17] == 0, source_iq520=rows[..., 17] == 520,
        source_cold=np.broadcast_to(logt < 3.2 - 1e-12, rows.shape[:2]))
    result = {"status": "offline feasibility control, not runtime, primary replacement, or direct-source accuracy",
        "provenance": {"source_sha256": checksum(args.table), "saved_query_sha256": checksum(args.archive),
            "core_sha256": checksum(HERE.parent / "freeeos_potential.py"),
            "sidecar_sha256": {path.name: checksum(path) for path in sorted(HERE.glob("*.py"))},
            "python": platform.python_version(), "numpy": np.__version__, "command": sys.argv,
            "CRAD": float(CRAD), "composition": "archived solar X=.7 Z=.02"},
        "tests": {"run": tested.testsRun, "failures": len(tested.failures), "errors": len(tested.errors)},
        "domain": {"saved_shape": shape, "saved_count": int(sample_q.size), "common_cells": int(common_cells.sum()),
            "common_own_use_queries": len(query_q), "common_material_own_use_queries": int(material.sum()),
            "coverage_by_model": coverage_by_model,
            "saved_eligibility_required": True, "material": "saved material AND source-corner Prad/Pgas<1e-3, logQ[-8,5], logT<=7.8",
            "donors": "six source_constraints only; iq=0, logT>=3.2, finite source/redundant columns, density-consistent, positive source Cv/chiRho; ALL centered donors including zero-weight center and ALL redundant routes required at EVERY stage",
            "donor_selection_cuts": "material and own-use are query cuts, NOT donor cuts; donors may be outside those cuts but never outside native accepted source mask",
            "support": "up to four successive first-derivative stages from supplied second derivatives; conservative source support |du_index|+|dtau_index|<=4r with each axis<=3r; actual masks saved per jet",
            "nodal_jets": {"supplied": list(LOW_ORDERS), "estimated": list(ESTIMATED_ORDERS), "graph": {str(key): list(value) for key, value in ROUTE_GRAPH.items()}},
            "departure_scaling": "values: (C3-BQ)/abs(BQ); Jacobians: (C3-BQ)/max(1,abs(BQ response),abs(BQ partial)); x=lnrho,t=lnT; nonfinite counts retained"},
        "estimator": estimator_diagnostics(models), "analytic_convergence": analytic_convergence()}
    result["nodal_16_jet_reproduction"] = nodal_audit(models, logq, logt, common_cells)
    result["saved_query_departures"] = matched_departures(models, logq, logt, query_q, query_t,
        {"own_use": np.ones(len(query_q), bool), "material_own_use": material}, args.chunk_size)
    print("exact common faces and coarse-only stride-2 heldouts", flush=True)
    result["exact_faces"], face_products = face_audit(models, logq, logt, rows, common_cells, args.max_faces)
    result["stride2_heldout"], heldout_products = heldout(logq, logt, rows, constraints["total"][1])
    products.update(heldout_products)
    np.savez_compressed(args.output / "domains.npz", **products)
    np.savez_compressed(args.output / "face_queries.npz", **face_products)
    result["elapsed_seconds"] = time.monotonic() - started
    (args.output / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    write_summary(args.output, result)
    print(json.dumps({"tests": result["tests"], "domain": result["domain"], "elapsed_seconds": result["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
