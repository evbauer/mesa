"""Table-only potential diagnostics at the inherited dense off-node queries."""

import argparse
import json
import os
import sys
import time
import warnings
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

from freeeos_potential import (
    LN10, checksum, construct_twists, evaluate, query_cells, read_table,
    source_constraints, thermodynamic_partials, valid_cells,
)


ROOT = Path(__file__).resolve().parents[4]
DEFAULT_OUTPUT = ROOT / "freeeos_biquintic_progress/experiments/004_dense_validation"
POSITIVE_FIELDS = ("Cv", "P_rho", "Pgas", "E", "S")
COMPARISON_FIELDS = (*POSITIVE_FIELDS, "P", "chiRho", "chiT", "Cp", "gamma1", "gamma3", "grad_ad")
IDENTITY_FIELDS = ("dpe", "dse", "dsp")
SUBSETS = ("all_valid", "broad_own_use", "material_own_use", "material_valid")
VARIANTS = {"radius2_mean_total": False, "radius2_mean_analytic_split": True}


def finite_number(value):
    return float(value) if np.isfinite(value) else None


def distribution(values):
    values = np.asarray(values, dtype=float).ravel()
    finite = values[np.isfinite(values)]
    result = {"count": int(values.size), "finite": int(finite.size),
              "nonfinite": int(values.size - finite.size)}
    if not finite.size:
        return {**result, "absolute": None, "signed": None}
    result["absolute"] = dict(zip(("median", "p90", "p99", "max"),
                                  map(float, np.quantile(abs(finite), (.5, .9, .99, 1)))))
    result["signed"] = dict(zip(("min", "p01", "median", "p99", "max"),
                                map(float, np.quantile(finite, (0, .01, .5, .99, 1)))))
    return result


def broad_own_use(logq, logt):
    return ((logt >= 3.1 - 1e-12)
            & ((logt >= 6.5 - 1e-12) | (logq <= -.5 + 1e-12)
               | (logq + 2 * logt - 12 <= -6.5 + 1e-12)))


def load_queries(path, logq, logt):
    with np.load(path, allow_pickle=False) as saved:
        required = ("logQ", "logT", "eligible", "material", "density_cell_index", "temperature_cell_index")
        missing = [field for field in required if field not in saved]
        if missing:
            raise ValueError(f"Missing saved query fields: {missing}")
        data = {field: saved[field].copy() for field in saved.files}
    shape = data["logQ"].shape
    if len(shape) != 2 or not all(shape):
        raise ValueError("Saved queries must have shape (samples_per_cell, cells)")
    for field in ("logT", "eligible", "material"):
        if data[field].shape != shape:
            raise ValueError(f"Saved {field} does not match query shape")
    for field in ("eligible", "material"):
        if data[field].dtype != np.dtype(bool):
            raise ValueError(f"Saved {field} must be a boolean mask")
    for field in ("density_cell_index", "temperature_cell_index"):
        if data[field].shape != (shape[1],) or not np.issubdtype(data[field].dtype, np.integer):
            raise ValueError(f"Saved {field} must be one integer per source cell")
    density_index, temperature_index, inside = query_cells(logq, logt, data["logQ"], data["logT"])
    if not np.all(inside):
        raise ValueError("Saved dense queries must be finite and inside the source grid")
    if (not np.all(density_index == data["density_cell_index"][None, :])
            or not np.all(temperature_index == data["temperature_cell_index"][None, :])):
        raise ValueError("Saved cell indices disagree with saved query coordinates")
    cell_pairs = np.stack((data["density_cell_index"], data["temperature_cell_index"]), axis=1)
    if np.unique(cell_pairs, axis=0).shape[0] != shape[1]:
        raise ValueError("Saved source cell columns must be unique")
    density_phase = (data["logQ"] - logq[density_index]) / (logq[1] - logq[0])
    temperature_phase = (data["logT"] - logt[temperature_index]) / (logt[1] - logt[0])
    if (np.any(density_phase <= 1e-9) or np.any(density_phase >= 1 - 1e-9)
            or np.any(temperature_phase <= 1e-9) or np.any(temperature_phase >= 1 - 1e-9)):
        raise ValueError("Dense archive must be strictly off-node in both coordinates")
    if not np.array_equal(data["eligible"], broad_own_use(data["logQ"], data["logT"])):
        raise ValueError("Saved eligible flag disagrees with the inherited broad own-use rule")
    for field, values in data.items():
        if field not in ("density_cell_index", "temperature_cell_index") and values.shape != shape:
            raise ValueError(f"Saved diagnostic {field} does not match query shape")
    data["density_phase"] = density_phase
    data["temperature_phase"] = temperature_phase
    return data


def subset_masks(valid, saved):
    eligible, material = saved["eligible"].ravel(), saved["material"].ravel()
    return {"all_valid": valid, "broad_own_use": valid & eligible,
            "material_own_use": valid & eligible & material, "material_valid": valid & material}


def identity_diagnostics(state):
    pressure, density, temperature = state["P"], state["rho"], state["T"]
    first = density * state["E_x"]
    second = state["P_t"]
    entropy_temperature = temperature * state["S_t"]
    energy_temperature = state["E_t"]
    entropy_density = density * temperature * state["S_x"]
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        functions = {"dpe": (first + second - pressure) / pressure,
                     "dse": entropy_temperature / energy_temperature - 1,
                     "dsp": -entropy_density / second - 1}
        columns = {"dpe": (density**2 * state["E_rho"] + temperature * state["P_T"] - pressure) / pressure,
                   "dse": temperature * state["S_T"] / state["E_T"] - 1,
                   "dsp": -density**2 * state["S_rho"] / state["P_T"] - 1}
        balanced = {"dpe": (first + second - pressure) / np.maximum.reduce((abs(first), abs(second), abs(pressure))),
                    "dse": (entropy_temperature - energy_temperature) / np.maximum(abs(entropy_temperature), abs(energy_temperature)),
                    "dsp": (entropy_density + second) / np.maximum(abs(entropy_density), abs(second))}
        condition = {"dpe": (abs(first) + abs(second) + abs(pressure)) / abs(pressure),
                     "dse": (abs(entropy_temperature) + abs(energy_temperature)) / abs(energy_temperature),
                     "dsp": (abs(entropy_density) + abs(second)) / abs(second)}
    return {"returned_function": functions, "physical_column": columns,
            "balanced_returned_function": balanced, "condition_number": condition}


def evaluate_chunks(nodes, complete, logq, logt, saved, valid, add_radiation, chunk_size):
    if chunk_size <= 0:
        raise ValueError("Chunk size must be positive")
    query_count = saved["logQ"].size
    selected = np.flatnonzero(valid)
    state_values = {field: np.full(query_count, np.nan) for field in COMPARISON_FIELDS}
    diagnostics = {route: {field: np.full(query_count, np.nan) for field in IDENTITY_FIELDS}
                   for route in ("returned_function", "physical_column", "balanced_returned_function", "condition_number")}
    nonfinite_state = {}
    sample_q, sample_t = saved["logQ"].ravel(), saved["logT"].ravel()
    for start in range(0, selected.size, chunk_size):
        indices = selected[start:start + chunk_size]
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            native = evaluate(nodes, complete, logq, logt, sample_q[indices], sample_t[indices], add_radiation=add_radiation)
            state = thermodynamic_partials(native, sample_q[indices], sample_t[indices])
            current = identity_diagnostics(state)
        for field in state_values:
            state_values[field][indices] = state[field]
        for route, fields in current.items():
            for field, values in fields.items():
                diagnostics[route][field][indices] = values
        for field, values in state.items():
            nonfinite_state[field] = nonfinite_state.get(field, 0) + int(np.count_nonzero(~np.isfinite(values)))
    return state_values, diagnostics, nonfinite_state


def query_location(index, saved, value=None):
    sample_index, cell_column = np.unravel_index(int(index), saved["logQ"].shape)
    logq, logt = float(saved["logQ"].ravel()[index]), float(saved["logT"].ravel()[index])
    result = {"flat_index": int(index), "sample_index": int(sample_index), "archive_cell_column": int(cell_column),
              "density_cell_index": int(saved["density_cell_index"][cell_column]),
              "temperature_cell_index": int(saved["temperature_cell_index"][cell_column]),
              "logQ": logq, "logT": logt, "logRho": logq + 2 * logt - 12,
              "saved_eligible": bool(saved["eligible"].ravel()[index]),
              "saved_material": bool(saved["material"].ravel()[index])}
    if "original_flat_index" in saved:
        result["flat_index"] = int(saved["original_flat_index"].ravel()[index])
        result["archive_cell_column"] = int(saved["original_cell_columns"][cell_column])
    if value is not None:
        result["value"] = finite_number(value)
        result["nonfinite_kind"] = None if np.isfinite(value) else ("NaN" if np.isnan(value) else str(value))
    return result


def positivity_summary(values, masks, saved):
    result = {}
    for subset, mask in masks.items():
        fields = {}
        combined = np.zeros(mask.size, dtype=bool)
        for field in POSITIVE_FIELDS:
            current = values[field]
            finite = np.isfinite(current)
            failed = mask & (~finite | (current <= 0))
            combined |= failed
            finite_indices = np.flatnonzero(mask & finite)
            failures = np.flatnonzero(failed)
            finite_failures = np.flatnonzero(failed & finite)
            worst = finite_failures[np.argsort(current[finite_failures])[:10]]
            minimum = finite_indices[np.argmin(current[finite_indices])] if finite_indices.size else None
            fields[field] = {"count": int(mask.sum()), "finite": int((mask & finite).sum()),
                             "positive_finite": int((mask & finite & (current > 0)).sum()),
                             "nonpositive_finite": int((mask & finite & (current <= 0)).sum()),
                             "negative_finite": int((mask & finite & (current < 0)).sum()),
                             "zero_finite": int((mask & finite & (current == 0)).sum()),
                             "nonfinite": int((mask & ~finite).sum()), "failed": int(failed.sum()),
                             "minimum_finite": query_location(minimum, saved, current[minimum]) if minimum is not None else None,
                             "failure_examples": [query_location(index, saved, current[index]) for index in failures[:10]],
                             "most_negative_failures": [query_location(index, saved, current[index]) for index in worst]}
        result[subset] = {"count": int(mask.sum()), "any_field_failed": int(combined.sum()),
                          "all_five_finite_positive": int((mask & ~combined).sum()), "fields": fields}
    return result


def summarize_identities(diagnostics, masks, condition_limit=1000.0):
    result = {}
    for subset, mask in masks.items():
        result[subset] = {route: {field: distribution(values[mask]) for field, values in fields.items()}
                          for route, fields in diagnostics.items()}
        result[subset]["well_conditioned_returned_function"] = {}
        for field in IDENTITY_FIELDS:
            values = diagnostics["returned_function"][field]
            condition = diagnostics["condition_number"][field]
            selected = mask & np.isfinite(condition) & (condition <= condition_limit)
            result[subset]["well_conditioned_returned_function"][field] = {
                **distribution(values[selected]), "condition_limit": condition_limit,
                "conditioning_excluded": int((mask & ~selected).sum()),
                "absolute_residual_above_1e-10": int(np.count_nonzero(abs(values[selected]) > 1e-10))}
    return result


def relative_departure(candidate, reference):
    result = np.full(candidate.shape, np.nan)
    accepted = np.isfinite(candidate) & np.isfinite(reference) & (reference != 0)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        result[accepted] = candidate[accepted] / reference[accepted] - 1
    return result


def legacy_context(saved, masks, state_values):
    diagnostics = {}
    for subset, mask in masks.items():
        diagnostics[subset] = {}
        for route in ("mathematical", "physical"):
            diagnostics[subset][route] = {field: distribution(saved[f"{route}_{field}"].ravel()[mask])
                                          for field in IDENTITY_FIELDS if f"{route}_{field}" in saved}
        diagnostics[subset]["Cv_mismatch_LEGACY_function_vs_LEGACY_column_ONLY"] = (
            distribution(saved["Cv_mismatch"].ravel()[mask]) if "Cv_mismatch" in saved else None)
    departures = {}
    unavailable = []
    for field in COMPARISON_FIELDS:
        key = next((key for key in (f"legacy_{field}", field) if key in saved), None)
        if key is None:
            unavailable.append(field)
            continue
        reference = saved[key].ravel()
        difference = relative_departure(state_values[field], reference)
        departures[field] = {"saved_reference_field": key,
                             "subsets": {subset: distribution(difference[mask]) for subset, mask in masks.items()}}
    return {"saved_diagnostics_on_identical_accepted_queries": diagnostics,
            "candidate_relative_departure_from_saved_legacy_states_NOT_source_accuracy": departures,
            "unavailable_saved_legacy_state_fields": unavailable,
            "Cv_mismatch_semantics": "Saved legacy returned dE/dT divided by legacy interpolated Cv column minus one; not candidate error or source accuracy. No cache is reevaluated."}


def cell_median(values, mask, shape, absolute=False):
    selected = np.where(mask.reshape(shape), values.reshape(shape), np.nan)
    if absolute:
        selected = abs(selected)
    selected[~np.isfinite(selected)] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmedian(selected, axis=0)


def cell_failure_maps(values, mask, shape):
    selected = mask.reshape(shape)
    current = values.reshape(shape)
    reference = cell_median(values, mask, shape, absolute=True)
    finite = np.isfinite(current)
    negative = selected & finite & (current < 0)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        severity = np.where(negative, -current / np.maximum(reference[None, :], np.finfo(float).tiny), 0)
    maximum = np.max(severity, axis=0)
    maximum[~np.any(negative, axis=0)] = np.nan
    return {"negative_max_relative_to_cell_median_absolute": maximum,
            "nonfinite_count": np.sum(selected & ~finite, axis=0),
            "zero_count": np.sum(selected & finite & (current == 0), axis=0),
            "failed_count": np.sum(selected & (~finite | (current <= 0)), axis=0),
            "selected_count": np.sum(selected, axis=0), "median_absolute_reference": reference}


def make_cell_maps(variants, masks, saved, logq, logt):
    shape = saved["logQ"].shape
    maps = {"logQ_grid": logq, "logT_grid": logt,
            "density_cell_index": saved["density_cell_index"], "temperature_cell_index": saved["temperature_cell_index"]}
    for subset in ("all_valid", "material_own_use"):
        mask = masks[subset]
        for variant, (values, diagnostics) in variants.items():
            for field in IDENTITY_FIELDS:
                maps[f"{variant}__{subset}__{field}__median_absolute"] = cell_median(
                    diagnostics["returned_function"][field], mask, shape, absolute=True)
            for field in POSITIVE_FIELDS:
                for name, result in cell_failure_maps(values[field], mask, shape).items():
                    maps[f"{variant}__{subset}__{field}__{name}"] = result
        for route in ("mathematical", "physical"):
            for field in IDENTITY_FIELDS:
                if f"{route}_{field}" in saved:
                    maps[f"legacy_{route}__{subset}__{field}__median_absolute"] = cell_median(
                        saved[f"{route}_{field}"].ravel(), mask, shape, absolute=True)
    return maps


def build_models(rows, logq, logt):
    models = {}
    for variant, split_radiation in VARIANTS.items():
        low, source_valid, state = source_constraints(rows, logq, logt, split_radiation=split_radiation)
        nodes, complete, _ = construct_twists(low, source_valid, (logq[1] - logq[0]) * LN10,
                                               (logt[1] - logt[0]) * LN10, radius=2, route="mean")
        models[variant] = {"low": low, "source_valid": source_valid, "nodes": nodes,
                           "complete": complete, "source_state": state, "split_radiation": split_radiation}
    return models


def domain_summary(rows, logq, logt, saved, models):
    source_density = models["radius2_mean_total"]["source_state"]["rho"]
    printed_bad = ~np.isfinite(rows[:, :, 18]) | (abs(rows[:, :, 18] - np.log10(source_density)) >= 1e-10)
    density_index, temperature_index = saved["density_cell_index"], saved["temperature_cell_index"]
    shape = saved["logQ"].shape

    def query_mask(node_mask):
        return np.broadcast_to(valid_cells(node_mask)[density_index, temperature_index], shape).ravel()

    above_floor = saved["logT"].ravel() >= 3.2
    native_corners = query_mask(rows[:, :, 17] == 0)
    printed_corners = query_mask(~printed_bad)
    source_corners = np.logical_and.reduce([query_mask(model["source_valid"]) for model in models.values()])
    available = {variant: query_mask(model["complete"]) & above_floor for variant, model in models.items()}
    common = np.logical_and.reduce(list(available.values()))
    remaining = np.ones(common.size, dtype=bool)
    exclusions = {}
    for reason, accepted in (("below_logT_3.2", above_floor), ("non_native_corner", native_corners),
                             ("printed_density_corner_mismatch", printed_corners),
                             ("invalid_source_corner_constraint", source_corners),
                             ("outside_common_complete_twist_stencil", common)):
        exclusions[reason] = int(np.count_nonzero(remaining & ~accepted))
        remaining &= accepted
    exclusions["accepted"] = int(common.sum())
    if sum(exclusions.values()) != common.size:
        raise AssertionError("Exclusive domain exclusions do not account for every saved query")
    result = {"saved_shape": list(shape), "saved_queries": int(common.size), "saved_cells": int(shape[1]),
              "strictly_off_node_queries": int(common.size),
              "density_phase_minmax": [float(saved["density_phase"].min()), float(saved["density_phase"].max())],
              "temperature_phase_minmax": [float(saved["temperature_phase"].min()), float(saved["temperature_phase"].max())],
              "saved_broad_own_use_queries": int(saved["eligible"].sum()),
              "saved_material_queries": int(saved["material"].sum()),
              "saved_material_own_use_queries": int((saved["material"] & saved["eligible"]).sum()),
              "exclusive_query_exclusions_in_order": exclusions,
              "overlapping_query_exclusions": {"below_logT_3.2": int((~above_floor).sum()),
                                                "non_native_corner": int((~native_corners).sum()),
                                                "printed_density_corner_mismatch": int((~printed_corners).sum()),
                                                "outside_common_complete_stencil": int((~common).sum())},
              "accepted_queries_by_variant_before_common_intersection": {name: int(mask.sum()) for name, mask in available.items()},
              "common_complete_stencil_queries": int(common.sum()),
              "source_grid_shape": list(rows.shape[:2]), "source_native_nodes": int((rows[:, :, 17] == 0).sum()),
              "source_printed_density_mismatch_nodes": int(printed_bad.sum()),
              "source_printed_density_max_absolute_log10_difference": finite_number(np.nanmax(abs(rows[:, :, 18] - np.log10(source_density)))),
              "source_valid_nodes": {name: int(model["source_valid"].sum()) for name, model in models.items()},
              "complete_nodes": {name: int(model["complete"].sum()) for name, model in models.items()},
              "complete_cells": {name: int(valid_cells(model["complete"]).sum()) for name, model in models.items()},
              "low_temperature_floor_logT": 3.2,
              "printed_density_tolerance_absolute_log10": 1e-10,
              "selection": "No new queries: intersection of total and split centered-radius-2 complete cells. Core source masks remove non-native, low-T, printed-coordinate mismatches and invalid low-order constraints before both twist stages; no extrapolation.",
              "material_definition": "Inherited saved mask: legacy Prad/Pgas < 1e-3; never reclassified using candidate pressure.",
              "broad_own_use_definition": "Inherited eligible mask: logT>=3.1 and (logT>=6.5 or logQ<=-0.5 or logRho<=-6.5), with inherited 1e-12 tolerance; candidate evaluation additionally requires logT>=3.2 and a complete stencil."}
    return common, result


def select_suspicious_cells(variants, masks, saved, maximum):
    if maximum < 0:
        raise ValueError("Maximum suspicious cells must be nonnegative")
    selected = {}
    for subset in ("material_own_use", "broad_own_use", "all_valid"):
        for variant, (values, _) in variants.items():
            for field in POSITIVE_FIELDS:
                current = values[field]
                failed = np.flatnonzero(masks[subset] & (~np.isfinite(current) | (current <= 0)))
                finite_failed = failed[np.isfinite(current[failed])]
                ordered = np.concatenate((failed[~np.isfinite(current[failed])][:8],
                                          finite_failed[np.argsort(current[finite_failed])[:8]]))
                for index in ordered:
                    column = int(index % saved["logQ"].shape[1])
                    selected.setdefault(column, []).append(f"{subset}:{variant}:{field} positivity failure")
    total = variants["radius2_mean_total"][0]
    split = variants["radius2_mean_analytic_split"][0]
    for subset in ("material_own_use", "broad_own_use", "all_valid"):
        change = abs(relative_departure(split["Cv"], total["Cv"]))
        indices = np.flatnonzero(masks[subset] & np.isfinite(change))
        ordered = indices[np.argsort(change[indices])[-maximum:][::-1]] if maximum else []
        for index in ordered:
            column = int(index % saved["logQ"].shape[1])
            selected.setdefault(column, []).append(f"{subset}:large total/split Cv representation departure")
    return [{"archive_cell_column": column, "reasons": reasons}
            for column, reasons in list(selected.items())[:maximum]]


def twist_sensitivity(models, variants, valid, saved, logq, logt, selections, chunk_size):
    if not selections:
        return {"selected_cells": [], "schemes": {}, "scope": "No suspicious-cell sensitivity requested."}
    columns = np.array([selection["archive_cell_column"] for selection in selections], dtype=int)
    compact = {field: values[:, columns] if values.ndim == 2 else values[columns]
               for field, values in saved.items()}
    compact["original_flat_index"] = np.arange(valid.size).reshape(saved["logQ"].shape)[:, columns]
    compact["original_cell_columns"] = columns
    compact_valid = valid.reshape(saved["logQ"].shape)[:, columns].ravel()
    schemes = {}
    for variant, model in models.items():
        central = {field: values.reshape(saved["logQ"].shape)[:, columns].ravel()
                   for field, values in variants[variant][0].items()}
        for radius, route in ((1, "mean"), (2, "first"), (2, "second"), (3, "mean")):
            nodes, complete, _ = construct_twists(model["low"], model["source_valid"],
                                                   (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10,
                                                   radius=radius, route=route)
            scheme_valid = compact_valid & np.broadcast_to(valid_cells(complete)[
                compact["density_cell_index"], compact["temperature_cell_index"]], compact["logQ"].shape).ravel()
            values, _, nonfinite = evaluate_chunks(nodes, complete, logq, logt, compact, scheme_valid,
                                                  model["split_radiation"], chunk_size)
            masks = subset_masks(scheme_valid, compact)
            name = f"{variant}__radius{radius}_{route}"
            schemes[name] = {"accepted_saved_queries": int(scheme_valid.sum()),
                             "outside_this_scheme_complete_stencil": int((compact_valid & ~scheme_valid).sum()),
                             "positivity": positivity_summary(values, masks, compact),
                             "nonfinite_returned_fields": nonfinite,
                             "relative_departure_from_same_radiation_radius2_mean_NOT_source_accuracy": {
                                 field: {subset: distribution(relative_departure(values[field], central[field])[mask])
                                         for subset, mask in masks.items()} for field in COMPARISON_FIELDS}}
    return {"selected_cells": selections, "selected_saved_queries": int(compact_valid.sum()), "schemes": schemes,
            "scope": "Only inherited dense queries in selected suspicious cells, no new phases or duplicate original 12000-point matched audit. Comparisons are representation sensitivity, not fresh-source error."}


def save_failures(output, variants, masks, saved):
    failures = {}
    for variant, (values, _) in variants.items():
        for field in POSITIVE_FIELDS:
            indices = np.flatnonzero(masks["all_valid"] & (~np.isfinite(values[field]) | (values[field] <= 0)))
            prefix = f"{variant}__{field}"
            columns = indices % saved["logQ"].shape[1]
            failures[prefix + "__flat_index"] = indices
            failures[prefix + "__value"] = values[field][indices]
            failures[prefix + "__logQ"] = saved["logQ"].ravel()[indices]
            failures[prefix + "__logT"] = saved["logT"].ravel()[indices]
            failures[prefix + "__density_cell_index"] = saved["density_cell_index"][columns]
            failures[prefix + "__temperature_cell_index"] = saved["temperature_cell_index"][columns]
            failures[prefix + "__saved_material"] = saved["material"].ravel()[indices]
            failures[prefix + "__saved_eligible"] = saved["eligible"].ravel()[indices]
    np.savez_compressed(output / "positivity_failures.npz", **failures)


def save_figures(output, maps, variants):
    os.environ["MPLCONFIGDIR"] = str(output / ".mplconfig")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    logq, logt = maps["logQ_grid"], maps["logT_grid"]
    density_edges = logq[:, None] + 2 * logt[None, :] - 12
    temperature_edges = np.broadcast_to(logt[None, :], density_edges.shape)
    cell_shape = (len(logq) - 1, len(logt) - 1)
    density_centers = (logq[:-1] + logq[1:]) / 2
    temperature_centers = (logt[:-1] + logt[1:]) / 2
    physical_centers = density_centers[:, None] + 2 * temperature_centers[None, :] - 12
    own_use = broad_own_use(density_centers[:, None], temperature_centers[None, :])
    norm = TwoSlopeNorm(vmin=-7, vcenter=-3, vmax=1)
    cmap = plt.colormaps["coolwarm"].copy()
    cmap.set_bad("0.82")

    def expand(values):
        cells = np.full(cell_shape, np.nan)
        cells[maps["density_cell_index"], maps["temperature_cell_index"]] = values
        return cells

    def panel(axis, values, title):
        cells = expand(values)
        image = axis.pcolormesh(density_edges.T, temperature_edges.T,
                               np.ma.masked_invalid(np.log10(np.maximum(cells.T, 1e-12))),
                               shading="flat", cmap=cmap, norm=norm, rasterized=True)
        axis.contour(physical_centers.T, np.broadcast_to(temperature_centers[None, :], cell_shape).T,
                     own_use.T, levels=[.5], colors="white", linewidths=.7)
        axis.set(title=title, xlabel="log10 rho [g/cm3]", ylabel="log10 T [K]", xlim=(-9, 5), ylim=(3.05, 7.8))
        axis.set_facecolor("0.82")
        return image

    def finish(figure, axes, image, title, label, name):
        figure.colorbar(image, ax=axes, label=label, ticks=np.arange(-7, 2), extend="both", shrink=.8)
        figure.suptitle(title, fontsize=11)
        figure.savefig(output / f"{name}.png", dpi=180)
        figure.savefig(output / f"{name}.pdf")
        plt.close(figure)

    for variant in variants:
        figure, axes = plt.subplots(2, 3, figsize=(13, 7), layout="constrained")
        for row, subset in enumerate(("all_valid", "material_own_use")):
            for column, field in enumerate(IDENTITY_FIELDS):
                image = panel(axes[row, column], maps[f"{variant}__{subset}__{field}__median_absolute"], f"{subset}: {field}")
        finish(figure, axes, image, f"{variant}: inherited dense off-node returned-function identities\nGrey = no accepted subset query; white = inherited broad own-use boundary",
               "log10 per-cell median absolute residual", f"{variant}_identity_medians")
        figure, axes = plt.subplots(2, 5, figsize=(18, 7), layout="constrained")
        for row, subset in enumerate(("all_valid", "material_own_use")):
            for column, field in enumerate(POSITIVE_FIELDS):
                prefix = f"{variant}__{subset}__{field}"
                image = panel(axes[row, column], maps[prefix + "__negative_max_relative_to_cell_median_absolute"], f"{subset}: {field}")
                exceptional = (maps[prefix + "__nonfinite_count"] + maps[prefix + "__zero_count"]) > 0
                if np.any(exceptional):
                    source_density = maps["density_cell_index"][exceptional]
                    source_temperature = maps["temperature_cell_index"][exceptional]
                    axes[row, column].scatter(physical_centers[source_density, source_temperature],
                                              temperature_centers[source_temperature], s=3, c="black", marker="s")
        finish(figure, axes, image, f"{variant}: negative-value failure maxima, same saved queries\nGrey = no finite negative failure; black squares = nonfinite or exact zero; this is not source accuracy",
               "log10 max(-value) / median(abs(value)) per cell and subset", f"{variant}_positivity_failure_maxima")
    if all(f"legacy_{route}__all_valid__{field}__median_absolute" in maps
           for route in ("mathematical", "physical") for field in IDENTITY_FIELDS):
        figure, axes = plt.subplots(4, 3, figsize=(13, 12), layout="constrained")
        for row, (subset, route) in enumerate((subset, route) for subset in ("all_valid", "material_own_use")
                                             for route in ("mathematical", "physical")):
            for column, field in enumerate(IDENTITY_FIELDS):
                image = panel(axes[row, column], maps[f"legacy_{route}__{subset}__{field}__median_absolute"], f"{subset}, legacy {route}: {field}")
        finish(figure, axes, image, "Saved legacy context on the identical candidate-accepted dense queries\nMathematical = legacy returned-function identities; physical = legacy column closure; neither is candidate error",
               "log10 per-cell median absolute legacy residual", "saved_legacy_context_medians")


def run_audit(table, archive, output, chunk_size=16384, suspicious_cells=32, figures=True, query_archive=True):
    started = time.perf_counter()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    header, logq, logt, rows = read_table(table)
    if not np.allclose(header[1:3], [.7, .02], rtol=0, atol=1e-12):
        raise ValueError("Dense audit requires the inherited X=0.70, Z=0.02 source table")
    saved = load_queries(archive, logq, logt)
    models = build_models(rows, logq, logt)
    valid, domain = domain_summary(rows, logq, logt, saved, models)
    if not np.any(valid):
        raise ValueError("No saved dense query has a common complete radius-2 stencil")
    masks = subset_masks(valid, saved)
    results = {"status": "table-only dense potential interpolation audit; no source-accuracy or runtime claims",
               "provenance": {"source": str(Path(table).resolve()), "source_sha256": checksum(table),
                              "saved_archive": str(Path(archive).resolve()), "saved_archive_sha256": checksum(archive),
                              "core_sha256": checksum(Path(__file__).with_name("freeeos_potential.py")),
                              "audit_sha256": checksum(__file__), "python": sys.version, "numpy": np.__version__,
                              "command": [sys.executable, *sys.argv], "chunk_size": chunk_size,
                              "composition": {"X": .7, "Z": .02}, "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())},
               "domain": domain, "subset_counts": {name: int(mask.sum()) for name, mask in masks.items()},
               "variants": {}, "interpretation": {
                   "positivity": "Failures are nonfinite or <=0 in the candidate functions on accepted saved queries; E and S sign tests are requested diagnostics, not a general thermodynamic stability theorem.",
                   "identities": "Both candidate routes use one potential. Balanced and condition-filtered returned-function residuals are reported separately; algebraic closure is not physical/source accuracy.",
                   "departures": "Relative departures exist only if the saved archive contains legacy state values. Cv_mismatch is legacy function-vs-column closure, never candidate accuracy.",
                   "map_normalization": "coolwarm, log10 display centered at -3, vmin=-7, vmax=+1; identity maps show per-cell median absolute residual. Failure maps show max(-value)/median(abs(value)) within the same cell/subset, only finite negative failures; zeros/nonfinite failures are black."}}
    variants = {}
    for variant, model in models.items():
        values, diagnostics, nonfinite = evaluate_chunks(model["nodes"], model["complete"], logq, logt, saved,
                                                       valid, model["split_radiation"], chunk_size)
        variants[variant] = (values, diagnostics)
        context = legacy_context(saved, masks, values)
        results["variants"][variant] = {
            "split_radiation_before_twists": model["split_radiation"],
            "add_analytic_radiation_after_evaluation": model["split_radiation"],
            "positivity": positivity_summary(values, masks, saved),
            "identities": summarize_identities(diagnostics, masks), "nonfinite_returned_fields_all_valid": nonfinite,
            "candidate_relative_departure_from_saved_legacy_states_NOT_source_accuracy": context[
                "candidate_relative_departure_from_saved_legacy_states_NOT_source_accuracy"]}
        if variant == "radius2_mean_total":
            results["legacy_context"] = context
        print(f"Evaluated {variant}: {int(valid.sum())} inherited queries", file=sys.stderr, flush=True)
    total = variants["radius2_mean_total"][0]
    split = variants["radius2_mean_analytic_split"][0]
    results["analytic_split_relative_departure_from_total_NOT_source_accuracy"] = {
        field: {subset: distribution(relative_departure(split[field], total[field])[mask])
                for subset, mask in masks.items()} for field in COMPARISON_FIELDS}
    selections = select_suspicious_cells(variants, masks, saved, suspicious_cells)
    sensitivity = twist_sensitivity(models, variants, valid, saved, logq, logt, selections, chunk_size)
    (output / "suspicious_cell_twist_sensitivity.json").write_text(json.dumps(sensitivity, indent=2, allow_nan=False) + "\n")
    results["twist_sensitivity"] = {"selected_cells": len(selections), "selected_saved_queries": sensitivity.get("selected_saved_queries", 0),
                                   "artifact": "suspicious_cell_twist_sensitivity.json", "scope": sensitivity["scope"]}
    maps = make_cell_maps(variants, masks, saved, logq, logt)
    np.savez_compressed(output / "cell_maps.npz", **maps)
    save_failures(output, variants, masks, saved)
    if query_archive:
        query_products = {field: saved[field] for field in ("logQ", "logT", "eligible", "material", "density_cell_index", "temperature_cell_index")}
        query_products["common_valid"] = valid.reshape(saved["logQ"].shape)
        for variant, (values, diagnostics) in variants.items():
            for field, current in values.items():
                query_products[f"{variant}__{field}"] = current.reshape(saved["logQ"].shape)
            for route, fields in diagnostics.items():
                for field, current in fields.items():
                    query_products[f"{variant}__{route}__{field}"] = current.reshape(saved["logQ"].shape)
        np.savez_compressed(output / "query_diagnostics.npz", **query_products)
    if figures:
        save_figures(output, maps, variants)
    results["elapsed_seconds"] = time.perf_counter() - started
    (output / "results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, default=ROOT / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data")
    parser.add_argument("--archive", type=Path, default=ROOT / "freeeos_report/plots/offgrid_sample_data.npz")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--chunk-size", type=int, default=16384)
    parser.add_argument("--suspicious-cells", type=int, default=32)
    parser.add_argument("--no-figures", action="store_true")
    parser.add_argument("--no-query-archive", action="store_true")
    arguments = parser.parse_args()
    if arguments.chunk_size <= 0 or arguments.suspicious_cells < 0:
        parser.error("chunk size must be positive and suspicious cells nonnegative")
    results = run_audit(arguments.table, arguments.archive, arguments.output, arguments.chunk_size,
                        arguments.suspicious_cells, not arguments.no_figures, not arguments.no_query_archive)
    summary = {"status": results["status"], "saved_queries": results["domain"]["saved_queries"],
               "subset_counts": results["subset_counts"], "elapsed_seconds": results["elapsed_seconds"],
               "positivity_failures": {variant: {subset: {field: summary["failed"] for field, summary in data["fields"].items()}
                                                for subset, data in details["positivity"].items()}
                                       for variant, details in results["variants"].items()},
               "output": str(arguments.output.resolve())}
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
