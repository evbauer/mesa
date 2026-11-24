"""Archived composition holdout feasibility; no fresh-source accuracy or runtime EOS."""

import argparse
import json
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np

from composition import (
    SUPPORTS, TARGET_X, combine, composition_partials, composition_weights, gauge_jet,
    grid_spacing,
)
from endpoint_audit import ENDPOINT_FIELDS, audit_printed_grid, comparison_header
from freeeos_potential import (
    CRAD, EVAL_ORDERS, checksum, construct_twists, evaluate, identity_residuals,
    query_cells, read_table, source_constraints as leadsource_constraints,
    thermodynamics, valid_cells,
)
from run_audit import own_use
from validate_table import material_nodes, table_state


FIELDS = ("P", "Pgas", "E", "S", "Cv", "Cp", "chiRho", "chiT", "gamma1", "gamma3", "grad_ad")
POSITIVE = ("P", "Pgas", "Cv", "Cp", "chiRho", "P_rho", "E", "S", "gamma1", "gamma3", "grad_ad")
WORKSPACE = Path(__file__).resolve().parents[5]
OUTPUT = WORKSPACE / "freeeos_biquintic_progress/experiments/012_composition_control"


class SourceGridError(ValueError):
    def __init__(self, metadata, failures):
        super().__init__("Source header/grid gate failed; no repair, trimming or resampling authorized")
        self.metadata = metadata
        self.failures = failures


def header_endpoint_failures(header, logq, logt):
    failures = []
    for name, declared, actual in (("logQ_max", header[9], logq[-1]), ("logT_max", header[5], logt[-1])):
        if not np.isclose(declared, actual, rtol=0, atol=1e-10):
            failures.append({"field": name, "declared": float(declared), "count_start_step_endpoint": float(actual)})
    return failures


def stats(values):
    values = np.asarray(values)
    finite = np.isfinite(values)
    result = {"count": int(values.size), "nonfinite": int((~finite).sum())}
    if finite.any():
        result.update(zip(("abs_median", "abs_p90", "abs_p99", "abs_max"),
                          map(float, np.quantile(abs(values[finite]), (.5, .9, .99, 1)))))
        result["signed_median"] = float(np.median(values[finite]))
    return result


def comparison(candidate, reference, sample_q, sample_t, selected=None):
    if selected is None:
        selected = np.ones(len(sample_q), bool)
    result = {}
    for field in FIELDS:
        first, second = candidate[field][selected], reference[field][selected]
        with np.errstate(divide="ignore", invalid="ignore"):
            residual = (first - second) / np.maximum(abs(second), 1e-300)
        result[field] = stats(residual)
        finite_indices = np.flatnonzero(np.isfinite(residual))
        if len(finite_indices):
            worst = finite_indices[np.argmax(abs(residual[finite_indices]))]
            result[field]["worst"] = {"logQ": float(sample_q[selected][worst]),
                                      "logT": float(sample_t[selected][worst]),
                                      "candidate": float(first[worst]), "reference": float(second[worst]),
                                      "signed_relative_difference": float(residual[worst])}
    return result


def diagnostics(state, sample_q, sample_t, selected=None):
    if selected is None:
        selected = np.ones(len(sample_q), bool)
    reduced = {field: values[selected] for field, values in state.items()}
    result = {"count": int(selected.sum()),
              "identities_balanced": {name: stats(values) for name, values in
                                      identity_residuals(reduced, balanced=True).items()},
              "identities_pressure_normalized": {name: stats(values) for name, values in
                                                 identity_residuals(reduced).items()},
              "positivity": {}}
    for field in POSITIVE:
        failed = (~np.isfinite(reduced[field])) | (reduced[field] <= 0)
        result["positivity"][field] = {"count": int(failed.sum()), "examples": [
            {"logQ": float(sample_q[selected][index]), "logT": float(sample_t[selected][index]),
             "value": float(reduced[field][index]) if np.isfinite(reduced[field][index]) else None}
            for index in np.flatnonzero(failed)[:5]]}
    return result


def printed_state(rows, sample_t):
    state = table_state(rows)
    state["P"] = state["Pgas"] + CRAD * (10.0**sample_t)**4 / 3
    return state


def column_control(source_states, weights):
    return {field: combine({coordinate: state[field] for coordinate, state in source_states.items()}, weights)
            for field in FIELDS}


def bilinear_columns(rows, logq, logt, sample_q, sample_t):
    density_index, temperature_index, _ = query_cells(logq, logt, sample_q, sample_t)
    density_fraction = (sample_q - logq[density_index]) / (logq[1] - logq[0])
    temperature_fraction = (sample_t - logt[temperature_index]) / (logt[1] - logt[0])
    result = {field: np.zeros(len(sample_q)) for field in FIELDS}
    for density_end in range(2):
        for temperature_end in range(2):
            weight = (density_fraction if density_end else 1 - density_fraction)
            weight = weight * (temperature_fraction if temperature_end else 1 - temperature_fraction)
            corner = printed_state(rows[density_index + density_end, temperature_index + temperature_end],
                                   logt[temperature_index + temperature_end])
            for field in FIELDS:
                if field != "P":
                    result[field] += weight * corner[field]
    result["P"] = result["Pgas"] + CRAD * (10.0**sample_t)**4 / 3
    return result


def output_suffix(arguments):
    return "_audited_endpoint" if arguments.allow_audited_endpoint_metadata else ""


def load_sources(mesa, metallicity, split, allow_audited_endpoint_metadata=False):
    rows_by_x, low_by_x, metadata = {}, {}, {}
    reference_grid = None
    failures = []
    audited_resolution = allow_audited_endpoint_metadata and metallicity == .02
    for coordinate in (.5, .6, .7, .8, .9):
        source = mesa / f"data/eosFreeEOS_data/mesa-FreeEOS_{round(100 * metallicity):02d}z{round(100 * coordinate):02d}x.data"
        digest = checksum(source)
        try:
            header, logq, logt, rows = read_table(source, validate_endpoints=not audited_resolution)
            endpoint_audit = None
            if audited_resolution:
                if not np.allclose(header[1:3], (coordinate, .02), rtol=0, atol=1e-12):
                    raise ValueError("Audited endpoint resolution requires the authorized X/Z=.02 header")
                with source.open() as stream:
                    endpoint_audit = audit_printed_grid(stream, header, logq, logt)
                endpoint_audit["loader"] = "lead read_table(validate_endpoints=False), bounded Z=.02 authorization and printed-grid audit before constraints use rows"
                endpoint_audit["allow_audited_endpoint_metadata"] = True
                endpoint_audit["core_validate_endpoints"] = False
        except ValueError as error:
            with source.open() as stream:
                next(stream)
                header = np.fromstring(next(stream), sep=" ")
                records = sum(len(columns) == 22 and columns[0][0] in "+-0123456789"
                              for columns in (line.split() for line in stream))
            reader_error = "audited_endpoint_reader_error" if audited_resolution else "lead_read_table_error"
            metadata[str(coordinate)] = {"path": str(source), "sha256": digest, "header": header.tolist(),
                                         reader_error: str(error), "printed_record_count": records,
                                         "header_expected_records": int(header[3] * header[7])}
            failures.append({"X": coordinate, "kind": "audited_endpoint_reader_rejected_source" if audited_resolution else "lead_read_table_rejected_source",
                             "error": str(error), "printed_record_count": records,
                             "header_expected_records": int(header[3] * header[7])})
            if reference_grid is None:
                reference_grid = (comparison_header(header, audited_resolution),
                                  header[8] + np.arange(int(header[7])) * header[10],
                                  header[4] + np.arange(int(header[3])) * header[6])
            continue
        if not np.allclose(header[1:3], (coordinate, metallicity), rtol=0, atol=1e-12):
            raise ValueError(f"Unexpected composition header: {source}")
        grid_header = comparison_header(header, audited_resolution)
        endpoint_failures = header_endpoint_failures(header, logq, logt)
        if endpoint_failures and endpoint_audit is None:
            failures.append({"X": coordinate, "kind": "internally_inconsistent_header_endpoint", "details": endpoint_failures})
        if reference_grid is not None:
            if not all(np.array_equal(first, second) for first, second in
                       zip(reference_grid, (grid_header, logq, logt))):
                failures.append({"X": coordinate, "kind": "header_or_coordinate_grid_differs_from_X_050",
                                 "header_grid_matches": bool(np.array_equal(reference_grid[0], grid_header)),
                                 "logQ_grid_matches": bool(np.array_equal(reference_grid[1], logq)),
                                 "logT_grid_matches": bool(np.array_equal(reference_grid[2], logt))})
        else:
            reference_grid = (grid_header, logq, logt)
        low, valid, _ = leadsource_constraints(rows, logq, logt, split_radiation=split)
        if endpoint_audit is not None:
            for exclusion in endpoint_audit["known_coordinate_exclusions"]:
                if np.any(valid[exclusion["iq"]]):
                    raise ValueError("Known iq=520 coordinate discrepancy was not excluded by the unchanged lead constraints")
                exclusion["excluded_by_unchanged_lead_source_mask"] = True
        rows_by_x[coordinate], low_by_x[coordinate] = rows, low
        metadata[str(coordinate)] = {"path": str(source), "sha256": digest, "header": header.tolist(),
                                     "clean_source_nodes": int(valid.sum()),
                                     "header_endpoint_failures": endpoint_failures,
                                     "derived_logQ_endpoint": float(logq[-1]), "derived_logT_endpoint": float(logt[-1])}
        metadata[str(coordinate)]["endpoint_metadata_audit"] = endpoint_audit
        metadata[str(coordinate)]["source_reader"] = "lead_reader_explicit_audited_endpoint_option" if audited_resolution else "lead_default_reader"
        metadata[str(coordinate)]["core_validate_endpoints"] = not audited_resolution
        metadata[str(coordinate)]["header_comparison_ignored_endpoint_fields"] = list(ENDPOINT_FIELDS) if audited_resolution else []
        metadata[str(coordinate)]["resolved_display_endpoint_discrepancies"] = endpoint_failures if endpoint_audit else []
        if endpoint_failures and endpoint_audit:
            print(f"X={coordinate:.2f}/Z={metallicity:.2f}: EXPLICIT AUDITED display-endpoint resolution: {endpoint_failures}", flush=True)
        metadata[str(coordinate)]["valid"] = valid
    if failures:
        for entry in metadata.values():
            entry.pop("valid", None)
            if checksum(entry["path"]) != entry["sha256"]:
                raise ValueError(f"Source changed during preflight: {entry['path']}")
        raise SourceGridError(metadata, failures)
    common_source = np.logical_and.reduce([entry.pop("valid") for entry in metadata.values()])
    common_material = np.logical_and.reduce([material_nodes(rows, logq, logt) for rows in rows_by_x.values()])
    models, complete_masks = {}, []
    for coordinate, low in low_by_x.items():
        nodes, complete, _ = construct_twists(low, common_source, *grid_spacing(logq, logt), radius=2, route="mean")
        models[coordinate] = nodes
        complete_masks.append(complete)
    complete = np.logical_and.reduce(complete_masks)
    return logq, logt, rows_by_x, models, complete, common_source, common_material, metadata


def finite_difference_identities(nodes, complete, logq, logt, sample_q, sample_t, split):
    from run_audit import LN10

    sample_q, sample_t = sample_q[:1000], sample_t[:1000]
    state = thermodynamics(evaluate(nodes, complete, logq, logt, sample_q, sample_t,
                                    add_radiation=split), sample_q, sample_t)
    result = {}
    for step in (1e-4, 3e-5, 1e-5):
        numerical = state.copy()
        for suffix, delta_q, delta_t, denominator in (
            ("_rho", step / LN10, 0, state["rho"]),
            ("_T", -2 * step / LN10, step / LN10, state["T"]),
        ):
            plus = thermodynamics(evaluate(nodes, complete, logq, logt, sample_q + delta_q, sample_t + delta_t,
                                           add_radiation=split), sample_q + delta_q, sample_t + delta_t)
            minus = thermodynamics(evaluate(nodes, complete, logq, logt, sample_q - delta_q, sample_t - delta_t,
                                            add_radiation=split), sample_q - delta_q, sample_t - delta_t)
            for field in ("P", "E", "S"):
                numerical[field + suffix] = (plus[field] - minus[field]) / (2 * step * denominator)
        result[str(step)] = {name: stats(values) for name, values in
                             identity_residuals(numerical, balanced=True).items()}
    return {"count": len(sample_q), "same_potential_returned_function_NOT_accuracy": result}


def partial_x_check(scheme, models, complete, logq, logt, sample_q, sample_t, split):
    weights, partial_weights = composition_weights(scheme, warn=False)
    native = evaluate(combine(models, weights), complete, logq, logt, sample_q, sample_t, add_radiation=split)
    partial_native = evaluate(combine(models, partial_weights), complete, logq, logt, sample_q, sample_t)
    state = thermodynamics(native, sample_q, sample_t)
    analytic = composition_partials(native, partial_native, state)
    steps = {}
    for step in (1e-3, 1e-4, 1e-5):
        plus_weights, _ = composition_weights(scheme, TARGET_X + step, warn=False)
        minus_weights, _ = composition_weights(scheme, TARGET_X - step, warn=False)
        plus = thermodynamics(evaluate(combine(models, plus_weights), complete, logq, logt, sample_q, sample_t,
                                      add_radiation=split), sample_q, sample_t)
        minus = thermodynamics(evaluate(combine(models, minus_weights), complete, logq, logt, sample_q, sample_t,
                                       add_radiation=split), sample_q, sample_t)
        steps[str(step)] = {field: stats(((plus[field] - minus[field]) / (2 * step) - analytic[field]) /
                                       np.maximum.reduce([np.ones(len(sample_q)), abs(state[field]), abs(analytic[field])]))
                            for field in FIELDS}
    return {"fixed_rho_T_not_physical_composition_validation": steps}, analytic


def gauge_check(models, weights, complete, logq, logt, sample_q, sample_t, split):
    constants = {coordinate: 1e10 * coordinate**2 for coordinate in weights}
    slopes = {coordinate: 1e6 * (1 + coordinate**3) for coordinate in weights}
    constant, slope = combine(constants, weights), combine(slopes, weights)
    nodes = combine(models, weights)
    native = evaluate(nodes, complete, logq, logt, sample_q, sample_t, add_radiation=split)
    state = thermodynamics(native, sample_q, sample_t)
    shifted_jets = {coordinate: gauge_jet(sample_q, sample_t, constants[coordinate], slopes[coordinate])
                    for coordinate in weights}
    shifted = thermodynamics(native + combine(shifted_jets, weights), sample_q, sample_t)
    result = {"exact_analytic_composition_gauge": {
        field: stats((shifted[field] - state[field] - expected) /
                     np.maximum.reduce([np.ones(len(sample_q)), abs(state[field]), np.full(len(sample_q), abs(expected))]))
        for field, expected in (("P", 0), ("Cv", 0), ("E", constant), ("S", -slope))}}
    node_temperature = np.broadcast_to(logt[None, :], complete.shape)
    shifted_nodes = nodes + gauge_jet(None, node_temperature, constant, slope)[:3, :3]
    nodal_shift = thermodynamics(evaluate(shifted_nodes, complete, logq, logt, sample_q, sample_t,
                                         add_radiation=split), sample_q, sample_t)
    result["BT_nodal_tau_polynomial_gauge_drift_NOT_exact_invariance"] = {
        field: stats((nodal_shift[field] - state[field] - expected) /
                     np.maximum.reduce([np.ones(len(sample_q)), abs(state[field]), np.full(len(sample_q), abs(expected))]))
        for field, expected in (("P", 0), ("Cv", 0), ("E", constant), ("S", -slope))}
    result["mixed_gauge_A_erg_per_g"] = float(constant)
    result["mixed_gauge_B_erg_per_g_K"] = float(slope)
    return result


def evaluate_metallicity(arguments, metallicity, samples):
    split = not arguments.total_potential
    try:
        logq, logt, rows, models, complete, common_source, material, metadata = load_sources(
            arguments.mesa, metallicity, split, arguments.allow_audited_endpoint_metadata)
    except SourceGridError as error:
        source = arguments.mesa / f"data/eosFreeEOS_data/mesa-FreeEOS_{round(100 * metallicity):02d}z70x.data"
        audited_resolution = arguments.allow_audited_endpoint_metadata and metallicity == .02
        target_header, logq, logt, target_rows = read_table(source, validate_endpoints=not audited_resolution)
        if audited_resolution:
            if not np.allclose(target_header[1:3], (TARGET_X, .02), rtol=0, atol=1e-12):
                raise ValueError("Audited target selection requires the authorized .70/.02 header")
            with source.open() as stream:
                audit_printed_grid(stream, target_header, logq, logt)
        target_low, target_valid, _ = leadsource_constraints(target_rows, logq, logt, split_radiation=split)
        _, target_complete, _ = construct_twists(target_low, target_valid, *grid_spacing(logq, logt), radius=2, route="mean")
        query_density, query_temperature, inside = query_cells(logq, logt, samples["logQ"], samples["logT"])
        direct_count = int((inside & own_use(samples["logQ"], samples["logT"]) &
                            valid_cells(target_complete)[query_density, query_temperature]).sum())
        if metallicity == .02 and direct_count != 6713:
            raise ValueError("Independent archived direct-.70 selection changed")
        result = {"Z": metallicity, "input_gate": "FAILED", "matching_header_grids": False,
                  "source_sha256_and_headers": error.metadata, "header_grid_failures": error.failures,
                  "all_source_sha256_reverified_after_evaluation": True,
                  "domain": {"saved_original_queries": len(samples),
                             "direct_07_saved_own_use_radius2_queries": direct_count,
                             "sampled_SOURCE_07_nodes": 0},
                  "holdout_status": "NOT EVALUATED: strict matching header/grid gate failed; no silent header repair, trimming or resampling",
                  "schemes": {}}
        (arguments.output / f"Z{round(100 * metallicity):02d}_input_gate{output_suffix(arguments)}.json").write_text(
            json.dumps(result, indent=2, allow_nan=False) + "\n")
        print(f"Z={metallicity:.2f}: HEADER/GRID GATE FAILED; composition holdout not evaluated", flush=True)
        return result
    cells = valid_cells(complete)
    material_cells = valid_cells(material)
    node_own_use = own_use(logq[:, None], logt[None, :])
    node_selection = cells & material_cells & node_own_use[:-1, :-1]
    density_index, temperature_index = np.nonzero(node_selection)
    generator = np.random.default_rng(20261028)
    choice = generator.choice(len(density_index), min(arguments.sample_nodes, len(density_index)), replace=False)
    density_index, temperature_index = density_index[choice], temperature_index[choice]
    node_q, node_t = logq[density_index], logt[temperature_index]
    if not len(node_q):
        raise ValueError("No shared material own-use source-node holdout domain")
    node_sources = {coordinate: printed_state(source[density_index, temperature_index], node_t)
                    for coordinate, source in rows.items()}
    original_q, original_t = samples["logQ"], samples["logT"]
    query_density, query_temperature, inside = query_cells(logq, logt, original_q, original_t)
    target_low, target_valid, _ = leadsource_constraints(rows[TARGET_X], logq, logt, split_radiation=split)
    _, target_complete, _ = construct_twists(
        target_low, target_valid, *grid_spacing(logq, logt), radius=2, route="mean")
    saved_selection = inside & own_use(original_q, original_t) & valid_cells(target_complete)[query_density, query_temperature]
    if metallicity == .02 and int(saved_selection.sum()) != 6713:
        raise ValueError("Inherited .02 direct-.7 saved own-use selection no longer has 6713 queries")
    common_queries = saved_selection & cells[query_density, query_temperature]
    query_q, query_t = original_q[common_queries], original_t[common_queries]
    query_material = material_cells[query_density[common_queries], query_temperature[common_queries]]
    if not len(query_q):
        raise ValueError("No shared saved off-node domain")
    direct_node = thermodynamics(evaluate(models[TARGET_X], complete, logq, logt, node_q, node_t,
                                         add_radiation=split), node_q, node_t)
    direct_query = thermodynamics(evaluate(models[TARGET_X], complete, logq, logt, query_q, query_t,
                                          add_radiation=split), query_q, query_t)
    query_columns = {coordinate: bilinear_columns(source, logq, logt, query_q, query_t)
                     for coordinate, source in rows.items()}
    audited_resolution = arguments.allow_audited_endpoint_metadata and metallicity == .02
    result = {"Z": metallicity, "input_gate": "PASSED_AUDITED_ENDPOINT_METADATA" if audited_resolution else "PASSED",
              "source_sha256_and_headers": metadata,
              "matching_header_grids": not any(entry["header_endpoint_failures"] for entry in metadata.values()),
              "matching_authoritative_header_fields_and_coordinate_grids": True,
              "endpoint_metadata_override_explicitly_requested": audited_resolution, "domain": {
                  "source_constraint_floor_logT": 3.2,
                  "common_complete_node_min_logT": float(logt[np.any(complete, axis=0)].min()),
                  "common_complete_cell_min_logT": float(logt[:-1][np.any(cells, axis=0)].min()),
                  "native_temperature_step": float(logt[1] - logt[0]),
                  "low_temperature_coverage_gap": "[3.1,3.28) is not covered by this centered radius2/two-stage BQ prototype; an upstream FreeEOS fix alone does not close it",
                  "common_clean_source_nodes": int(common_source.sum()), "common_complete_nodes": int(complete.sum()),
                  "common_complete_cells": int(cells.sum()), "common_material_own_use_source_nodes": int(node_selection.sum()),
                  "sampled_SOURCE_07_nodes": len(node_q), "saved_original_queries": len(samples),
                  "direct_07_saved_own_use_radius2_queries": int(saved_selection.sum()),
                  "saved_queries_common_composition_domain": len(query_q),
                  "saved_queries_excluded_by_common_domain": int((saved_selection & ~common_queries).sum()),
                  "common_material_corner_queries": int(query_material.sum())},
              "direct_07_nodal_redundant_column_closure_NOT_composition_holdout": comparison(
                  direct_node, node_sources[TARGET_X], node_q, node_t),
              "direct_07_offnode_diagnostics": diagnostics(direct_query, query_q, query_t), "schemes": {}}
    artifact = {"source_node_density_index": density_index, "source_node_temperature_index": temperature_index,
                "source_node_logQ": node_q, "source_node_logT": node_t,
                "saved_direct_07_selection": saved_selection, "saved_common_selection": common_queries,
                "offnode_logQ": query_q, "offnode_logT": query_t, "offnode_material_corners": query_material,
                "common_source_mask": common_source, "common_complete_mask": complete, "common_material_mask": material}
    for field in FIELDS:
        artifact["source_07__" + field] = node_sources[TARGET_X][field]
        artifact["direct_07_offnode__" + field] = direct_query[field]
    for scheme in SUPPORTS:
        weights, partial_weights = composition_weights(scheme)
        mixed_nodes = combine(models, weights)
        node_native = evaluate(mixed_nodes, complete, logq, logt, node_q, node_t, add_radiation=split)
        query_native = evaluate(mixed_nodes, complete, logq, logt, query_q, query_t, add_radiation=split)
        node_state, query_state = thermodynamics(node_native, node_q, node_t), thermodynamics(query_native, query_q, query_t)
        node_control, query_control = column_control(node_sources, weights), column_control(query_columns, weights)
        schemes = {"weights": {str(key): float(value) for key, value in weights.items()},
                   "weight_sum": float(sum(weights.values())), "partial_weight_sum": float(sum(partial_weights.values())),
                   "negative_weights_warning": "NONCONVEX: positivity is not inherited" if min(weights.values()) < 0 else None,
                   "SOURCE_07_archived_composition_holdout": {
                       "mixed_potential_vs_printed_07": comparison(node_state, node_sources[TARGET_X], node_q, node_t),
                       "same_X_weights_physical_columns_vs_printed_07_NOT_MESA_X_scheme": comparison(
                           node_control, node_sources[TARGET_X], node_q, node_t),
                       "mixed_potential_vs_same_X_weights_columns": comparison(node_state, node_control, node_q, node_t),
                       "mixed_potential_diagnostics": diagnostics(node_state, node_q, node_t)},
                   "saved_offnode_NOT_accuracy": {}}
        for subset, selected in (("common_own_use", np.ones(len(query_q), bool)), ("common_material_corners", query_material)):
            schemes["saved_offnode_NOT_accuracy"][subset] = {
                "mixed_vs_direct_07_biquintic": comparison(query_state, direct_query, query_q, query_t, selected),
                "same_X_weights_bilinear_columns_vs_direct_07_biquintic_NOT_MESA": comparison(
                    query_control, direct_query, query_q, query_t, selected),
                "mixed_diagnostics": diagnostics(query_state, query_q, query_t, selected)}
        schemes["finite_difference_rho_T_identities"] = finite_difference_identities(
            mixed_nodes, complete, logq, logt, query_q, query_t, split)
        source_jets = {coordinate: evaluate(models[coordinate], complete, logq, logt, query_q[:1000], query_t[:1000],
                                            add_radiation=split) for coordinate in weights}
        composed_jets = combine(source_jets, weights)
        schemes["nodal_mix_vs_mix_of_evaluated_source_jets"] = {
            str(order): stats((query_native[order][:1000] - composed_jets[order]) /
                              np.maximum.reduce([np.ones(1000 if len(query_q) >= 1000 else len(query_q)),
                                                 abs(query_native[order][:1000]), abs(composed_jets[order])]))
            for order in EVAL_ORDERS}
        schemes["same_radiation_coefficient"] = float(CRAD * sum(weights.values()))
        schemes["composition_gauge_control"] = gauge_check(
            models, weights, complete, logq, logt, query_q[:1000], query_t[:1000], split)
        if arguments.derivative_x:
            check, partials = partial_x_check(scheme, models, complete, logq, logt, query_q[:1000], query_t[:1000], split)
            schemes["optional_numerical_composition_partial"] = check
            for field, values in partials.items():
                artifact[scheme + "__partial_X__" + field] = values
        for field in FIELDS:
            artifact[scheme + "__source_node__" + field] = node_state[field]
            artifact[scheme + "__source_node_column_control__" + field] = node_control[field]
            artifact[scheme + "__offnode__" + field] = query_state[field]
            artifact[scheme + "__offnode_bilinear_column_control__" + field] = query_control[field]
        result["schemes"][scheme] = schemes
    for entry in metadata.values():
        if checksum(entry["path"]) != entry["sha256"]:
            raise ValueError(f"Source changed during audit: {entry['path']}")
    result["all_source_sha256_reverified_after_evaluation"] = True
    result["sample_artifact"] = f"Z{round(100 * metallicity):02d}_samples{output_suffix(arguments)}.npz"
    np.savez_compressed(arguments.output / result["sample_artifact"], **artifact)
    return result


def report(results):
    lines = ["# Offline composition holdout feasibility", "", results["status"], "",
             "No solver calls, new tables, runtime composition route, cache writes, commits or pushes.", "",
             "Intended primary reference: printed X=.70 SOURCE nodes, withheld from both mixtures.",
             "Off-node differences against the direct .70 biquintic are NOT accuracy estimates.", "",
             "Cubic weights (-1/6, 2/3, 2/3, -1/6) are NEGATIVE/nonconvex: stability is not inherited.", "",
             "## Input Gates And Explicit Metadata Resolution", ""]
    for case in results["metallicities"]:
        lines.append(f"- Z={case['Z']:.2f}: {case['input_gate']}.")
        if case["input_gate"] == "FAILED":
            lines.append(f"  {case['holdout_status']}.")
            for failure in case["header_grid_failures"]:
                lines.append(f"  X={failure['X']:.2f}: {json.dumps(failure, sort_keys=True)}")
        elif case.get("endpoint_metadata_override_explicitly_requested"):
            lines.extend(["  Explicit flag: `--allow-audited-endpoint-metadata`, authorized ONLY at Z=.02.",
                          "  Count/min/step, independent monotonic Q-block labels, exact block/record counts, and printed rho/T are authoritative.",
                          "  Current primary loading uses core `read_table(validate_endpoints=False)` ONLY after bounded Z=.02 authorization; independent `audit_printed_grid` must pass before constraints use rows.",
                          "  Default/non-audited reads retain `validate_endpoints=True`. This current core-keyword path supersedes the earlier independent-loader primary path; the separate strict sidecar loader remains optional.",
                          "  Only redundant display `logQ_max`/`logT_max` fields are excluded from header comparison; all other header fields/grids must match.",
                          "  No edits or resampling. X=.70 iq=520's known -1e-4 density discrepancy remains excluded by the unchanged lead mask."])
            for coordinate, entry in case["source_sha256_and_headers"].items():
                for discrepancy in entry["resolved_display_endpoint_discrepancies"]:
                    lines.append(f"  X={coordinate}: audited discrepancy {json.dumps(discrepancy, sort_keys=True)}.")
    if results.get("preserved_strict_artifacts"):
        lines.extend(["", "Initial strict failure and Z=.00 auxiliary artifacts remain unchanged; their SHA-256s are recorded in the audited JSON.",
                      "See the original `RESULTS.md`, `results.json`, and `Z02_input_gate.json` for that historical gate failure."])
    primary = next((case for case in results["metallicities"] if case["Z"] == .02 and case["input_gate"] != "FAILED"), None)
    if primary is not None:
        lines.extend(["", "## Primary Z=.02 Assessment", "",
                      f"{primary['domain']['sampled_SOURCE_07_nodes']} withheld SOURCE .70 nodes and "
                      f"{primary['domain']['saved_queries_common_composition_domain']} common saved off-node queries; "
                      f"{primary['domain']['common_material_corner_queries']} off-node queries share material corners."])
        for scheme, data in primary["schemes"].items():
            holdout = data["SOURCE_07_archived_composition_holdout"]["mixed_potential_vs_printed_07"]
            diagnostics_data = data["saved_offnode_NOT_accuracy"]["common_own_use"]["mixed_diagnostics"]
            maximum = max(entry.get("abs_max", 0) for entry in diagnostics_data["identities_balanced"].values())
            lines.append(f"- {scheme}: SOURCE-node absolute-relative p99 Pgas={100 * holdout['Pgas']['abs_p99']:.6g}%, "
                         f"Cv={100 * holdout['Cv']['abs_p99']:.6g}%; off-node balanced identity max={maximum:.6g}.")
        lines.extend(["", "P/E/S/Cv are linear in the fixed-X potential constraints at SOURCE nodes: mixing potentials does not",
                      "improve their archived composition fidelity over the same-weight physical-column control.",
                      "Nonlinear Cp/chi/Gamma/grad_ad closure differs, but is not uniformly better than that control.",
                      "No sampled positivity failures are observed here; negative cubic weights still give no general stability guarantee.",
                      "The archived holdout is scientifically useful evidence at fixed Z=.02, not fresh physical accuracy or runtime readiness."])
        lines.extend(["", "## Low-Temperature Coverage", "",
                      f"Source constraints require logT>=3.20. Native delta_logT={primary['domain']['native_temperature_step']:.2f}; "
                      "two centered radius-2 twist stages require four more temperature steps.",
                      f"The first common complete BQ node/cell begins at logT={primary['domain']['common_complete_node_min_logT']:.2f} "
                      f"(cell floor {primary['domain']['common_complete_cell_min_logT']:.2f}).",
                      "The inherited own-use criterion logT>=3.1 is not representation coverage. The interval [3.1,3.28) remains uncovered.",
                      "Correcting an upstream FreeEOS low-temperature/off-by-3.1 issue alone does not close that prototype gap.",
                      "No extrapolation, one-sided stencils, runtime fallback or coverage extension is implemented here."])
    lines.extend(["", "Z=.02 archived-field fidelity is assessed only after its explicitly audited metadata gate passes.",
                  "The Z=.00 auxiliary results cannot substitute for the requested Z=.02 primary assessment.", "",
             "## Archived source-node holdout (passed gates only)", "",
             "Relative differences use max(abs(reference),1e-300); entries are absolute p99.", "",
             "| Z | Scheme | Nodes | Field | Potential | Same-weight columns |", "|---|---|---:|---|---:|---:|"])
    for case in results["metallicities"]:
        for scheme, data in case["schemes"].items():
            holdout = data["SOURCE_07_archived_composition_holdout"]
            for field in FIELDS:
                first = holdout["mixed_potential_vs_printed_07"][field].get("abs_p99", float("nan"))
                second = holdout["same_X_weights_physical_columns_vs_printed_07_NOT_MESA_X_scheme"][field].get("abs_p99", float("nan"))
                lines.append(f"| {case['Z']:.2f} | {scheme} | {case['domain']['sampled_SOURCE_07_nodes']} | {field} | {first:.6g} | {second:.6g} |")
    lines.extend(["", "## Saved off-node comparisons (not accuracy)", ""])
    for case in results["metallicities"]:
        domain = case["domain"]
        if case["input_gate"] == "FAILED":
            lines.append(f"Z={case['Z']:.2f}: independent direct-.70 selection still has "
                         f"{domain['direct_07_saved_own_use_radius2_queries']} archived queries, but no composition queries are evaluated.")
            continue
        lines.append(f"Z={case['Z']:.2f}: {domain['direct_07_saved_own_use_radius2_queries']} direct-.70 saved queries; "
                     f"{domain['saved_queries_common_composition_domain']} in common domain; "
                     f"{domain['common_material_corner_queries']} share material corners.")
        for scheme, data in case["schemes"].items():
            diagnostics_data = data["saved_offnode_NOT_accuracy"]["common_own_use"]["mixed_diagnostics"]
            maximum = max(entry.get("abs_max", 0) for entry in diagnostics_data["identities_balanced"].values())
            failures = ", ".join(f"{field}={entry['count']}" for field, entry in diagnostics_data["positivity"].items() if entry["count"])
            lines.append(f"- {scheme}: balanced identity max {maximum:.6g}; positivity/nonfinite failures: {failures or 'none'}.")
    lines.extend(["", "## Interpretation", "",
                  "Identities test a common potential, not the historical source mixture or physical composition accuracy.",
                  "Column controls use these exact X weights on physical values, not logs and not MESA's actual X scheme.",
                  "The off-node column control additionally uses bilinear rho/T-grid interpolation, explicitly not MESA.",
                  "Gauge alignment and a common historical metal pattern are assumptions, not verified provenance.",
                  "No additive-volume mixing and no change to FreeEOS mixture physics are made.",
                  "E/S positivity is gauge dependent. Radiation is the same current-CRAD analytic reference; its remainder",
                  "is not proven pure material. Cubic weights can amplify source noise or create instabilities.",
                  "This is a bounded next-stage risk assessment, NOT a runtime-ready implementation.", "",
                  "Full extrema, masks, SHA-256s, identities, finite differences and optional X partials: results.json and Z*_samples.npz.", "",
                  "## Reproduction", "", "```bash", results["command_shell"], "```", ""])
    return "\n".join(lines)


def main():
    import shlex

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesa", type=Path, default=WORKSPACE / "mesa")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--sample-nodes", type=int, default=12000)
    parser.add_argument("--other-z", action="store_true")
    parser.add_argument("--derivative-x", action="store_true")
    parser.add_argument("--total-potential", action="store_true")
    parser.add_argument("--allow-audited-endpoint-metadata", action="store_true",
                        help="Explicit Z=.02-only resolution of redundant Q/T maximum metadata after independent printed-grid audit; never waive counts or coordinates")
    arguments = parser.parse_args()
    if arguments.output.resolve() != OUTPUT.resolve():
        parser.error("Exclusive outputs must stay in experiments/012_composition_control")
    if not 1 <= arguments.sample_nodes <= 12000:
        parser.error("Bounded node sample requires 1..12000 nodes")
    arguments.output.mkdir(parents=True, exist_ok=True)
    query_path = WORKSPACE / "freeeos_report/sample_residuals.csv"
    inherited = json.loads((WORKSPACE / "freeeos_report/audit_results.json").read_text())
    if checksum(arguments.mesa / "data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data") != inherited["sha256"]:
        raise ValueError("Withheld .02/.70 source SHA differs from the independent archived audit")
    samples = np.genfromtxt(query_path, delimiter=",", names=True)
    source_code = [Path(__file__), Path(__file__).with_name("composition.py"),
                   Path(__file__).with_name("endpoint_audit.py"),
                   Path(__file__).with_name("test_composition.py"),
                   Path(__file__).parents[1] / "freeeos_potential.py",
                   Path(__file__).parents[1] / "run_audit.py", Path(__file__).parents[1] / "validate_table.py"]
    fingerprints = {str(path): checksum(path) for path in source_code}
    result = {"status": "Bounded OFFLINE archived composition holdout; not fresh physical accuracy; NOT runtime-ready.",
              "command": [sys.executable, *sys.argv], "command_shell": shlex.join([sys.executable, *sys.argv]),
              "python": platform.python_version(), "numpy": np.__version__, "radiation_constant": float(CRAD),
              "representation": "radius2/mean analytic-radiation-reference split" if not arguments.total_potential else "radius2/mean total-potential control",
              "allow_audited_endpoint_metadata": arguments.allow_audited_endpoint_metadata,
              "source_code_sha256": fingerprints, "saved_query_path": str(query_path), "saved_query_sha256": checksum(query_path),
              "mesa_revision": subprocess.check_output(["git", "-C", str(arguments.mesa), "rev-parse", "HEAD"], text=True).strip(),
              "mesa_status": subprocess.check_output(["git", "-C", str(arguments.mesa), "status", "--short"], text=True),
              "domain_conventions": {
                  "clean": "All five native sources, finite six constraints, logT>=3.2, coordinate checks, centered radius2 both twist stages",
                  "coverage": "Source cut logT>=3.2 plus four native .02 temperature steps means centered complete BQ begins at logT=3.28; [3.1,3.28) remains uncovered even with an upstream FreeEOS fix alone",
                  "holdout": "12000 deterministic sampled lower-left native SOURCE .70 nodes, common clean complete cells, shared material node/corners and own-use; .70 only reference/mask",
                  "own_use": "logT>=3.1 and (logT>=6.5 or logQ<=-.5 or logRho<=-6.5), with inherited 1e-12 tolerance",
                  "material": "Intersection across all five sources of current Prad/printed Pgas<1e-3, -8<=logQ<=5, logT<=7.8; query corners not candidate pressure",
                  "weights": "Composition-only; all nine NODAL jets mixed BEFORE rho/T evaluation; no .70 contribution; radiation restored once because sum w=1",
                  "assumptions_NOT_verified": "Historically aligned E/S gauges and fixed metal pattern at each Z, helium=1-X-Z; filenames/headers do not verify mixture provenance",
                  "gauge": "Exact f->f+A(X)+B(X)T leaves P/Cv fixed-X unchanged; E+=A, S-=B. Nodal tau-polynomial representation of BT is not exact off-node; see gauge tests/README.",
                  "not_implemented": "Runtime selection/composition routes, additive-volume mixing, new source tables, solver or changed FreeEOS physics"},
              "metallicities": []}
    result["preserved_strict_artifacts"] = {}
    if arguments.allow_audited_endpoint_metadata:
        for name in ("results.json", "RESULTS.md", "Z02_input_gate.json", "Z04_input_gate.json", "Z00_samples.npz", "Z00_samples_audited_endpoint.npz", "run.log", "tests.log", "validation.log"):
            artifact = arguments.output / name
            if artifact.exists():
                result["preserved_strict_artifacts"][str(artifact)] = checksum(artifact)
    for metallicity in ([.02, .0, .04] if arguments.other_z else [.02]):
        print(f"Auditing fixed Z={metallicity:.2f}", flush=True)
        result["metallicities"].append(evaluate_metallicity(arguments, metallicity, samples))
    if result["metallicities"][0]["input_gate"] == "FAILED":
        result["status"] += " Primary Z=.02 HEADER/GRID GATE FAILED; its composition holdout is NOT evaluated."
    elif arguments.allow_audited_endpoint_metadata:
        result["status"] += " Primary Z=.02 evaluated with EXPLICIT independently audited display-endpoint metadata resolution."
    for path, digest in fingerprints.items():
        if checksum(path) != digest:
            raise ValueError(f"Shared/sidecar code changed during evaluation: {path}; rerun")
    if checksum(query_path) != result["saved_query_sha256"]:
        raise ValueError("Saved queries changed during evaluation")
    for path, digest in result["preserved_strict_artifacts"].items():
        if checksum(path) != digest:
            raise ValueError(f"Initial strict artifact changed: {path}")
    result_path = arguments.output / f"results{output_suffix(arguments)}.json"
    report_path = arguments.output / f"RESULTS{output_suffix(arguments)}.md"
    result_path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    report_path.write_text(report(result))
    print(f"Wrote {result_path}", flush=True)


if __name__ == "__main__":
    main()
