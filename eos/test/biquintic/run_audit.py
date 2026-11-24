"""Compare a native-only potential prototype on the saved legacy queries."""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

from freeeos_potential import (
    CRAD, EVAL_ORDERS, LN10, LOW_ORDERS, checksum, construct_twists, evaluate,
    identity_residuals, local_derivative, query_cells, read_table,
    source_constraints, thermodynamics, to_physical, valid_cells,
)


def distribution(values):
    values = np.abs(np.asarray(values))
    if values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("Cannot summarize empty or nonfinite data")
    return dict(zip(("median", "p90", "p99", "max"),
                    map(float, np.quantile(values, (.5, .9, .99, 1.)))))


def scaled_difference(first, second):
    return (first - second) / np.maximum.reduce([np.ones_like(first), abs(first), abs(second)])


def own_use(logq, logt):
    return (logt >= 3.1 - 1e-12) & ((logt >= 6.5 - 1e-12) | (logq <= -.5 + 1e-12)
                                               | (logq + 2 * logt - 12 <= -6.5 + 1e-12))


def nodal_integrability(low, valid, du, dtau, selected):
    result = {}
    for radius in (1, 2, 3):
        comparisons = {}
        for base, axis, target in (((0, 0), 0, (1, 0)), ((0, 0), 1, (0, 1)),
                                   ((1, 0), 0, (2, 0)), ((1, 0), 1, (1, 1)),
                                   ((0, 1), 0, (1, 1)), ((0, 1), 1, (0, 2))):
            derivative, accepted = local_derivative(low[base], valid, (du, dtau)[axis], axis, radius)
            accepted &= selected
            comparisons[f"d{base}_axis{axis}_vs_{target}"] = distribution(
                scaled_difference(derivative[accepted], low[target][accepted]))
        result[f"radius{radius}"] = comparisons
    return result


def finite_difference_audit(nodes, complete, logq, logt, sample_q, sample_t, state):
    results = {}
    for step in (1e-4, 1e-5, 1e-6):
        shift = step / LN10
        density_plus = thermodynamics(evaluate(nodes, complete, logq, logt, sample_q + shift, sample_t), sample_q + shift, sample_t)
        density_minus = thermodynamics(evaluate(nodes, complete, logq, logt, sample_q - shift, sample_t), sample_q - shift, sample_t)
        temperature_plus = thermodynamics(evaluate(nodes, complete, logq, logt, sample_q - 2 * shift, sample_t + shift), sample_q - 2 * shift, sample_t + shift)
        temperature_minus = thermodynamics(evaluate(nodes, complete, logq, logt, sample_q + 2 * shift, sample_t - shift), sample_q + 2 * shift, sample_t - shift)
        numerical = state.copy()
        for name in ("P", "E", "S"):
            numerical[name + "_rho"] = (density_plus[name] - density_minus[name]) / (2 * step * state["rho"])
            numerical[name + "_T"] = (temperature_plus[name] - temperature_minus[name]) / (2 * step * state["T"])
        results[str(step)] = {"identities": {name: distribution(values)
                                               for name, values in identity_residuals(numerical).items()},
                              "Cv_relative_mismatch": distribution(numerical["E_T"] / state["Cv"] - 1)}
    return results


def save_figures(output, legacy, physical, candidate, changes, sample_q, sample_t):
    os.environ.setdefault("MPLCONFIGDIR", str(output / ".mplconfig"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    figure, axes = plt.subplots(1, 3, figsize=(10.5, 3.5), constrained_layout=True)
    for axis, name in zip(axes, ("dpe", "dse", "dsp")):
        for source, label in ((legacy, "Legacy returned-function identities"),
                              (physical, "Legacy physical-column closure"),
                              (candidate, "Potential identities / column closure")):
            sorted_values = np.sort(np.maximum(abs(source[name]), 1e-17))
            axis.semilogx(sorted_values, np.arange(1, len(sorted_values) + 1) / len(sorted_values), label=label)
        axis.set_title(name)
        axis.set_xlabel("Absolute normalized residual")
        axis.grid(alpha=.25)
    axes[0].set_ylabel("Cumulative fraction")
    axes[1].legend(fontsize=6, loc="upper left")
    figure.suptitle(f"Same {len(sample_q):,} saved off-node queries: complete native stencils and own-use cuts\nPotential routes coincide analytically; fresh-source accuracy has not been measured", fontsize=10)
    for extension in ("pdf", "png"):
        figure.savefig(output / f"matched_consistency.{extension}", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(9, 3.8), constrained_layout=True)
    norm = TwoSlopeNorm(vmin=-7, vcenter=-3, vmax=1)
    for axis, values, title in ((axes[0], legacy["dpe"], "Legacy returned-function Dpe"),
                                (axes[1], changes["Cv"], "Potential vs legacy Cv (NOT source accuracy)")):
        image = axis.scatter(sample_q + 2 * sample_t - 12, sample_t,
                             c=np.log10(np.maximum(abs(values), 1e-17)), s=3, cmap="coolwarm", norm=norm, rasterized=True)
        axis.set_title(title, fontsize=9)
        axis.set_xlabel("log10 rho [g cm^-3]")
        axis.set_ylabel("log10 T [K]")
    figure.colorbar(image, ax=axes, label="log10 absolute residual / relative difference", extend="both")
    figure.suptitle("Individual matched off-node queries; no nodal or cell-median aggregation", fontsize=10)
    for extension in ("pdf", "png"):
        figure.savefig(output / f"matched_diagnostics.{extension}", dpi=180)
    plt.close(figure)


def main():
    workspace = Path(__file__).resolve().parents[4]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesa", type=Path, default=workspace / "mesa")
    parser.add_argument("--baseline", type=Path, default=workspace / "freeeos_report")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(arguments.baseline.resolve()))
    import audit_freeeos as baseline

    source = arguments.mesa / "data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
    cache = arguments.mesa / "data/eosDT_data/cache/mesa-FreeEOS_02z70x.bin"
    queries = arguments.baseline / "sample_residuals.csv"
    expected = json.loads((arguments.baseline / "audit_results.json").read_text())
    header, logq, logt, rows = read_table(source)
    if checksum(source) != expected["sha256"]:
        raise ValueError("Source checksum differs from inherited baseline")
    if not np.allclose(header[1:3], [.7, .02], rtol=0, atol=1e-12):
        raise ValueError("Expected X=.70 Z=.02 source table")
    samples = np.genfromtxt(queries, delimiter=",", names=True)
    _, coefficients = baseline.read_cache(cache, len(logq), len(logt))
    sample_q, sample_t = samples["logQ"], samples["logT"]
    values = baseline.evaluate_spline(coefficients, logq, logt, sample_q, sample_t)
    density_derivative = baseline.evaluate_spline(coefficients, logq, logt, sample_q, sample_t, (1, 0)) / LN10
    temperature_derivative = baseline.evaluate_spline(coefficients, logq, logt, sample_q, sample_t, (0, 1)) / LN10 - 2 * density_derivative
    density, temperature = 10**(sample_q + 2 * sample_t - 12), 10**sample_t
    legacy, columns, mismatch = baseline.residuals(values, density_derivative, temperature_derivative, density, temperature)
    saved_error = {name: float(np.max(abs(legacy[name] - samples[name]))) for name in legacy}
    if max(saved_error.values()) > 1e-9:
        raise ValueError("Legacy evaluator does not reproduce saved queries")
    for name, expected_key in ((legacy, "full_freeeos_eligibility_mathematical_residuals"),
                               (columns, "full_freeeos_eligibility_physical_residuals")):
        for field, result in name.items():
            actual = distribution(result[own_use(sample_q, sample_t)])
            if max(abs(actual[quantile] - expected[expected_key][field][quantile]) for quantile in actual) > 1e-10:
                raise ValueError("Inherited baseline summary changed")

    low, source_valid, source_state = source_constraints(rows, logq, logt)
    du, dtau = (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10
    schemes = {}
    for radius, route in ((1, "mean"), (2, "mean"), (3, "mean"), (2, "first"), (2, "second")):
        schemes[f"radius{radius}_{route}"] = construct_twists(low, source_valid, du, dtau, radius, route)
    density_index, temperature_index, inside = query_cells(logq, logt, sample_q, sample_t)
    common = inside.copy()
    for _, complete, _ in schemes.values():
        common &= valid_cells(complete)[density_index, temperature_index]
    eligible = common & own_use(sample_q, sample_t)
    if not np.any(eligible):
        raise ValueError("No common own-use query domain")
    selected_corners = np.zeros_like(source_valid)
    for density_end in range(2):
        for temperature_end in range(2):
            selected_corners[density_index[eligible] + density_end, temperature_index[eligible] + temperature_end] = True
    all_native = np.all(np.stack([rows[density_index + density_end, temperature_index + temperature_end, 17] == 0
                                  for density_end in range(2) for temperature_end in range(2)]), axis=0)
    pressure_rad = CRAD * temperature**4 / 3
    legacy_state = {"P": np.exp(values[:, 0]) + pressure_rad, "Pgas": np.exp(values[:, 0]),
                    "E": np.exp(values[:, 1]), "S": np.exp(values[:, 2]),
                    "Cv": values[:, baseline.FIELDS["Cv"]], "Cp": values[:, baseline.FIELDS["Cp"]],
                    "chiRho": values[:, baseline.FIELDS["chiRho"]], "chiT": values[:, baseline.FIELDS["chiT"]]}
    results = {"status": "offline prototype; not runtime or fresh-source validated",
               "provenance": {"mesa_revision": subprocess.check_output(["git", "-C", str(arguments.mesa), "rev-parse", "HEAD"], text=True).strip(),
                              "mesa_branch": subprocess.check_output(["git", "-C", str(arguments.mesa), "branch", "--show-current"], text=True).strip(),
                              "mesa_status": subprocess.check_output(["git", "-C", str(arguments.mesa), "status", "--short"], text=True),
                              "source_sha256": checksum(source), "cache_sha256": checksum(cache), "query_sha256": checksum(queries),
                              "python": sys.version, "numpy": np.__version__, "radiation_constant": float(CRAD),
                              "composition": {"X": .7, "Z": .02, "metal_pattern": "builder mass_frac_solar; direct source not yet available"},
                              "script_sha256": checksum(__file__), "preprocessing_sha256": checksum(Path(__file__).with_name("freeeos_potential.py")),
                              "command": sys.argv},
               "domain": {"original_saved_queries": len(samples), "original_own_use_queries": int(own_use(sample_q, sample_t).sum()),
                          "native_corners": int(all_native.sum()), "common_complete_stencil_queries": int(common.sum()),
                          "common_own_use_queries": int(eligible.sum()), "low_temperature_floor_logT": 3.2,
                          "printed_density_grid_disagreement_nodes": int(np.sum(abs(rows[:, :, 18] - np.log10(source_state["rho"])) > 1e-10)),
                          "printed_density_grid_max_absolute_log10_difference": float(np.max(abs(rows[:, :, 18] - np.log10(source_state["rho"])))),
                          "stencils": "centered radii 1,2,3, two-stage twists; all input nodes native, positive Cv/chiRho, logT>=3.2",
                          "sampling": "original saved independent phases .1-.9 in both coordinates, not new random queries"},
               "baseline_saved_residual_max_absolute_difference": saved_error,
               "baseline_common_own_use": {"function": {name: distribution(result[eligible]) for name, result in legacy.items()},
                                            "columns": {name: distribution(result[eligible]) for name, result in columns.items()}},
               "nodal_integrability_fd_estimates_not_source_error": nodal_integrability(low, source_valid, du, dtau, selected_corners),
               "schemes": {}}
    source_selected = selected_corners & source_valid
    source_density, source_temperature, source_pressure = source_state["rho"], source_state["T"], source_state["P"]
    mixed_pressure = source_pressure * rows[:, :, 5] / source_density
    mixed_energy = source_pressure / source_density - source_density * rows[:, :, 8]
    mixed_entropy = -source_density * source_temperature * rows[:, :, 10]
    heat_entropy = source_temperature * rows[:, :, 9]
    def source_difference(first, second):
        scale = np.maximum(abs(first), abs(second))
        valid_scale = source_selected & (scale > 0)
        difference = (first[valid_scale] - second[valid_scale]) / scale[valid_scale]
        flat_indices = np.flatnonzero(valid_scale)
        worst_flat = flat_indices[int(np.argmax(abs(difference)))]
        worst_density, worst_temperature = np.unravel_index(worst_flat, source_selected.shape)
        return {"count": int(valid_scale.sum()), "undefined_zero_scale": int((source_selected & (scale == 0)).sum()),
                **distribution(difference), "worst": {"logQ": float(logq[worst_density]), "logT": float(logt[worst_temperature]),
                                                       "first": float(first[worst_density, worst_temperature]), "second": float(second[worst_density, worst_temperature])}}
    results["redundant_source_routes_at_matched_own_use_corners"] = {
        "mixed_pressure_vs_energy": source_difference(mixed_pressure, mixed_energy),
        "mixed_pressure_vs_entropy": source_difference(mixed_pressure, mixed_entropy),
        "Cv_vs_T_S_T": source_difference(rows[:, :, 7], heat_entropy)}
    archive = {"logQ": sample_q, "logT": sample_t, "common": common, "own_use": eligible}
    central_state = None
    states = {}
    for name, (nodes, complete, routes) in schemes.items():
        native = evaluate(nodes, complete, logq, logt, sample_q[common], sample_t[common])
        state = thermodynamics(native, sample_q[common], sample_t[common])
        own = own_use(sample_q[common], sample_t[common])
        identities = identity_residuals(state)
        balanced_identities = identity_residuals(state, balanced=True)
        differences = {field: state[field] / legacy_state[field][common] - 1 for field in legacy_state}
        violations = {}
        for field in ("Cv", "P_rho", "Pgas", "E", "S"):
            failed = (~np.isfinite(state[field])) | (state[field] <= 0)
            violations[field] = {"all_common": int(failed.sum()), "own_use": int((failed & own).sum()),
                                 "examples": [{"logQ": float(sample_q[common][index]), "logT": float(sample_t[common][index]),
                                               "value": float(state[field][index])} for index in np.flatnonzero(failed)[:10]]}
        results["schemes"][name] = {"complete_nodes": int(complete.sum()), "complete_cells": int(valid_cells(complete).sum()),
                                    "identities_own_use": {field: distribution(result[own]) for field, result in identities.items()},
                                    "balanced_identities_own_use": {field: distribution(result[own]) for field, result in balanced_identities.items()},
                                    "relative_departure_from_legacy_columns_NOT_source_accuracy": {field: distribution(result[own]) for field, result in differences.items()},
                                    "positivity": violations,
                                    "twist_route_scaled_disagreement": {str(orders): distribution(scaled_difference(first[valid], second[valid]))
                                                                         for orders, (first, second, valid) in routes.items()},
                                    "twist_route_disagreement_at_matched_own_use_corners": {str(orders): distribution(scaled_difference(first[selected_corners], second[selected_corners]))
                                                                         for orders, (first, second, _) in routes.items()}}
        states[name] = state
        if name == "radius2_mean":
            central_state = state
            central_nodes, central_complete = nodes, complete
            candidate_residuals = {field: result[own] for field, result in identities.items()}
            central_differences = {field: result[own] for field, result in differences.items()}
        archive[name + "_Cv"] = state["Cv"]
    for name, state in states.items():
        if name != "radius2_mean":
            results["schemes"][name]["relative_change_from_radius2_mean_own_use"] = {
                field: distribution(state[field][own] / central_state[field][own] - 1)
                for field in ("Pgas", "E", "S", "Cv", "chiRho", "chiT")}
    own_common_indices = np.flatnonzero(eligible)[:1000]
    finite_q, finite_t = sample_q[own_common_indices], sample_t[own_common_indices]
    finite_state = thermodynamics(evaluate(central_nodes, central_complete, logq, logt, finite_q, finite_t), finite_q, finite_t)
    results["finite_difference_returned_functions"] = {"count": len(finite_q), "log_step_is_natural": True,
        "step_sweep": finite_difference_audit(central_nodes, central_complete, logq, logt, finite_q, finite_t, finite_state)}
    print(json.dumps({"domain": results["domain"], "central_scheme": results["schemes"]["radius2_mean"]}, indent=2))
    (output / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    np.savez_compressed(output / "matched_queries.npz", **archive)
    save_figures(output, {name: result[eligible] for name, result in legacy.items()},
                 {name: result[eligible] for name, result in columns.items()}, candidate_residuals,
                 central_differences, sample_q[eligible], sample_t[eligible])


if __name__ == "__main__":
    main()
