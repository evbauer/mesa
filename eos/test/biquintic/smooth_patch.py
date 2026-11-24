"""Restricted C4 quintic B-spline potential fit; not a replacement backend."""

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import BSpline
from scipy.linalg import lstsq
from scipy.sparse import csr_matrix, kron, vstack
from scipy.sparse.linalg import lsqr

from freeeos_potential import EVAL_ORDERS, LN10, LOW_ORDERS, construct_twists, evaluate, identity_residuals, read_table, source_constraints, thermodynamic_partials
from run_audit import distribution, own_use
from validate_table import relative_summary, table_state


def spline_basis(grid, query, order=0):
    knots = np.r_[np.repeat(grid[0], 6), grid[1:-1], np.repeat(grid[-1], 6)]
    count = len(knots) - 6
    spline = BSpline(knots, np.eye(count), 5, extrapolate=False)
    return spline(query, nu=order)


def evaluate_fit(coefficients, grid_u, grid_tau, sample_u, sample_tau):
    result = np.full((4, 4, len(sample_u)), np.nan)
    basis_u = {order: spline_basis(grid_u, sample_u, order) for order in range(4)}
    basis_tau = {order: spline_basis(grid_tau, sample_tau, order) for order in range(4)}
    for u_order, tau_order in EVAL_ORDERS:
        result[u_order, tau_order] = np.einsum("qi,ij,qj->q", basis_u[u_order], coefficients, basis_tau[tau_order])
    return result


def plot_slice(output, fits, grid_u, grid_tau, nodes, complete, logq, logt, workspace):
    os.environ.setdefault("MPLCONFIGDIR", str(output / ".mplconfig"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sample_q = np.linspace(2.25, 2.75, 2001)
    sample_t = np.full_like(sample_q, 6.5117106292763225)
    reference = thermodynamic_partials(evaluate(nodes, complete, logq, logt, sample_q, sample_t), sample_q, sample_t)
    sys.path.insert(0, str(workspace / "freeeos_report"))
    import audit_freeeos as baseline
    _, cache = baseline.read_cache(workspace / "mesa/data/eosDT_data/cache/mesa-FreeEOS_02z70x.bin", len(logq), len(logt))
    legacy_value = baseline.evaluate_spline(cache, logq, logt, sample_q, sample_t)[:, baseline.FIELDS["Cv"]]
    legacy_u = baseline.evaluate_spline(cache, logq, logt, sample_q, sample_t, (1, 0))[:, baseline.FIELDS["Cv"]] / LN10
    legacy_tau = baseline.evaluate_spline(cache, logq, logt, sample_q, sample_t, (0, 1))[:, baseline.FIELDS["Cv"]] / LN10
    legacy = {"Cv": legacy_value, "Cv_t": legacy_tau - 2 * legacy_u}
    figure, axes = plt.subplots(2, 1, figsize=(8, 5.5), sharex=True, constrained_layout=True)
    for label, state in [("Exact-constraint biquintic", reference), ("Legacy independent Cv column", legacy)] + [
        (name, thermodynamic_partials(evaluate_fit(coefficients, grid_u, grid_tau, sample_q * LN10, sample_t * LN10), sample_q, sample_t))
        for name, coefficients in fits.items() if name in ("native_1.0", "thermo_0.1")]:
        axes[0].plot(sample_q + 2 * sample_t - 12, state["Cv"] / 1e8, label=label)
        axes[1].plot(sample_q + 2 * sample_t - 12, state["Cv_t"] / 1e8, label=label)
    for coordinate in logq[(logq >= 2.25) & (logq <= 2.75)]:
        for axis in axes:
            axis.axvline(coordinate + 2 * sample_t[0] - 12, color="gray", alpha=.15, linewidth=.7)
    axes[0].set_ylabel("Cv / 1e8 [erg g^-1 K^-1]")
    axes[1].set_ylabel("dCv/dlnT / 1e8 (at fixed rho)")
    axes[1].set_xlabel("log10 rho [g cm^-3]")
    axes[0].legend(fontsize=8)
    figure.suptitle("Restricted C4 potential-fit control at log10 T=6.51171\nSame table inputs; curves are numerical alternatives, not independent source accuracy", fontsize=10)
    for extension in ("png", "pdf"):
        figure.savefig(output / f"response_slice.{extension}", dpi=180)
    plt.close(figure)


def fit_constraints(low, grid_u, grid_tau, curvature_weight):
    spacing_u, spacing_tau = grid_u[1] - grid_u[0], grid_tau[1] - grid_tau[0]
    center_u, center_tau = np.meshgrid(grid_u, grid_tau, indexing="ij")
    value_scale = float(np.median(abs(low[0, 0])))
    trend = np.array([np.mean(low[0, 0]), np.mean(low[1, 0]), np.mean(low[0, 1])])
    trend_grid = trend[0] + trend[1] * (center_u - np.mean(grid_u)) + trend[2] * (center_tau - np.mean(grid_tau))
    residual = low.copy()
    residual[0, 0] -= trend_grid
    residual[1, 0] -= trend[1]
    residual[0, 1] -= trend[2]
    matrices, targets = [], []
    for u_order, tau_order in LOW_ORDERS:
        first = csr_matrix(spline_basis(grid_u, grid_u, u_order))
        second = csr_matrix(spline_basis(grid_tau, grid_tau, tau_order))
        weight = spacing_u**u_order * spacing_tau**tau_order
        if u_order + tau_order == 2:
            weight *= curvature_weight
        matrices.append(kron(second, first, format="csr") * weight)
        targets.append((residual[u_order, tau_order] / value_scale * weight).ravel(order="F"))
    matrix, target = vstack(matrices, format="csr"), np.concatenate(targets)
    solution = lsqr(matrix, target, atol=2e-13, btol=2e-13, iter_lim=20000)
    if solution[1] not in (1, 2):
        raise RuntimeError(f"Constrained patch solve did not converge, istop={solution[1]}")
    coefficients = solution[0].reshape((matrices[0].shape[1] // spline_basis(grid_tau, grid_tau).shape[1],
                                        spline_basis(grid_tau, grid_tau).shape[1]), order="F") * value_scale
    greville_u = np.array([np.mean(np.r_[np.repeat(grid_u[0], 6), grid_u[1:-1], np.repeat(grid_u[-1], 6)][index + 1:index + 6])
                          for index in range(coefficients.shape[0])])
    greville_tau = np.array([np.mean(np.r_[np.repeat(grid_tau[0], 6), grid_tau[1:-1], np.repeat(grid_tau[-1], 6)][index + 1:index + 6])
                            for index in range(coefficients.shape[1])])
    coefficients += trend[0] + trend[1] * (greville_u[:, None] - np.mean(grid_u)) + trend[2] * (greville_tau[None, :] - np.mean(grid_tau))
    return coefficients, {"istop": int(solution[1]), "iterations": int(solution[2]), "residual_norm": float(solution[3]),
                          "condition_estimate": float(solution[6]), "unknowns": matrix.shape[1], "constraints": matrix.shape[0],
                          "weights": "h_u^a h_tau^b; curvature constraints multiplied by declared weight; common potential scale; affine detrend",
                          "curvature_weight": curvature_weight}


def fit_thermodynamic_constraints(low, grid_u, grid_tau, response_weight):
    matrices = {}
    for orders in LOW_ORDERS:
        matrices[orders] = kron(csr_matrix(spline_basis(grid_tau, grid_tau, orders[1])),
                                csr_matrix(spline_basis(grid_u, grid_u, orders[0])), format="csr")
    physical_temperature = matrices[0, 1] - 2 * matrices[1, 0]
    physical_curvature = matrices[0, 2] - 4 * matrices[1, 1] + 4 * matrices[2, 0]
    combinations = {"E": matrices[0, 0] - physical_temperature, "TS": -physical_temperature,
                    "P_over_rho": matrices[1, 0], "P_rho": matrices[1, 0] + matrices[2, 0],
                    "P_over_rho_chiT": matrices[1, 1] - 2 * matrices[2, 0],
                    "TCv": physical_temperature - physical_curvature}
    ft = low[0, 1] - 2 * low[1, 0]
    ftt = low[0, 2] - 4 * low[1, 1] + 4 * low[2, 0]
    target_values = {"E": low[0, 0] - ft, "TS": -ft, "P_over_rho": low[1, 0],
                     "P_rho": low[1, 0] + low[2, 0], "P_over_rho_chiT": low[1, 1] - 2 * low[2, 0],
                     "TCv": ft - ftt}
    value_scale = float(np.median(abs(low[0, 0])))
    matrices_weighted, targets = [], []
    for field, matrix in combinations.items():
        value = target_values[field].ravel(order="F")
        relative_scale = np.maximum(abs(value), .01 * np.median(abs(value)))
        importance = response_weight if field in ("P_rho", "P_over_rho_chiT", "TCv") else 1
        weights = importance * value_scale / relative_scale
        matrices_weighted.append(matrix.multiply(weights[:, None]))
        targets.append(value / value_scale * weights)
    matrix, target = vstack(matrices_weighted, format="csr"), np.concatenate(targets)
    column_scale = np.sqrt(np.asarray(matrix.power(2).sum(axis=0)).ravel())
    column_scale = np.maximum(column_scale, 1)
    equilibrated = matrix.multiply(1 / column_scale).toarray()
    solution, _, rank, _ = lstsq(equilibrated, target, cond=1e-13, lapack_driver="gelsy")
    if rank != matrix.shape[1]:
        raise RuntimeError(f"Thermodynamic patch solve rank deficient, rank={rank}")
    shape = (spline_basis(grid_u, grid_u).shape[1], spline_basis(grid_tau, grid_tau).shape[1])
    coefficients = (solution / column_scale).reshape(shape, order="F") * value_scale
    return coefficients, {"method": "dense rank-revealing QR after column equilibration", "rank": int(rank), "iterations": 0,
                          "residual_norm": float(np.linalg.norm(equilibrated @ solution - target)),
                          "unknowns": matrix.shape[1], "constraints": matrix.shape[0],
                          "weights": "relative linear E,TS,P/rho,P_rho,P/rho*chiT,T*Cv; declared response multiplier; column equilibration",
                          "response_weight": response_weight}


def main():
    workspace = Path(__file__).resolve().parents[4]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesa", type=Path, default=workspace / "mesa")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    _, logq, logt, rows = read_table(arguments.mesa / "data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data")
    low, valid, _ = source_constraints(rows, logq, logt)
    selected_q = (logq >= 1.8) & (logq <= 2.9)
    selected_t = (logt >= 6.5 - 1e-12) & (logt <= 6.9 + 1e-12)
    if np.flatnonzero(selected_q)[0] % 2:
        selected_q[np.flatnonzero(selected_q)[0]] = False
    if np.flatnonzero(selected_t)[0] % 2:
        selected_t[np.flatnonzero(selected_t)[0]] = False
    patch_q, patch_t = logq[selected_q], logt[selected_t]
    patch_low = low[:, :, selected_q][:, :, :, selected_t]
    patch_valid = valid[selected_q][:, selected_t]
    if not np.all(patch_valid):
        raise ValueError("Restricted smooth patch must have all native valid nodes")
    grid_u, grid_tau = patch_q * LN10, patch_t * LN10
    node_q, node_t = np.meshgrid(patch_q, patch_t, indexing="ij")
    node_q, node_t = node_q.ravel(), node_t.ravel()
    node_rows = rows[selected_q][:, selected_t].reshape(-1, 22)
    generator = np.random.default_rng(20261031)
    sample_q = generator.uniform(patch_q[1], patch_q[-2], 5000)
    sample_t = generator.uniform(patch_t[1], patch_t[-2], 5000)
    nodes, complete, _ = construct_twists(low, valid, (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10)
    reference = thermodynamic_partials(evaluate(nodes, complete, logq, logt, sample_q, sample_t), sample_q, sample_t)
    result = {"purpose": "restricted C4 potential fitting tradeoff, not production EOS or verified physical accuracy",
              "domain": {"logQ": [float(patch_q[0]), float(patch_q[-1])], "logT": [float(patch_t[0]), float(patch_t[-1])],
                         "nodes": len(node_q), "offnode_queries": len(sample_q), "boundary_caveat": "fit-patch ends are not external EOS or neighboring-patch blends"},
              "basis": "degree5 tensor B-splines with simple interior knots at retained native coordinates (C4), six soft nodal constraints",
              "candidates": {}}
    configurations = [(f"native_{weight}", fit_constraints, weight) for weight in (.1, 1., 10.)]
    configurations += [(f"thermo_{weight}", fit_thermodynamic_constraints, weight) for weight in (.01, .1, 1.)]
    result["failed_attempts"] = [{"method": "thermodynamic-weights column-equilibrated LSQR", "response_weight": .01,
                                 "iteration_limit": 40000, "failure": "istop7 iteration limit; changed restricted patch to dense rank-revealing QR"}]
    fits = {}
    for name, fitter, importance in configurations:
        coefficients, solve = fitter(patch_low, grid_u, grid_tau, importance)
        fits[name] = coefficients
        native = evaluate_fit(coefficients, grid_u, grid_tau, sample_q * LN10, sample_t * LN10)
        state = thermodynamic_partials(native, sample_q, sample_t)
        node_native = evaluate_fit(coefficients, grid_u, grid_tau, node_q * LN10, node_t * LN10)
        node_state = thermodynamic_partials(node_native, node_q, node_t)
        candidate = {"solve": solve, "archived_node_relative_departure": relative_summary(node_state, table_state(node_rows), node_q, node_t),
                     "offnode_difference_vs_exact_constraint_biquintic_NOT_source_accuracy": relative_summary(state, reference, sample_q, sample_t),
                     "identities": {field: distribution(values) for field, values in identity_residuals(state).items()},
                     "positivity": {field: int(np.sum((state[field] <= 0) | ~np.isfinite(state[field]))) for field in ("Cv", "P_rho", "Pgas", "E", "S")},
                     "response_jacobian_difference_scaled_by_response_or_partial": {}}
        for field in ("Cv", "Cp", "chiRho", "chiT", "gamma1", "grad_ad"):
            for suffix in ("_x", "_t"):
                scale = np.maximum.reduce([abs(reference[field]), abs(reference[field + suffix]), np.ones(len(sample_q))])
                candidate["response_jacobian_difference_scaled_by_response_or_partial"][field + suffix] = distribution((state[field + suffix] - reference[field + suffix]) / scale)
        result["candidates"][name] = candidate
    retained_q, retained_t = patch_q[::2], patch_t[::2]
    retained_low = patch_low[:, :, ::2, ::2]
    withheld_q_index, withheld_t_index = np.meshgrid(np.arange(1, len(patch_q) - 2, 2), np.arange(1, len(patch_t) - 2, 2), indexing="ij")
    withheld_q, withheld_t = patch_q[withheld_q_index.ravel()], patch_t[withheld_t_index.ravel()]
    withheld_rows = rows[selected_q][:, selected_t][withheld_q_index.ravel(), withheld_t_index.ravel()]
    result["stride2_withheld_patch"] = {"count": len(withheld_q), "sampling": "source native nodes withheld in both axes; coarse patch fit from retained source constraints only; not source accuracy", "candidates": {}}
    for name, fitter, importance in configurations:
        coefficients, solve = fitter(retained_low, retained_q * LN10, retained_t * LN10, importance)
        native = evaluate_fit(coefficients, retained_q * LN10, retained_t * LN10, withheld_q * LN10, withheld_t * LN10)
        state = thermodynamic_partials(native, withheld_q, withheld_t)
        result["stride2_withheld_patch"]["candidates"][name] = {
            "solve": solve, "archived_withheld_node_relative_departure": relative_summary(state, table_state(withheld_rows), withheld_q, withheld_t)}
    coarse_q, coarse_t, coarse_rows = logq[::2], logt[::2], rows[::2, ::2]
    coarse_low, coarse_valid, _ = source_constraints(coarse_rows, coarse_q, coarse_t)
    coarse_nodes, coarse_complete, _ = construct_twists(coarse_low, coarse_valid, (coarse_q[1] - coarse_q[0]) * LN10, (coarse_t[1] - coarse_t[0]) * LN10)
    coarse_state = thermodynamic_partials(evaluate(coarse_nodes, coarse_complete, coarse_q, coarse_t, withheld_q, withheld_t), withheld_q, withheld_t)
    result["stride2_withheld_patch"]["biquintic_control"] = relative_summary(coarse_state, table_state(withheld_rows), withheld_q, withheld_t)
    result["command"] = sys.argv
    result["limitations"] = "same patch data used for fitting and nodal validation; no fit-accuracy extrapolation or source oracle; internal simple knots C4, patch boundaries not blended"
    plot_slice(output, fits, grid_u, grid_tau, nodes, complete, logq, logt, workspace)
    (output / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({name: {"iterations": data["solve"]["iterations"], "nodal_Cv_p99": data["archived_node_relative_departure"]["Cv"]["p99"],
                            "offnode_Cv_p99": data["offnode_difference_vs_exact_constraint_biquintic_NOT_source_accuracy"]["Cv"]["p99"]}
                      for name, data in result["candidates"].items()}, indent=2))


if __name__ == "__main__":
    main()
