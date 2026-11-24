"""Bounded multislice reconciliation/holdout checks of the lead C4 soft fit."""

import os

os.environ["OPENBLAS_NUM_THREADS"] = "8"
os.environ["OMP_NUM_THREADS"] = "8"
os.environ["MKL_NUM_THREADS"] = "8"

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import multiprocessing
from pathlib import Path
import resource
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from freeeos_potential import (
    CRAD, EVAL_ORDERS, LN10, LOW_ORDERS, checksum, construct_twists, evaluate,
    identity_residuals, query_cells, read_table, source_constraints,
    thermodynamic_partials, valid_cells,
)
from smooth_patch import evaluate_fit, fit_constraints, fit_thermodynamic_constraints, spline_basis
from validate_table import material_nodes, table_state


PATCHES = {"ionization": ((-3.55, -2.53), (3.5, 4.0)),
           "helium_ionization": ((-2.5, -1.48), (4.4, 4.9)),
           "dense_tail": ((1.85, 2.87), (6.5, 6.9)),
           "active_chiT": ((.35, 1.4), (6.5, 7.0))}
CONFIGURATIONS = (("thermo_0.01", "thermo", .01), ("thermo_0.1", "thermo", .1),
                  ("thermo_1", "thermo", 1.), ("native_1", "native", 1.))
FIELDS = ("P", "Pgas", "E", "S", "Cv", "Cp", "chiRho", "chiT", "gamma1", "gamma3", "grad_ad")
RESPONSES = ("Cv", "Cp", "chiRho", "chiT", "gamma1", "gamma3", "grad_ad")
JACOBIANS = tuple(field + suffix for field in RESPONSES for suffix in ("_x", "_t"))


def distribution(values):
    values = np.abs(np.asarray(values))
    finite = values[np.isfinite(values)]
    result = {"count": int(values.size), "finite_count": int(finite.size), "nonfinite_count": int(values.size - finite.size)}
    if finite.size:
        result.update(dict(zip(("median", "p99", "max"), map(float, np.quantile(finite, (.5, .99, 1))))))
    else:
        result.update({key: None for key in ("median", "p99", "max")})
    return result


def choose_patch(quality, logq, logt, bounds):
    qbounds, tbounds = bounds
    density_indices = np.flatnonzero((logq >= qbounds[0] - 1e-12) & (logq <= qbounds[1] + 1e-12))
    temperature_indices = np.flatnonzero((logt >= tbounds[0] - 1e-12) & (logt <= tbounds[1] + 1e-12))
    best = None
    for start_q in density_indices[density_indices % 2 == 0]:
        for count_q in range(17, 36, 2):
            if start_q + count_q - 1 > density_indices[-1]:
                continue
            for start_t in temperature_indices[temperature_indices % 2 == 0]:
                for count_t in range(17, 36, 2):
                    if start_t + count_t - 1 > temperature_indices[-1] or (count_q // 2) * (count_t // 2) < 144:
                        continue
                    if not quality[start_q:start_q + count_q, start_t:start_t + count_t].all():
                        continue
                    score = count_q * count_t
                    if best is None or score > best[0]:
                        best = score, int(start_q), int(start_t), count_q, count_t
    if best is None:
        raise ValueError("No all-valid clean-donor rectangle with >=144 odd/odd holdouts inside requested bounds")
    return best[1:]


def reference_state(rows):
    result = table_state(rows)
    result["P"] = result["Pgas"] + CRAD * (10**rows[..., 0])**4 / 3
    return result


def comparison(candidate, reference, sample_q, sample_t, fields=FIELDS, jacobians=False):
    result = {}
    for field in fields:
        if jacobians:
            base = field[:-2]
            scale = np.maximum(np.maximum(abs(reference[field]), abs(reference[base])), np.finfo(float).tiny)
        else:
            scale = np.maximum(abs(reference[field]), np.finfo(float).tiny)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            difference = (candidate[field] - reference[field]) / scale
        metric = distribution(difference)
        finite = np.isfinite(difference)
        if finite.any():
            index = int(np.argmax(np.where(finite, abs(difference), -1)))
            metric["worst"] = {"logQ": float(sample_q[index]), "logT": float(sample_t[index]),
                               "candidate": float(candidate[field][index]), "reference": float(reference[field][index]),
                               "scale": float(scale[index]), "signed_difference": float(difference[index])}
        result[field] = metric
    return result


def positivity(state, sample_q, sample_t):
    result = {}
    for field in ("P", "Pgas", "E", "S", "Cv", "Cp", "P_rho", "chiRho"):
        failed = ~np.isfinite(state[field]) | (state[field] <= 0)
        result[field] = {"count": int(failed.sum()), "total": len(sample_q), "examples": [
            {"logQ": float(sample_q[index]), "logT": float(sample_t[index]),
             "value": float(state[field][index]) if np.isfinite(state[field][index]) else None}
            for index in np.flatnonzero(failed)[:5]]}
    return result


def state_subset(state, selected):
    return {field: values[selected] for field, values in state.items()}


def with_tails(candidate, reference, sample_q, sample_t, grid_q, grid_t, margin):
    boundary = ((sample_q < grid_q[0] + margin * (grid_q[1] - grid_q[0])) |
                (sample_q > grid_q[-1] - margin * (grid_q[1] - grid_q[0])) |
                (sample_t < grid_t[0] + margin * (grid_t[1] - grid_t[0])) |
                (sample_t > grid_t[-1] - margin * (grid_t[1] - grid_t[0])))
    return {label: {"count": int(selected.sum()), "relative_departure": comparison(
        state_subset(candidate, selected), state_subset(reference, selected), sample_q[selected], sample_t[selected])}
        for label, selected in (("all", np.ones(len(sample_q), bool)), ("boundary_band", boundary), ("interior", ~boundary))}


def native_reconciliation(native, low, sample_q, sample_t):
    result = {}
    for orders in LOW_ORDERS:
        reference = low[orders].ravel()
        scale = np.maximum(abs(reference), .01 * np.median(abs(reference)))
        scale = np.maximum(scale, np.finfo(float).tiny)
        difference = (native[orders] - reference) / scale
        index = int(np.argmax(abs(difference)))
        result[str(orders)] = {**distribution(difference), "worst": {
            "logQ": float(sample_q[index]), "logT": float(sample_t[index]),
            "signed_difference": float(difference[index])}}
    return result


def sampled_faces(grid_q, grid_t, count, seed):
    generator = np.random.default_rng(seed)
    result = {}
    for axis, direction in ((0, "u"), (1, "tau")):
        boundary = generator.integers(1, len(grid_q if axis == 0 else grid_t) - 1, count)
        tangential_grid = grid_t if axis == 0 else grid_q
        tangential = generator.integers(0, len(tangential_grid) - 1, count)
        fraction = generator.uniform(.1, .9, count)
        position = tangential_grid[tangential] + fraction * (tangential_grid[1] - tangential_grid[0])
        result[direction] = (grid_q[boundary], position) if axis == 0 else (position, grid_t[boundary])
    return result


def fourth_partials(coefficients, grid_u, grid_tau, sample_u, sample_tau):
    return {str((density_order, 4 - density_order)): np.einsum(
        "qi,ij,qj->q", spline_basis(grid_u, sample_u, density_order), coefficients,
        spline_basis(grid_tau, sample_tau, 4 - density_order)) for density_order in range(5)}


def continuity(coefficients, grid_q, grid_t, count, seed):
    grid_u, grid_tau = grid_q * LN10, grid_t * LN10
    result = {"theorem": "degree5 simple interior knots: potential C4, not generally C5; response Jacobians continuous away from singular denominators",
              "method": "finite one-sided native natural-log steps; gaps include O(step) variation, not exact same-coordinate jumps; two-step linear intercept also reported", "faces": {}}
    for axis, (direction, (sample_q, sample_t)) in enumerate(sampled_faces(grid_q, grid_t, count, seed).items()):
        spacing = (grid_u[1] - grid_u[0]) if axis == 0 else (grid_tau[1] - grid_tau[0])
        steps, raw_gaps, scales = {}, [], []
        for fraction in (1e-3, 1e-4, 1e-5):
            step = fraction * spacing
            left_q, right_q = sample_q.copy(), sample_q.copy()
            left_t, right_t = sample_t.copy(), sample_t.copy()
            if axis == 0:
                left_q -= step / LN10
                right_q += step / LN10
            else:
                left_t -= step / LN10
                right_t += step / LN10
            left = thermodynamic_partials(evaluate_fit(coefficients, grid_u, grid_tau, left_q * LN10, left_t * LN10), left_q, left_t)
            right = thermodynamic_partials(evaluate_fit(coefficients, grid_u, grid_tau, right_q * LN10, right_t * LN10), right_q, right_t)
            fourth_left = fourth_partials(coefficients, grid_u, grid_tau, left_q * LN10, left_t * LN10)
            fourth_right = fourth_partials(coefficients, grid_u, grid_tau, right_q * LN10, right_t * LN10)
            gaps, scale = {}, {}
            for field in JACOBIANS:
                gaps[field] = left[field] - right[field]
                scale[field] = np.maximum.reduce([abs(left[field]), abs(right[field]), abs(left[field[:-2]]),
                                                 abs(right[field[:-2]]), np.full(count, np.finfo(float).tiny)])
            for field in fourth_left:
                key = "potential_fourth_" + field
                gaps[key] = fourth_left[field] - fourth_right[field]
                scale[key] = np.maximum.reduce([abs(fourth_left[field]), abs(fourth_right[field]),
                    np.full(count, max(np.median(abs(fourth_left[field])) * 1e-12, np.finfo(float).tiny))])
            steps[str(fraction)] = {"natural_step": float(step), "scaled_gap": {field: distribution(values / scale[field]) for field, values in gaps.items()}}
            raw_gaps.append(gaps)
            scales.append(scale)
        intercept = {field: (10 * raw_gaps[-1][field] - raw_gaps[-2][field]) / 9 for field in raw_gaps[-1]}
        result["faces"][direction] = {"count": count, "step_sweep": steps,
            "estimated_zero_step_intercept_scaled": {field: distribution(values / scales[-1][field]) for field, values in intercept.items()},
            "roundoff_caveat": "intercept cancellation and large potential coefficients can set a numerical floor; theorem is separate from sampled arithmetic"}
    return result


def process_patch(task):
    started = time.perf_counter()
    name, selection, source_path, output_path, query_count, face_count, seed, fingerprints = task
    _, logq, logt, rows = read_table(source_path)
    low, valid, _ = source_constraints(rows, logq, logt)
    nodes, complete, _ = construct_twists(low, valid, (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10)
    start_q, start_t, count_q, count_t = selection
    selected_q, selected_t = np.arange(start_q, start_q + count_q), np.arange(start_t, start_t + count_t)
    patch_q, patch_t = logq[selected_q], logt[selected_t]
    patch_rows = rows[np.ix_(selected_q, selected_t)]
    patch_low = low[:, :, selected_q][:, :, :, selected_t]
    if not (valid[np.ix_(selected_q, selected_t)] & material_nodes(patch_rows, patch_q, patch_t)).all():
        raise ValueError("Patch provenance/material changed after selection")
    coarse_q, coarse_t, coarse_rows = logq[::2], logt[::2], rows[::2, ::2]
    coarse_low, coarse_valid, _ = source_constraints(coarse_rows, coarse_q, coarse_t)
    coarse_nodes, coarse_complete, _ = construct_twists(coarse_low, coarse_valid, (coarse_q[1] - coarse_q[0]) * LN10, (coarse_t[1] - coarse_t[0]) * LN10)
    retained_q, retained_t, retained_low = patch_q[::2], patch_t[::2], patch_low[:, :, ::2, ::2]
    local_coarse_nodes, local_coarse_complete, _ = construct_twists(retained_low, np.ones(retained_low.shape[2:], bool),
        (retained_q[1] - retained_q[0]) * LN10, (retained_t[1] - retained_t[0]) * LN10)
    node_q, node_t = np.meshgrid(patch_q, patch_t, indexing="ij")
    node_q, node_t = node_q.ravel(), node_t.ravel()
    node_reference = reference_state(patch_rows.reshape(-1, 22))
    generator = np.random.default_rng(seed)
    density_index = generator.integers(0, count_q - 1, query_count)
    temperature_index = generator.integers(0, count_t - 1, query_count)
    sample_q = patch_q[density_index] + generator.uniform(.1, .9, query_count) * (patch_q[1] - patch_q[0])
    sample_t = patch_t[temperature_index] + generator.uniform(.1, .9, query_count) * (patch_t[1] - patch_t[0])
    biquintic = thermodynamic_partials(evaluate(nodes, complete, logq, logt, sample_q, sample_t), sample_q, sample_t)
    withheld_density, withheld_temperature = np.meshgrid(np.arange(1, count_q, 2), np.arange(1, count_t, 2), indexing="ij")
    withheld_q, withheld_t = patch_q[withheld_density.ravel()], patch_t[withheld_temperature.ravel()]
    withheld_rows = patch_rows[withheld_density.ravel(), withheld_temperature.ravel()]
    withheld_reference = reference_state(withheld_rows)
    if len(withheld_q) < 144 or start_q % 2 or start_t % 2:
        raise AssertionError("Holdout count/alignment contract violated")
    coarse_biquintic = thermodynamic_partials(evaluate(coarse_nodes, coarse_complete, coarse_q, coarse_t, withheld_q, withheld_t), withheld_q, withheld_t)
    local_density, local_temperature, _ = query_cells(retained_q, retained_t, withheld_q, withheld_t)
    local_subset = valid_cells(local_coarse_complete)[local_density, local_temperature]
    local_biquintic = thermodynamic_partials(evaluate(local_coarse_nodes, local_coarse_complete, retained_q, retained_t,
        withheld_q[local_subset], withheld_t[local_subset]), withheld_q[local_subset], withheld_t[local_subset])
    result = {"patch": name, "requested_bounds": PATCHES[name], "selected_indices": selection,
        "actual_bounds": {"logQ": [float(patch_q[0]), float(patch_q[-1])], "logT": [float(patch_t[0]), float(patch_t[-1])]},
        "provenance": fingerprints, "domain": {"all_source_nodes_native_valid_material": True, "full_global_and_coarse_donors_complete_clean": True,
            "node_count": len(node_q), "queries": query_count, "withheld_count": len(withheld_q), "local_coarse_common_withheld_count": int(local_subset.sum()),
            "sampling": "independent fractions [.1,.9] strictly inside both native coordinates, one random cell per draw", "seed": seed},
        "coarse_biquintic_global_donor_control": {"archived_withheld_departure": with_tails(coarse_biquintic, withheld_reference, withheld_q, withheld_t, retained_q, retained_t, 2),
            "positivity": positivity(coarse_biquintic, withheld_q, withheld_t)},
        "coarse_biquintic_patch_only_control": {"relative_departure": comparison(local_biquintic, state_subset(withheld_reference, local_subset), withheld_q[local_subset], withheld_t[local_subset]),
            "positivity": positivity(local_biquintic, withheld_q[local_subset], withheld_t[local_subset])},
        "global_control_on_patch_only_common_subset": {"relative_departure": comparison(state_subset(coarse_biquintic, local_subset), state_subset(withheld_reference, local_subset), withheld_q[local_subset], withheld_t[local_subset]),
            "global_vs_local_donor_departure": comparison(local_biquintic, state_subset(coarse_biquintic, local_subset), withheld_q[local_subset], withheld_t[local_subset])},
        "controls_caveat": "global biquintic uses clean donor rows outside fit patch; patch-only coarse control uses the SAME retained rows but stencil erosion reduces its common holdout set",
        "full_biquintic_positivity": positivity(biquintic, sample_q, sample_t), "candidates": {}}
    archived = {"patch_logQ": patch_q, "patch_logT": patch_t, "logQ": sample_q, "logT": sample_t,
                "withheld_logQ": withheld_q, "withheld_logT": withheld_t, "local_coarse_common": local_subset,
                "metrics": np.array(FIELDS + JACOBIANS), "biquintic_values": np.stack([biquintic[field] for field in FIELDS + JACOBIANS])}
    for candidate_name, kind, weight in CONFIGURATIONS:
        fitter = fit_thermodynamic_constraints if kind == "thermo" else fit_constraints
        candidate = {"full": {}, "coarse": {}}
        for resolution, fit_low, grid_q, grid_t in (("full", patch_low, patch_q, patch_t), ("coarse", retained_low, retained_q, retained_t)):
            fit_started = time.perf_counter()
            try:
                coefficients, solve = fitter(fit_low, grid_q * LN10, grid_t * LN10, weight)
                solve["elapsed_seconds"] = time.perf_counter() - fit_started
                diagnostics = {"status": "solved", "solve": solve}
                if resolution == "full":
                    native = evaluate_fit(coefficients, grid_q * LN10, grid_t * LN10, sample_q * LN10, sample_t * LN10)
                    state = thermodynamic_partials(native, sample_q, sample_t)
                    nodal_native = evaluate_fit(coefficients, grid_q * LN10, grid_t * LN10, node_q * LN10, node_t * LN10)
                    nodal_state = thermodynamic_partials(nodal_native, node_q, node_t)
                    diagnostics.update({"archived_node_reconciliation_NOT_accuracy": with_tails(nodal_state, node_reference, node_q, node_t, patch_q, patch_t, 2),
                        "native_six_constraint_reconciliation": native_reconciliation(nodal_native, patch_low, node_q, node_t),
                        "offnode_departure_vs_biquintic_NOT_accuracy": with_tails(state, biquintic, sample_q, sample_t, patch_q, patch_t, 2),
                        "response_jacobian_departure_field_scaled": comparison(state, biquintic, sample_q, sample_t, JACOBIANS, True),
                        "query_positivity": positivity(state, sample_q, sample_t), "nodal_positivity": positivity(nodal_state, node_q, node_t),
                        "identities": {field: distribution(values) for field, values in identity_residuals(state).items()}})
                    archived[candidate_name + "_values"] = np.stack([state[field] for field in FIELDS + JACOBIANS])
                else:
                    state = thermodynamic_partials(evaluate_fit(coefficients, grid_q * LN10, grid_t * LN10, withheld_q * LN10, withheld_t * LN10), withheld_q, withheld_t)
                    diagnostics.update({"archived_withheld_departure_NOT_accuracy": with_tails(state, withheld_reference, withheld_q, withheld_t, retained_q, retained_t, 2),
                        "matched_patch_only_donor_subset_departure": comparison(state_subset(state, local_subset), state_subset(withheld_reference, local_subset), withheld_q[local_subset], withheld_t[local_subset]),
                        "withheld_positivity": positivity(state, withheld_q, withheld_t)})
                    archived[candidate_name + "_withheld_values"] = np.stack([state[field] for field in FIELDS])
                diagnostics["continuity"] = continuity(coefficients, grid_q, grid_t, face_count, seed + 1)
                diagnostics["total_validation_seconds"] = time.perf_counter() - fit_started
                archived[candidate_name + "_" + resolution + "_coefficients"] = coefficients
                candidate[resolution] = diagnostics
                print(f"{name}: {candidate_name} {resolution} solved in {solve['elapsed_seconds']:.1f}s", flush=True)
            except (RuntimeError, ValueError, np.linalg.LinAlgError) as error:
                candidate[resolution] = {"status": "failed", "error": str(error), "elapsed_seconds": time.perf_counter() - fit_started}
                print(f"{name}: {candidate_name} {resolution} FAILED: {error}", flush=True)
        result["candidates"][candidate_name] = candidate
    result["elapsed_seconds"] = time.perf_counter() - started
    result["peak_rss_kib_linux"] = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    result["limitations"] = "table-only soft reconciliation, no verified gauge or source accuracy; no patch blending or global C4 across patch ends; sampling is not a positivity bound"
    for path, expected in fingerprints["files"].items():
        if checksum(path) != expected:
            raise RuntimeError("Read-only input changed during worker execution")
    output = Path(output_path) / name
    output.mkdir(parents=True, exist_ok=True)
    (output / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    np.savez_compressed(output / "fits_and_queries.npz", **archived)
    return result


def summary(results):
    lines = ["# Multislice C4 soft-potential validation", "", "Table reconciliation/holdout fidelity, NOT physical accuracy. Smoothing is not necessarily truth.", "",
             "Degree5 simple interior knots give C4 within each restricted patch; patch ends are not blended. Four workers x eight BLAS threads, no fallback/global smoother fit.", "",
             "| Patch | Selected q / logT | Nodes | Withheld / local-donor common |", "|---|---|---:|---:|"]
    for name, patch in results["patches"].items():
        if "error" in patch:
            lines.append(f"| {name} | FAILED: {patch['error']} | — | — |")
            continue
        bounds, domain = patch["actual_bounds"], patch["domain"]
        lines.append(f"| {name} | {bounds['logQ']} / {bounds['logT']} | {domain['node_count']} | {domain['withheld_count']} / {domain['local_coarse_common_withheld_count']} |")
    lines.extend(["", "Ionization bounds are explicitly shrunk to remove excluded source/donor nodes; starts are globally even and dimensions are odd, 17–35 nodes per axis.", "",
                  "| Patch / fit | Full nodal Cv p99 | Offnode Cv p99 vs biquintic | Stride2 withheld Cv p99 | Withheld Cv max |", "|---|---:|---:|---:|---:|"])
    for name, patch in results["patches"].items():
        if "error" in patch:
            continue
        for candidate_name, candidate in patch["candidates"].items():
            full, coarse = candidate["full"], candidate["coarse"]
            node = full.get("archived_node_reconciliation_NOT_accuracy", {}).get("all", {}).get("relative_departure", {}).get("Cv", {}).get("p99")
            offnode = full.get("offnode_departure_vs_biquintic_NOT_accuracy", {}).get("all", {}).get("relative_departure", {}).get("Cv", {}).get("p99")
            withheld = coarse.get("archived_withheld_departure_NOT_accuracy", {}).get("all", {}).get("relative_departure", {}).get("Cv", {})
            lines.append(f"| {name} / {candidate_name} | {node} | {offnode} | {withheld.get('p99')} | {withheld.get('max')} |")
        control = patch["coarse_biquintic_global_donor_control"]["archived_withheld_departure"]["all"]["relative_departure"]["Cv"]
        lines.append(f"| {name} / coarse biquintic (global donors) | — | — | {control['p99']} | {control['max']} |")
    lines.extend(["", "## Interpretation and next checks", "",
        "- Increasing the thermodynamic response multiplier trades state reconciliation against P_rho/chiT/Cv reconciliation; the native-weight control uses spacing-scaled potential constraints, not the same relative objective.",
        "- Withheld-node comparisons use identical native odd/odd rows and coarse-only training. Global coarse biquintic can evaluate boundary holdouts using clean outside donors; patch-only stencils cannot. Compare the SAME common interior subset, not its p99 against the global all-holdout p99. Both donor controls on that subset and their differences are recorded.",
        "- Inspect boundary/interior p99/max and positivity before selecting a weight. Better nodal agreement or C4 continuity alone does not imply better holdout fidelity, stability, or solver behavior.",
        "- Finite one-sided natural steps test convergence toward continuity; they include smooth O(step) variation and a roundoff-sensitive zero-step intercept. C4 is a basis property, not a claim of unchanged source curvature or C4 patch blending.",
        "- Source flags, masks, radiation and Gamma3 conventions remain unchanged. Gauge alignment, source accuracy, compositions, active fractions and Newton inversion are separate gates.",
        "", "Full distributions, worst coordinates, solve weights/rank/residuals, nonfinite accounting, timings/RSS and coefficient/query archives are in each patch directory."])
    return "\n".join(lines) + "\n"


def main():
    workspace = Path(__file__).resolve().parents[5]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--queries", type=int, default=5000)
    parser.add_argument("--faces", type=int, default=128)
    parser.add_argument("--patches", nargs="+", choices=tuple(PATCHES), default=list(PATCHES))
    arguments = parser.parse_args()
    if not 1 <= arguments.workers <= 4 or arguments.queries < 1 or arguments.faces < 1:
        parser.error("Require 1–4 workers and positive query/face counts")
    source_path = workspace / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
    output = workspace / "freeeos_biquintic_progress/experiments/015_smooth_validation"
    code = Path(__file__).resolve().parents[1]
    inputs = [source_path, code / "freeeos_potential.py", code / "smooth_patch.py", code / "validate_table.py", Path(__file__)]
    fingerprints = {"files": {str(path): checksum(path) for path in inputs}, "python": sys.version, "numpy": np.__version__, "command": sys.argv}
    _, logq, logt, rows = read_table(source_path)
    low, valid, _ = source_constraints(rows, logq, logt)
    _, complete, _ = construct_twists(low, valid, (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10)
    coarse_q, coarse_t, coarse_rows = logq[::2], logt[::2], rows[::2, ::2]
    coarse_low, coarse_valid, _ = source_constraints(coarse_rows, coarse_q, coarse_t)
    _, coarse_complete, _ = construct_twists(coarse_low, coarse_valid, (coarse_q[1] - coarse_q[0]) * LN10, (coarse_t[1] - coarse_t[0]) * LN10)
    mesh_q, mesh_t = np.meshgrid(logq, logt, indexing="ij")
    full_density, full_temperature, _ = query_cells(logq, logt, mesh_q, mesh_t)
    coarse_density, coarse_temperature, _ = query_cells(coarse_q, coarse_t, mesh_q, mesh_t)
    quality = valid & material_nodes(rows, logq, logt)
    quality &= valid_cells(complete)[full_density, full_temperature] & valid_cells(coarse_complete)[coarse_density, coarse_temperature]
    results = {"status": "bounded table-only multislice C4 soft-fit validation", "provenance": fingerprints,
        "parallelism": {"workers": arguments.workers, "BLAS_threads_each": 8, "maximum_solver_threads": arguments.workers * 8},
        "selection": "largest all-valid material clean-full/coarse-donor rectangle inside requested bounds; even global starts; odd 17–35 dimensions; >=144 odd/odd holdouts", "patches": {}}
    tasks = []
    for position, name in enumerate(arguments.patches):
        try:
            selection = choose_patch(quality, logq, logt, PATCHES[name])
            print(name, "selection", selection, flush=True)
            tasks.append((name, selection, str(source_path), str(output), arguments.queries, arguments.faces, 2026103100 + position, fingerprints))
        except ValueError as error:
            results["patches"][name] = {"error": str(error)}
    started = time.perf_counter()
    with ProcessPoolExecutor(max_workers=arguments.workers, mp_context=multiprocessing.get_context("spawn")) as executor:
        futures = {executor.submit(process_patch, task): task[0] for task in tasks}
        for future in as_completed(futures):
            name = futures[future]
            try:
                results["patches"][name] = future.result()
            except Exception as error:
                results["patches"][name] = {"error": str(error), "exception": type(error).__name__}
            print(name, "completed", flush=True)
    results["elapsed_seconds"] = time.perf_counter() - started
    results["patches"] = {name: results["patches"][name] for name in arguments.patches}
    output.mkdir(parents=True, exist_ok=True)
    (output / "results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")
    (output / "SUMMARY.md").write_text(summary(results))


if __name__ == "__main__":
    main()
