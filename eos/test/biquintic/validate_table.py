"""Table-only fidelity, real-table faces, and response-Jacobian validation."""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

from freeeos_potential import (
    CRAD, EVAL_ORDERS, LN10, LOW_ORDERS, checksum, construct_twists, evaluate,
    identity_residuals, query_cells, read_table, source_constraints,
    thermodynamic_partials, thermodynamics, valid_cells,
)
from run_audit import distribution, own_use


THERMO_FIELDS = ("Pgas", "E", "S", "Cv", "Cp", "chiRho", "chiT", "gamma1", "gamma3", "grad_ad")


def material_nodes(rows, logq, logt):
    ratio = CRAD * (10**logt[None, :])**4 / (3 * 10**rows[:, :, 1])
    return (ratio < 1e-3) & (logq[:, None] >= -8) & (logq[:, None] <= 5) & (logt[None, :] <= 7.8)


def table_state(rows):
    return {"Pgas": 10**rows[..., 1], "E": 10**rows[..., 2], "S": 10**rows[..., 3],
            "chiRho": rows[..., 4], "chiT": rows[..., 5], "Cp": rows[..., 6],
            "Cv": rows[..., 7], "gamma1": rows[..., 13], "gamma3": rows[..., 14],
            "grad_ad": rows[..., 15]}


def relative_summary(candidate, reference, sample_q, sample_t):
    result = {}
    for field in THERMO_FIELDS:
        denominator = np.maximum(abs(reference[field]), 1e-300)
        residual = (candidate[field] - reference[field]) / denominator
        index = int(np.argmax(abs(residual)))
        result[field] = {**distribution(residual), "worst": {
            "logQ": float(sample_q[index]), "logT": float(sample_t[index]),
            "candidate": float(candidate[field][index]), "reference": float(reference[field][index]),
            "signed_relative_difference": float(residual[index])}}
    return result


def positivity(state, sample_q, sample_t):
    result = {}
    for field in ("Cv", "P_rho", "Pgas", "E", "S"):
        failed = (~np.isfinite(state[field])) | (state[field] <= 0)
        result[field] = {"count": int(failed.sum()), "examples": [
            {"logQ": float(sample_q[index]), "logT": float(sample_t[index]), "value": float(state[field][index])}
            for index in np.flatnonzero(failed)[:12]]}
    return result


def nodal_reproduction(nodes, complete, logq, logt, rows, sample_count=12000):
    accepted = valid_cells(complete) & valid_cells(material_nodes(rows, logq, logt))
    density_index, temperature_index = np.nonzero(accepted)
    generator = np.random.default_rng(20261029)
    choice = generator.choice(len(density_index), min(sample_count, len(density_index)), replace=False)
    density_index, temperature_index = density_index[choice], temperature_index[choice]
    sample_q, sample_t = logq[density_index], logt[temperature_index]
    native = evaluate(nodes, complete, logq, logt, sample_q, sample_t)
    constraint_error = {}
    for orders in LOW_ORDERS:
        reference = nodes[orders][density_index, temperature_index]
        scale = np.maximum(abs(reference), 1)
        constraint_error[str(orders)] = distribution((native[orders] - reference) / scale)
    state = thermodynamics(native, sample_q, sample_t)
    return {"count": len(sample_q), "constraints_relative_or_absolute_floor1": constraint_error,
            "returned_fields_vs_printed_nodes": relative_summary(state, table_state(rows[density_index, temperature_index]), sample_q, sample_t),
            "identities": {name: distribution(value) for name, value in identity_residuals(state).items()}}


def faces(nodes, complete, logq, logt, rows, material=True, max_faces=12000, split_radiation=False):
    cells = valid_cells(complete)
    if material:
        cells &= valid_cells(material_nodes(rows, logq, logt))
    output = {}
    generator = np.random.default_rng(20261029)
    for axis, direction in ((0, "u"), (1, "tau")):
        adjacent = cells[:-1, :] & cells[1:, :] if axis == 0 else cells[:, :-1] & cells[:, 1:]
        density_index, temperature_index = np.nonzero(adjacent)
        choice = generator.choice(len(density_index), min(max_faces, len(density_index)), replace=False)
        density_index, temperature_index = density_index[choice], temperature_index[choice]
        fraction = generator.uniform(.1, .9, len(choice))
        if axis == 0:
            sample_q = logq[density_index + 1]
            sample_t = logt[temperature_index] + fraction * (logt[1] - logt[0])
            right_indices = density_index + 1, temperature_index
        else:
            sample_q = logq[density_index] + fraction * (logq[1] - logq[0])
            sample_t = logt[temperature_index + 1]
            right_indices = density_index, temperature_index + 1
        selected = own_use(sample_q, sample_t)
        sample_q, sample_t = sample_q[selected], sample_t[selected]
        left_indices = density_index[selected], temperature_index[selected]
        right_indices = right_indices[0][selected], right_indices[1][selected]
        left = evaluate(nodes, complete, logq, logt, sample_q, sample_t, cell_indices=left_indices, add_radiation=split_radiation)
        right = evaluate(nodes, complete, logq, logt, sample_q, sample_t, cell_indices=right_indices, add_radiation=split_radiation)
        left_state = thermodynamic_partials(left, sample_q, sample_t)
        right_state = thermodynamic_partials(right, sample_q, sample_t)
        result = {"count": len(sample_q), "sampling": "two exact one-sided cell polynomials at shared face; random tangential fraction .1-.9; own-use cuts"}
        result["potential_partials_relative_or_floor1"] = {str(orders): distribution((left[orders] - right[orders]) /
            np.maximum.reduce([abs(left[orders]), abs(right[orders]), np.ones(len(sample_q))])) for orders in EVAL_ORDERS}
        result["thermodynamic_value_jump_relative"] = {field: distribution((left_state[field] - right_state[field]) /
            np.maximum(abs(left_state[field]), abs(right_state[field]))) for field in THERMO_FIELDS}
        result["thermodynamic_value_jump_absolute_cgs"] = {field: distribution(left_state[field] - right_state[field]) for field in THERMO_FIELDS}
        result["response_log_coordinate_jacobian_jump_scaled_by_response_or_partial"] = {}
        for field in ("Cv", "Cp", "chiRho", "chiT", "gamma1", "gamma3", "grad_ad"):
            for suffix in ("_x", "_t"):
                first, second = left_state[field + suffix], right_state[field + suffix]
                scale = np.maximum.reduce([abs(first), abs(second), abs(left_state[field]), abs(right_state[field]), np.ones(len(sample_q))])
                residual = (first - second) / scale
                worst = int(np.argmax(abs(residual)))
                result["response_log_coordinate_jacobian_jump_scaled_by_response_or_partial"][field + suffix] = {
                    **distribution(residual), "worst": {"logQ": float(sample_q[worst]), "logT": float(sample_t[worst]),
                                                        "left": float(first[worst]), "right": float(second[worst]), "scale": float(scale[worst])}}
        output[direction] = result
    return output


def jacobian_check(nodes, complete, logq, logt, sample_q, sample_t):
    sample_q, sample_t = sample_q[:1500], sample_t[:1500]
    state = thermodynamic_partials(evaluate(nodes, complete, logq, logt, sample_q, sample_t), sample_q, sample_t)
    result = {"count": len(sample_q), "derivatives": "physical x=lnrho,t=lnT; scaled by max(1,|field|,|analytic derivative|)", "steps": {}}
    for step in (1e-4, 3e-5, 1e-5):
        shift = step / LN10
        errors = {}
        for suffix, delta_q, delta_t in (("_x", shift, 0), ("_t", -2 * shift, shift)):
            plus = thermodynamics(evaluate(nodes, complete, logq, logt, sample_q + delta_q, sample_t + delta_t), sample_q + delta_q, sample_t + delta_t)
            minus = thermodynamics(evaluate(nodes, complete, logq, logt, sample_q - delta_q, sample_t - delta_t), sample_q - delta_q, sample_t - delta_t)
            for field in THERMO_FIELDS:
                numerical = (plus[field] - minus[field]) / (2 * step)
                analytic = state[field + suffix]
                scale = np.maximum.reduce([np.ones(len(sample_q)), abs(state[field]), abs(analytic)])
                errors[field + suffix] = distribution((numerical - analytic) / scale)
        result["steps"][str(step)] = errors
    return result


def heldout(logq, logt, rows, baseline, stride, material=True, split_radiation=False):
    coarse_q, coarse_t, coarse_rows = logq[::stride], logt[::stride], rows[::stride, ::stride]
    low, valid, _ = source_constraints(coarse_rows, coarse_q, coarse_t, split_radiation=split_radiation)
    du, dtau = (coarse_q[1] - coarse_q[0]) * LN10, (coarse_t[1] - coarse_t[0]) * LN10
    schemes = {f"radius{radius}_{route}": construct_twists(low, valid, du, dtau, radius, route)
               for radius, route in ((1, "mean"), (2, "mean"), (3, "mean"), (2, "first"), (2, "second"))}
    _, full_valid, _ = source_constraints(rows, logq, logt)
    density_index, temperature_index = np.nonzero(full_valid)
    sample_q, sample_t = logq[density_index], logt[temperature_index]
    coarse_density, coarse_temperature, inside = query_cells(coarse_q, coarse_t, sample_q, sample_t)
    selected = inside & own_use(sample_q, sample_t)
    selected &= (density_index % stride != 0) & (temperature_index % stride != 0)
    for _, complete, _ in schemes.values():
        selected &= valid_cells(complete)[coarse_density, coarse_temperature]
    if material:
        selected &= material_nodes(rows, logq, logt)[density_index, temperature_index]
        selected &= valid_cells(material_nodes(coarse_rows, coarse_q, coarse_t))[coarse_density, coarse_temperature]
    density_index, temperature_index = density_index[selected], temperature_index[selected]
    sample_q, sample_t = sample_q[selected], sample_t[selected]
    reference = table_state(rows[density_index, temperature_index])
    output = {"stride": stride, "count": len(sample_q), "material_domain": material, "split_radiation": split_radiation,
              "sampling": "existing native nodes withheld in both coordinates, common complete coarse-stencil domain, own-use cuts",
              "limitation": "coarsened interpolation fidelity versus printed nodes, NOT independent full-resolution off-node source accuracy", "schemes": {}}
    for name, (nodes, complete, _) in schemes.items():
        state = thermodynamics(evaluate(nodes, complete, coarse_q, coarse_t, sample_q, sample_t, add_radiation=split_radiation), sample_q, sample_t)
        output["schemes"][name] = {"relative_difference_from_printed_nodes": relative_summary(state, reference, sample_q, sample_t),
                                    "positivity": positivity(state, sample_q, sample_t)}
    legacy_values = np.stack([coarse_rows[:, :, 1] * LN10, coarse_rows[:, :, 2] * LN10, coarse_rows[:, :, 3] * LN10,
                             coarse_rows[:, :, 7], coarse_rows[:, :, 6], coarse_rows[:, :, 4], coarse_rows[:, :, 5],
                             coarse_rows[:, :, 13], coarse_rows[:, :, 14], coarse_rows[:, :, 15]], axis=-1)
    coefficients = baseline.spline_coefficients(legacy_values, coarse_q, coarse_t).transpose(0, 3, 1, 2)
    evaluated = baseline.evaluate_spline(coefficients, coarse_q, coarse_t, sample_q, sample_t)
    legacy_state = {field: np.exp(evaluated[:, index]) if index < 3 else evaluated[:, index]
                    for index, field in enumerate(THERMO_FIELDS)}
    output["synthetic_independent_column_bicubic"] = {
        "method": "global not-a-knot bicubic reconstructed on SAME coarsened printed grid; not installed full-resolution cache; fallback influence possible",
        "relative_difference_from_printed_nodes": relative_summary(legacy_state, reference, sample_q, sample_t)}
    return output


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
    header, logq, logt, rows = read_table(source)
    low, valid, _ = source_constraints(rows, logq, logt)
    nodes, complete, _ = construct_twists(low, valid, (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10)
    split_low, split_valid, _ = source_constraints(rows, logq, logt, split_radiation=True)
    split_nodes, split_complete, _ = construct_twists(split_low, split_valid, (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10)
    samples = np.genfromtxt(arguments.baseline / "sample_residuals.csv", delimiter=",", names=True)
    density_index, temperature_index, inside = query_cells(logq, logt, samples["logQ"], samples["logT"])
    selected = inside & valid_cells(complete)[density_index, temperature_index] & own_use(samples["logQ"], samples["logT"])
    results = {"status": "table-only standalone validation; runtime unchanged",
               "provenance": {"mesa_revision": subprocess.check_output(["git", "-C", str(arguments.mesa), "rev-parse", "HEAD"], text=True).strip(),
                              "source_sha256": checksum(source), "script_sha256": checksum(__file__),
                              "preprocessing_sha256": checksum(Path(__file__).with_name("freeeos_potential.py")),
                              "composition": {"X": float(header[1]), "Z": float(header[2])},
                              "radiation_constant": float(CRAD), "command": sys.argv},
               "nodal_reproduction": nodal_reproduction(nodes, complete, logq, logt, rows),
               "main_domain": "native complete stencils; q[-8,5], logT[3.2,7.8], source-corner Prad/Pgas<1e-3; own-use cuts for faces/heldouts",
               "faces": faces(nodes, complete, logq, logt, rows),
               "broad_own_use_stress_faces": faces(nodes, complete, logq, logt, rows, material=False),
               "split_radiation_broad_own_use_stress_faces": faces(split_nodes, split_complete, logq, logt, rows, material=False, split_radiation=True),
               "jacobian_finite_difference": jacobian_check(nodes, complete, logq, logt, samples["logQ"][selected], samples["logT"][selected]),
               "heldout": {str(stride): heldout(logq, logt, rows, baseline, stride) for stride in (2, 3, 4)},
               "broad_own_use_stress_heldout": {str(stride): heldout(logq, logt, rows, baseline, stride, material=False) for stride in (2, 3, 4)},
               "split_radiation_broad_own_use_stress_heldout": {str(stride): heldout(logq, logt, rows, baseline, stride, material=False, split_radiation=True) for stride in (2, 3, 4)}}
    (output / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps({"nodal_count": results["nodal_reproduction"]["count"], "faces": {name: value["count"] for name, value in results["faces"].items()},
                      "heldout_counts": {name: value["count"] for name, value in results["heldout"].items()}}, indent=2))


if __name__ == "__main__":
    main()
