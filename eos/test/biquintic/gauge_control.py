"""Measure finite-grid gauge covariance without changing a source table."""

import argparse
import json
from pathlib import Path

import numpy as np

from freeeos_potential import EVAL_ORDERS, LN10, LOW_ORDERS, checksum, construct_twists, evaluate, query_cells, read_table, source_constraints, thermodynamic_partials, valid_cells
from run_audit import distribution, own_use
from validate_table import material_nodes


def gauge_partials(sample_q, sample_t, constant, temperature_coefficient, size=4):
    temperature = 10.0**np.asarray(sample_t)
    shape = np.broadcast_arrays(sample_q, sample_t)[0].shape
    result = np.zeros((size, size, *shape))
    result[0, 0] = constant + temperature_coefficient * temperature
    for order in range(1, size):
        result[0, order] = temperature_coefficient * temperature
    return result


def apply_gauge(low, sample_q, sample_t, constant, temperature_coefficient):
    shifted = low.copy()
    offset = gauge_partials(sample_q, sample_t, constant, temperature_coefficient, size=3)
    for orders in LOW_ORDERS:
        shifted[orders] += offset[orders]
    return shifted


def main():
    workspace = Path(__file__).resolve().parents[4]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=workspace / "freeeos_biquintic_progress/experiments/026_gauge_control")
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source = workspace / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
    _, logq, logt, rows = read_table(source)
    low, valid, _ = source_constraints(rows, logq, logt)
    spacing = (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10
    nodes, complete, _ = construct_twists(low, valid, *spacing)
    saved = np.genfromtxt(workspace / "freeeos_report/sample_residuals.csv", delimiter=",", names=True)
    sample_q, sample_t = saved["logQ"], saved["logT"]
    iq, it, inside = query_cells(logq, logt, sample_q, sample_t)
    selected = inside & valid_cells(complete)[iq, it] & valid_cells(material_nodes(rows, logq, logt))[iq, it] & own_use(sample_q, sample_t)
    sample_q, sample_t = sample_q[selected], sample_t[selected]
    reference_native = evaluate(nodes, complete, logq, logt, sample_q, sample_t)
    reference = thermodynamic_partials(reference_native, sample_q, sample_t)
    entropy_scale = float(np.median(reference["S"]))
    result = {"purpose": "Numerical gauge covariance of existing-grid potential interpolation, not source gauge alignment or physical accuracy",
              "gauge": "f -> f + a + b*T leaves P,Cv,chi invariant analytically; E shifts+a and S shifts-b",
              "sampling": "same saved offnode queries with clean centered stencils, material corners and own-use cuts",
              "not_changed": "Source rows/validity/header/cache unchanged. Added jets are artificial analytic controls; source gauges are not fitted or repaired.",
              "reference_entropy_scale": entropy_scale, "count": len(sample_q),
              "provenance": {"source_sha256": checksum(source), "script_sha256": checksum(__file__), "preprocessing_sha256": checksum(Path(__file__).with_name("freeeos_potential.py"))}, "controls": {}}
    archives = {"logQ": sample_q, "logT": sample_t}
    grid_q, grid_t = logq[:, None], logt[None, :]
    for name, constant, temperature_coefficient in (("constant_1e12", 1e12, 0.), ("constant_1e16", 1e16, 0.),
                                                   ("entropy_scale", 0., entropy_scale), ("entropy_scale_x100", 0., 100 * entropy_scale)):
        shifted = apply_gauge(low, grid_q, grid_t, constant, temperature_coefficient)
        gauged_nodes, gauged_complete, _ = construct_twists(shifted, valid, *spacing)
        np.testing.assert_array_equal(gauged_complete, complete)
        native = evaluate(gauged_nodes, gauged_complete, logq, logt, sample_q, sample_t)
        analytic = gauge_partials(sample_q, sample_t, constant, temperature_coefficient)
        for orders in EVAL_ORDERS:
            native[orders] -= analytic[orders]
        state = thermodynamic_partials(native, sample_q, sample_t)
        control = {"a_erg_per_g": constant, "b_erg_per_g_K": temperature_coefficient, "values": {}, "response_jacobians": {}}
        for field in ("Pgas", "E", "S", "Cv", "Cp", "chiRho", "chiT", "gamma1", "gamma3", "grad_ad"):
            residual = (state[field] - reference[field]) / np.maximum(abs(reference[field]), 1e-300)
            worst = int(np.argmax(abs(residual)))
            control["values"][field] = {**distribution(residual), "worst": {"logQ": float(sample_q[worst]), "logT": float(sample_t[worst]), "signed_difference": float(residual[worst])}}
            archives[f"{name}_{field}_relative"] = residual
        for field in ("Cv", "Cp", "chiRho", "chiT", "gamma1", "gamma3", "grad_ad"):
            for suffix in ("_x", "_t"):
                scale = np.maximum.reduce([abs(reference[field]), abs(reference[field + suffix]), abs(state[field]), abs(state[field + suffix]), np.ones(len(sample_q))])
                residual = (state[field + suffix] - reference[field + suffix]) / scale
                control["response_jacobians"][field + suffix] = distribution(residual)
                archives[f"{name}_{field + suffix}_scaled"] = residual
        result["controls"][name] = control
    result["interpretation"] = "Constant shifts are representable but can amplify floating cancellation; b*T is exponential in native tau and only approximately reproduced by the polynomial. Coherent analytic reference subtraction/restoration avoids interpolation of a chosen known gauge, but no unknown historical gauge is inferred here. E/S departures subtract the imposed analytic gauge before comparison."
    np.savez_compressed(output / "gauge_queries.npz", **archives)
    (output / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({name: {"Cv": data["values"]["Cv"], "Cv_t": data["response_jacobians"]["Cv_t"]} for name, data in result["controls"].items()}, indent=2))


if __name__ == "__main__":
    main()
