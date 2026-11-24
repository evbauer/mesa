"""Cross-check streaming real/AD Fortran derivatives against the Python cell."""

import argparse
import json
import os
import subprocess
from pathlib import Path

import numpy as np

from freeeos_potential import (
    EVAL_ORDERS, LN10, checksum, construct_twists, evaluate, query_cells,
    read_table, source_constraints, thermodynamic_partials, valid_cells,
)
from run_audit import distribution, own_use


DRIVER_ORDERS = ((0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2), (3, 0), (2, 1), (1, 2), (0, 3))


def main():
    workspace = Path(__file__).resolve().parents[4]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--driver", type=Path, default=Path(__file__).with_name("test_fortran") / "cell_driver")
    parser.add_argument("--mesa", type=Path, default=workspace / "mesa")
    parser.add_argument("--baseline", type=Path, default=workspace / "freeeos_report")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source = arguments.mesa / "data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
    _, logq, logt, rows = read_table(source)
    low, valid, _ = source_constraints(rows, logq, logt)
    nodes, complete, _ = construct_twists(low, valid, (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10)
    samples = np.genfromtxt(arguments.baseline / "sample_residuals.csv", delimiter=",", names=True)
    density_index, temperature_index, inside = query_cells(logq, logt, samples["logQ"], samples["logT"])
    selected = inside & valid_cells(complete)[density_index, temperature_index] & own_use(samples["logQ"], samples["logT"])
    sample_q, sample_t = samples["logQ"][selected], samples["logT"][selected]
    density_index, temperature_index = density_index[selected], temperature_index[selected]
    input_lines = []
    for query, (density_cell, temperature_cell) in enumerate(zip(density_index, temperature_index)):
        endpoints = [logq[density_cell], logq[density_cell + 1], logt[temperature_cell], logt[temperature_cell + 1]]
        input_lines.append(" ".join(f"{value * LN10:.17e}" for value in endpoints))
        input_lines.append(f"{sample_q[query] * LN10:.17e} {sample_t[query] * LN10:.17e}")
        for temperature_end in range(2):
            for density_end in range(2):
                corner = nodes[:, :, density_cell + density_end, temperature_cell + temperature_end]
                input_lines.append(" ".join(f"{value:.17e}" for value in corner.ravel(order="F")))
    environment = os.environ.copy()
    sdk_lib = workspace / "mesasdk-x86_64-linux-26.6.1/lib"
    environment["LD_LIBRARY_PATH"] = str(sdk_lib) + ":" + environment.get("LD_LIBRARY_PATH", "")
    process = subprocess.run([str(arguments.driver.resolve())], input="\n".join(input_lines) + "\n",
                             text=True, capture_output=True, env=environment, check=True)
    lines = process.stdout.splitlines()
    if len(lines) != 3 * len(sample_q):
        raise ValueError("Fortran streaming driver returned unexpected record count")
    real = np.stack([np.fromstring(line, sep=" ") for line in lines[1::3]])
    automatic = np.stack([np.fromstring(line, sep=" ") for line in lines[2::3]])
    if real.shape != (len(sample_q), 10) or automatic.shape != real.shape:
        raise ValueError("Invalid Fortran derivative record shape")
    python = evaluate(nodes, complete, logq, logt, sample_q, sample_t)
    python_columns = np.stack([python[orders] for orders in DRIVER_ORDERS], axis=1)
    result = {"count": len(sample_q), "sampling": "saved off-node queries, own-use cuts, complete radius2 native stencils",
              "source_sha256": checksum(source), "driver_sha256": checksum(arguments.driver),
              "comparison": "same nodal values; Fortran HELM mirrored basis versus Python polynomial basis; independent real/AD evaluations",
              "derivatives_relative_or_absolute_floor1": {}, "derivatives_scaled_by_potential_cell_condition": {}}
    free_scale = abs(python[0, 0])
    spacing_u, spacing_tau = (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10
    for index, orders in enumerate(DRIVER_ORDERS):
        scale = np.maximum.reduce([abs(real[:, index]), abs(automatic[:, index]), abs(python_columns[:, index]), np.ones(len(sample_q))])
        result["derivatives_relative_or_absolute_floor1"][str(orders)] = {
            "real_python": distribution((real[:, index] - python_columns[:, index]) / scale),
            "ad_python": distribution((automatic[:, index] - python_columns[:, index]) / scale),
            "real_ad": distribution((real[:, index] - automatic[:, index]) / scale)}
        conditioned_scale = np.maximum(scale, free_scale / (spacing_u**orders[0] * spacing_tau**orders[1]))
        result["derivatives_scaled_by_potential_cell_condition"][str(orders)] = {
            "real_python": distribution((real[:, index] - python_columns[:, index]) / conditioned_scale),
            "ad_python": distribution((automatic[:, index] - python_columns[:, index]) / conditioned_scale),
            "real_ad": distribution((real[:, index] - automatic[:, index]) / conditioned_scale)}
    python_state = thermodynamic_partials(python, sample_q, sample_t)
    result["thermodynamic_values_and_log_jacobians_relative_or_response_scale"] = {}
    for label, columns in (("real", real), ("ad", automatic)):
        native = np.full_like(python, np.nan)
        for index, orders in enumerate(DRIVER_ORDERS):
            native[orders] = columns[:, index]
        state = thermodynamic_partials(native, sample_q, sample_t)
        comparisons = {}
        for field in ("Pgas", "E", "S", "Cv", "Cp", "chiRho", "chiT", "gamma1", "gamma3", "grad_ad"):
            for suffix in ("", "_x", "_t"):
                scale = np.maximum.reduce([abs(python_state[field]), abs(python_state[field + suffix]), np.ones(len(sample_q))])
                comparisons[field + suffix] = distribution((state[field + suffix] - python_state[field + suffix]) / scale)
        result["thermodynamic_values_and_log_jacobians_relative_or_response_scale"][label] = comparisons
    (output / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    (output / "driver_stderr.txt").write_text(process.stderr)
    print(json.dumps({"count": result["count"], "real_ad_max_f": result["derivatives_relative_or_absolute_floor1"]["(0, 0)"]["real_ad"]["max"],
                      "real_python_p99_Cv": result["thermodynamic_values_and_log_jacobians_relative_or_response_scale"]["real"]["Cv"]["p99"]}, indent=2))


if __name__ == "__main__":
    main()
