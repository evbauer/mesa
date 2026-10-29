"""End-to-end strict standalone cache/physical-coordinate Fortran checks."""

import argparse
import json
import os
import subprocess
from pathlib import Path

import numpy as np

from freeeos_potential import LN10, checksum, construct_twists, evaluate, query_cells, read_table, source_constraints, thermodynamic_partials, to_physical, valid_cells
from run_audit import distribution, own_use


ORDERS = ((0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2), (3, 0), (2, 1), (1, 2), (0, 3))


def main():
    workspace = Path(__file__).resolve().parents[4]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--driver", type=Path, default=Path(__file__).with_name("table_fortran") / "build/table_query")
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=workspace / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data")
    parser.add_argument("--baseline", type=Path, default=workspace / "freeeos_report")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    metadata = json.loads(arguments.cache.with_suffix(arguments.cache.suffix + ".json").read_text())
    if metadata["source_sha256"] != checksum(arguments.source):
        raise ValueError("Exported source checksum mismatch")
    header, logq, logt, rows = read_table(arguments.source)
    low, valid, _ = source_constraints(rows, logq, logt, split_radiation=metadata["split_radiation"])
    nodes, complete, _ = construct_twists(low, valid, (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10,
                                        metadata["radius"], metadata["route"])
    samples = np.genfromtxt(arguments.baseline / "sample_residuals.csv", delimiter=",", names=True)
    density_index, temperature_index, inside = query_cells(logq, logt, samples["logQ"], samples["logT"])
    selected = inside & valid_cells(complete)[density_index, temperature_index] & own_use(samples["logQ"], samples["logT"])
    sample_q, sample_t = samples["logQ"][selected], samples["logT"][selected]
    logdensity, logtemperature = (sample_q + 2 * sample_t - 12) * LN10, sample_t * LN10
    input_text = "\n".join(f"{density:.17e} {temperature:.17e}" for density, temperature in zip(logdensity, logtemperature)) + "\n"
    environment = os.environ.copy()
    sdk = workspace / "mesasdk-x86_64-linux-26.6.1"
    environment["LD_LIBRARY_PATH"] = f"{sdk}/lib:{sdk}/lib64:" + environment.get("LD_LIBRARY_PATH", "")
    command = [str(arguments.driver.resolve()), str(arguments.cache.resolve()), str(header[1]), str(header[2]), metadata["source_sha256"]]
    process = subprocess.run(command, input=input_text, text=True, capture_output=True, check=True, env=environment)
    lines = process.stdout.splitlines()
    if len(lines) != len(sample_q) + 1:
        raise ValueError("Unexpected table-driver output count")
    fields = lines[0].split()
    numeric = np.stack([np.fromstring(line, sep=" ") for line in lines[1:]])
    if numeric.shape != (len(sample_q), len(fields)):
        raise ValueError("Unexpected table-driver output shape")
    returned = {field: numeric[:, index] for index, field in enumerate(fields)}
    for field in ("status_real", "status_ad", "status_state"):
        if np.any(returned[field] != 0):
            raise ValueError(f"Unexpected valid-query failure in {field}")
    native = evaluate(nodes, complete, logq, logt, sample_q, sample_t, add_radiation=metadata["split_radiation"])
    physical = to_physical(native, ORDERS)
    state = thermodynamic_partials(native, sample_q, sample_t)
    result = {"status": "strict standalone adapter only; no MESA runtime/default change", "count": len(sample_q),
              "provenance": {"command": command, "source_sha256": checksum(arguments.source), "cache_sha256": checksum(arguments.cache),
                             "driver_sha256": checksum(arguments.driver), "script_sha256": checksum(__file__), "cache_metadata": metadata},
              "potential_physical_partials": {}, "thermodynamics": {}}
    potential_names = ("f", "f_x", "f_t", "f_xx", "f_xt", "f_tt", "f_xxx", "f_xxt", "f_xtt", "f_ttt")
    for name, orders in zip(potential_names, ORDERS):
        scale = np.maximum.reduce([abs(physical[orders]), abs(returned[name]), abs(returned["ad_" + name]), np.ones(len(sample_q))])
        result["potential_physical_partials"][name] = {"real_python": distribution((returned[name] - physical[orders]) / scale),
                                                       "ad_python": distribution((returned["ad_" + name] - physical[orders]) / scale),
                                                       "real_ad": distribution((returned[name] - returned["ad_" + name]) / scale)}
    for field in ("Pgas", "E", "S", "Cv", "Cp", "chiRho", "chiT", "gamma1", "gamma3", "grad_ad"):
        for suffix in ("", "_x", "_t"):
            if field + suffix not in returned:
                raise ValueError(f"Driver missing required state field {field + suffix}")
            scale = np.maximum.reduce([abs(state[field]), abs(state[field + suffix]), np.ones(len(sample_q))])
            result["thermodynamics"][field + suffix] = distribution((returned[field + suffix] - state[field + suffix]) / scale)
    result["metadata_rejections"] = {}
    for label, expected_x, expected_z, expected_hash in (("composition_X", .6, header[2], metadata["source_sha256"]),
                                                        ("composition_Z", header[1], .04, metadata["source_sha256"]),
                                                        ("source_hash", header[1], header[2], "0" * 64)):
        rejection = subprocess.run([str(arguments.driver.resolve()), str(arguments.cache.resolve()), str(expected_x), str(expected_z), expected_hash],
                                   input="", text=True, capture_output=True, env=environment)
        if rejection.returncode == 0:
            raise ValueError(f"Metadata mismatch not rejected: {label}")
        result["metadata_rejections"][label] = {"exit_code": rejection.returncode, "stderr": rejection.stderr.strip()}
    (output / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"count": result["count"], "Cv_p99": result["thermodynamics"]["Cv"]["p99"],
                      "metadata_rejections": list(result["metadata_rejections"])}, indent=2))


if __name__ == "__main__":
    main()
