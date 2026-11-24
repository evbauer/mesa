"""Declared-reference heat-capacity cancellation indicators at frozen faces."""

import argparse
import json
from pathlib import Path

import numpy as np

from build_manifest import CODE, WORKSPACE, checksum
from freeeos_potential import LN10, construct_twists, evaluate, radiation_potential, read_table, source_constraints, to_physical


def indicators(physical, temperature):
    free_energy, first, second = physical[0, 0], physical[0, 1], physical[0, 2]
    heat_capacity = (first - second) / temperature
    energy = free_energy - first
    entropy = -first / temperature
    quotient_first = -energy / temperature
    quotient_second = energy / temperature - heat_capacity
    if not np.all(np.isfinite(heat_capacity)) or np.any(heat_capacity <= 0):
        raise ValueError("Require finite positive total Cv for conditioning indicators")
    return {"Cv": heat_capacity, "entropy_over_Cv": entropy / heat_capacity,
            "energy_over_TCv": energy / (temperature * heat_capacity),
            "kappa_f": (abs(first) + abs(second)) / (temperature * heat_capacity),
            "kappa_f_over_T": (abs(quotient_first) + abs(quotient_second)) / heat_capacity}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=WORKSPACE / "freeeos_biquintic_progress/experiments/044_local_conditioning")
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest_path = WORKSPACE / "freeeos_biquintic_progress/experiments/038_solver_manifest/queries_v2.json"
    if checksum(manifest_path) != "afc6fe676ba13bf858fc0a611268402b3663f4ad1e0ff991049713fd78788f19":
        raise ValueError("Frozen manifest identity changed")
    manifest = json.loads(manifest_path.read_text())
    table = WORKSPACE / manifest["provenance"]["source_table"]
    if checksum(table) != manifest["provenance"]["source_sha256"] or checksum(CODE / "freeeos_potential.py") != manifest["provenance"]["core_sha256"]:
        raise ValueError("Frozen source/core identity changed")
    header, logq, logt, rows = read_table(table)
    if not np.allclose(header[[1, 2]], [.7, .02], rtol=0, atol=1e-12):
        raise ValueError("Unexpected source composition")
    low, valid, _ = source_constraints(rows, logq, logt, split_radiation=True)
    nodes, complete, _ = construct_twists(low, valid, (logq[1] - logq[0]) * LN10,
                                         (logt[1] - logt[0]) * LN10)
    records = []
    for case in manifest["cases"]:
        sample_q, sample_t = np.array([case["logQ"]]), np.array([case["logT"]])
        temperature = 10.0**sample_t
        cell = tuple(np.array([index]) for index in case["left_cell"])
        remainder = to_physical(evaluate(nodes, complete, logq, logt, sample_q, sample_t,
                                         cell_indices=cell, add_radiation=False))
        analytic = to_physical(radiation_potential(sample_q, sample_t))
        total = remainder + analytic
        values = indicators(total, temperature)
        values["kappa_radiation_subtracted_over_total_TCv"] = (
            abs(remainder[0, 1]) + abs(remainder[0, 2])) / (temperature * values["Cv"])
        records.append({"case_id": case["case_id"], "direction": case["direction"], "category": case["category"],
                        "logQ": case["logQ"], "logT": case["logT"],
                        **{name: float(value[0]) for name, value in values.items()}})
    fields = ("entropy_over_Cv", "energy_over_TCv", "kappa_f", "kappa_f_over_T", "kappa_radiation_subtracted_over_total_TCv")
    summaries = {}
    for group in ("all", "ordinary", "tail"):
        selected = [record for record in records if group == "all" or record["category"] == group]
        summaries[group] = {name: {"min": min(record[name] for record in selected),
                                  "median": float(np.median([record[name] for record in selected])),
                                  "max": max(record[name] for record in selected)} for name in fields}
    result = {"scope": "Arithmetic cancellation indicators for the existing potential, not a new f/T fit or a source-accuracy/conditioning bound on all operations",
              "reference_limit": "Indicators depend on energy/entropy gauge; unknown offsets are preserved. Radiation subtraction/restoration does not choose or repair an unknown zero point.",
              "definitions": {"kappa_f": "(|f_t|+|f_tt|)/|T Cv|", "kappa_f_over_T": "(|q_t|+|q_tt|)/|Cv| with q=f/T, q_t=-E/T, q_tt=E/T-Cv",
                              "kappa_radiation_subtracted_over_total_TCv": "(|f_remainder,t|+|f_remainder,tt|)/|T Cv_total|; analytic radiation contribution is restored exactly"},
              "provenance": {"manifest_sha256": checksum(manifest_path), "source_sha256": checksum(table),
                             "core_sha256": checksum(CODE / "freeeos_potential.py"), "script_sha256": checksum(__file__)},
              "cases": records, "summary": summaries}
    (output / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
