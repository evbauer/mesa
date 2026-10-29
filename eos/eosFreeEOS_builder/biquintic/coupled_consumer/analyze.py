"""Summarize paired coupled trials without changing frozen inputs or outcomes."""

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np


WORKSPACE = Path(__file__).resolve().parents[5]
OUTPUT = WORKSPACE / "freeeos_biquintic_progress/experiments/042_coupled_consumer"


def distribution(values):
    numbers = np.abs(np.asarray(values, dtype=float))
    if not len(numbers):
        return {"count": 0}
    return {"count": len(numbers), "median": float(np.median(numbers)),
            "p90": float(np.quantile(numbers, .9)), "p99": float(np.quantile(numbers, .99)),
            "max": float(np.max(numbers))}


def csv_count(path):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        required = ["case_id", "trial_id", "model", "stage", "logRho", "logT"]
        if reader.fieldnames != required:
            raise ValueError("Unexpected path-probe schema")
        records = list(reader)
    if not all(np.isfinite(float(record[field])) for record in records for field in ("logRho", "logT")):
        raise ValueError("Nonfinite exported coordinate")
    return {"rows": len(records), "finite_rows": len(records),
            "stages": dict(Counter(record["stage"] for record in records))}


def plot_histories(paired):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 2, figsize=(8.5, 3.4), constrained_layout=True)
    for axis, case_id in zip(axes, ("u_tail_01", "tau_tail_01")):
        pair = paired[(case_id, 1, 1e-2, -1)]
        for name, color in (("candidate", "#277da1"), ("legacy", "#eb8a44")):
            history = pair[name]["newton"]["history"]
            axis.semilogy(range(len(history)),
                          [max(item["residual_inf"], 1e-16) for item in history],
                          marker="o", label=name, color=color)
        axis.set_title(case_id)
        axis.set_xlabel("accepted Newton iterate")
        axis.set_ylabel("scaled residual infinity norm")
        axis.grid(alpha=.2)
    axes[0].legend()
    figure.savefig(OUTPUT / "tail_residual_histories.pdf")
    plt.close(figure)


def main():
    trials = [json.loads(line) for line in (OUTPUT / "trials.jsonl").read_text().splitlines()]
    checks = json.loads((OUTPUT / "jacobian_checks.json").read_text())
    common = [trial for trial in trials if trial["family"] == "common"]
    paired = {}
    for trial in common:
        key = (trial["case_id"], trial["sign"], trial["scale"],
               trial["density_over_temperature_increment"])
        paired.setdefault(key, {})[trial["model"]] = trial
    if len(paired) != 288 or not all(set(pair) == {"candidate", "legacy"} for pair in paired.values()):
        raise ValueError("Unexpected matched forcing population")
    iteration_difference = Counter(pair["candidate"]["newton"]["iterations"]
                                   - pair["legacy"]["newton"]["iterations"]
                                   for pair in paired.values())
    root_differences = []
    for pair in paired.values():
        candidate = np.asarray(pair["candidate"]["reference"]["root"])
        legacy = np.asarray(pair["legacy"]["reference"]["root"])
        root_differences.append(float(np.max(np.abs(candidate - legacy))))
    jacobian = {}
    for group in ("within_cell", "one_sided_face"):
        records = [item for case in checks for item in case[group]
                   if item.get("status", "checked") == "checked"]
        jacobian[group] = {"count": len(records),
                           "scaled_error": distribution([item["scaled_error"] for item in records]),
                           "by_step": {str(step): distribution(
                               [item["scaled_error"] for item in records if item["step"] == step])
                                       for step in sorted({item["step"] for item in records})}}
    jacobian["within_cell_skipped_crossing"] = sum(
        item.get("status") == "crosses_another_cell"
        for case in checks for item in case["within_cell"])
    face_jumps = {}
    for direction in ("u", "tau"):
        for category in ("ordinary", "tail"):
            selected = [case for case in checks if case["case_id"].startswith(direction + "_" + category)]
            fields = {}
            for equation in range(2):
                for axis in range(2):
                    jumps = []
                    for case in selected:
                        values = {(item["side"], item["equation"], item["axis"]): item["analytic"]
                                  for item in case["one_sided_face"] if item["step"] == 1e-6}
                        left = values["left", equation, axis]
                        right = values["right", equation, axis]
                        jumps.append((left - right) / max(1, abs(left), abs(right)))
                    fields[f"R{equation + 1}_{'x' if axis == 0 else 't'}"] = distribution(jumps)
            face_jumps[direction + "_" + category] = fields
    statuses = Counter(trial["status"] for trial in trials)
    references = Counter(trial["reference"]["status"] for trial in trials)
    result = {
        "provenance": {"trials_sha256": hashlib.sha256((OUTPUT / "trials.jsonl").read_bytes()).hexdigest(),
                       "jacobian_checks_sha256": hashlib.sha256((OUTPUT / "jacobian_checks.json").read_bytes()).hexdigest(),
                       "analysis_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        "status_counts": dict(statuses), "reference_counts": dict(references),
        "trial_count": len(trials), "common_model_trials": len(common),
        "matched_common_forcing_pairs": len(paired),
        "paired_iteration_difference_candidate_minus_legacy": {str(key): value for key, value in sorted(iteration_difference.items())},
        "iterations_by_model": {name: dict(Counter(trial["newton"]["iterations"]
                                                   for trial in common if trial["model"] == name))
                                for name in ("candidate", "legacy")},
        "max_root_error_against_own_reference": max(trial.get("root_error_inf", 0) for trial in trials),
        "model_reference_root_difference_inf": distribution(root_differences),
        "max_manufactured_target_error_inf": max(trial.get("manufactured_target_error_inf", 0) for trial in trials),
        "total_rejections": sum(trial["newton"].get("rejected", 0) for trial in trials),
        "total_damped_steps": sum(trial["newton"].get("damped_steps", 0) for trial in trials),
        "max_jacobian_condition": max(item["condition"] for trial in trials for item in trial["newton"]["history"]),
        "sampled_monotone_reference_count": sum(trial["reference"]["diagnostic"].get("sampled_monotone", False) for trial in trials),
        "sampled_one_sign_change_count": sum(trial["reference"]["diagnostic"].get("sampled_sign_changes_or_zeros") == 1 for trial in trials),
        "selected_face_crossings_straight_common": sum(trial["crossings"].get("selected_face_crossed_straight", False) for trial in common),
        "selected_face_crossings_accepted_common": sum(trial["crossings"].get("selected_face_crossed_accepted", False) for trial in common),
        "accepted_native_crossings_common": {
            direction: sum(trial["crossings"].get("accepted_" + direction + "_crossings", 0) for trial in common)
            for direction in ("u", "tau")},
        "jacobian_checks": jacobian, "scaled_exact_face_jumps": face_jumps,
        "path_queries": csv_count(OUTPUT / "path_queries.csv"),
        "span_queries": csv_count(OUTPUT / "span_queries.csv"),
        "limitations": "Sampled monotonicity does not prove global uniqueness; native face crossings and finite span samples do not establish MESA selection or stellar behavior. Model root differences are not solver errors or source-accuracy measures.",
    }
    (OUTPUT / "analysis.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    plot_histories(paired)
    print(json.dumps({key: result[key] for key in ("status_counts", "matched_common_forcing_pairs",
                                                 "paired_iteration_difference_candidate_minus_legacy",
                                                 "max_root_error_against_own_reference")}, indent=2))


if __name__ == "__main__":
    main()
