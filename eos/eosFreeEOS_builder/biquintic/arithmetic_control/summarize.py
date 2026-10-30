"""Summarize immutable053 records and draw the report's comparison figure."""

from collections import Counter
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

DIRECTORY = Path(__file__).resolve().parent
WORKSPACE = DIRECTORY.parents[4]
OUTPUT = WORKSPACE / "freeeos_biquintic_progress/experiments/053_arithmetic_control"
BASELINE = OUTPUT.parent / "051_runtime_consumers"


def read(path):
    return json.loads(path.read_text())


def main():
    native = read(OUTPUT / "validation/native_records.json")
    data = read(OUTPUT / "validation/native_results.json")
    summary = {}
    for name in ("baseline", "variant"):
        errors = np.abs([record[name + "_scaled_error"] for record in native])
        worst_row, worst_partial = np.unravel_index(np.argmax(errors), errors.shape)
        summary[name + "_worst_native"] = dict(sample=native[worst_row], partial=int(worst_partial),
                                                scaled_error=float(errors[worst_row, worst_partial]))
    summary["native_variant_improved_entries"] = int(np.count_nonzero(
        np.abs([record["variant_scaled_error"] for record in native]) <
        np.abs([record["baseline_scaled_error"] for record in native])))
    face_changes = []
    for index in range(0, 96, 4):
        first, second = native[index], native[index + 3]
        order = 6 if first["label"].startswith("u_") else 9
        jumps = {name: second[name][order] - first[name][order] for name in ("baseline", "variant", "reference")}
        face_changes.append(dict(case_id=first["label"], order=order, **jumps,
                                 relative_change=abs(jumps["variant"] - jumps["baseline"]) / max(1, abs(jumps["reference"]))))
    summary["normal_third_jump_comparisons"] = face_changes
    summary["max_relative_jump_change"] = max(record["relative_change"] for record in face_changes)
    summary["normal_jumps_same_sign"] = all(np.sign(record["baseline"]) == np.sign(record["variant"]) for record in face_changes)
    matched = read(OUTPUT / "validation/matched_api_records.json")
    summary["matched_status_counts"] = dict(Counter(str(record["status"]) for record in matched))
    successes = [record for record in matched if record["status"] == 0]
    summary["matched_populations"] = dict(Counter(
        record["scenario"] + "/" + record["route"] + "/" +
        ("inactive" if record["values_baseline"][21] == 0 else "pure" if record["values_baseline"][21] == 1 else "blend")
        for record in successes))
    changes = np.abs([record["scaled_changes"] for record in matched])
    summary["matched_maximum_scaled_value_change"] = float(np.max(changes[:, :27]))
    summary["matched_maximum_scaled_derivative_change"] = float(np.max(changes[:, 27:81]))
    summary["matched_abundance_unchanged"] = bool(np.all(changes[:, 81:] == 0))
    component = read(OUTPUT / "consumer_variant/full/analysis.json")
    summary["variant_fd_summary"] = component["fd_summary"] if "fd_summary" in component else component.get("derivative_checks")
    summary["variant_actual_query_count"] = component["actual_query_count"]
    summary["variant_groups"] = component["groups"]
    summary["variant_scalar_by_tolerance"] = component["scalar_by_tolerance"]
    summary["variant_actual_populations"] = component["actual_query_populations"]
    summary["baseline_actual_query_count"] = read(BASELINE / "full/analysis.json")["actual_query_count"]
    (OUTPUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    samples = [record for record in read(OUTPUT / "validation/neighborhood_records.json") if record["route"] == "component"]
    figure, axes = plt.subplots(1, 2, figsize=(10, 3.6), layout="constrained")
    positions = np.arange(len(samples)) - len(samples) // 2
    for name, label in (("baseline", "Frozen051"), ("variant", "Isolated053")):
        axes[0].plot(positions, [record[name + "_residual"] for record in samples], label=label, linewidth=1.)
    axes[0].axhspan(-1e-13, 1e-13, color="gray", alpha=.2, label="Frozen tolerance")
    axes[0].set(xlabel="Binary64 neighbors of stalled051 state", ylabel="Scalar residual", title="T0980: same forcing and coordinates")
    axes[0].legend(fontsize=8)
    positions = np.arange(10)
    for name, label in (("baseline", "Frozen kernel"), ("variant", "Arithmetic variant")):
        axes[1].semilogy(positions, data[name]["rms_scaled_by_partial"], marker="o", label=label)
    axes[1].set(xticks=positions, xticklabels=["f", "u", "t", "uu", "ut", "tt", "uuu", "uut", "utt", "ttt"],
                ylabel="RMS |error| / max(1, |80-digit reference|)", title="160 native states; identical stored polynomial")
    axes[1].legend(fontsize=8)
    figure.savefig(OUTPUT / "arithmetic_comparison.pdf")
    figure.savefig(OUTPUT / "arithmetic_comparison.png", dpi=150)
    print(json.dumps({key: summary[key] for key in ("native_variant_improved_entries", "max_relative_jump_change", "normal_jumps_same_sign",
          "matched_status_counts", "matched_populations", "matched_maximum_scaled_value_change", "matched_maximum_scaled_derivative_change",
          "variant_actual_query_count")}, indent=2))


if __name__ == "__main__":
    main()
