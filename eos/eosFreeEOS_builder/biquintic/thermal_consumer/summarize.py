"""Summarize frozen-consumer trials without rerunning numerical models."""

import json
from collections import Counter, defaultdict

import numpy as np

from models import WORKSPACE
from freeeos_potential import checksum


OUTPUT = WORKSPACE / "freeeos_biquintic_progress/experiments/039_thermal_consumer"


def quantiles(values):
    values = np.asarray(values, dtype=float)
    return {"median": float(np.median(values)), "p90": float(np.quantile(values, .9)),
            "max": float(np.max(values))}


def main():
    records = [json.loads(line) for line in (OUTPUT / "trials.jsonl").read_text().splitlines()]
    checks = json.loads((OUTPUT / "jacobian_checks.json").read_text())
    baseline = json.loads((OUTPUT / "results.json").read_text())
    if baseline["manifest_sha256"] != "afc6fe676ba13bf858fc0a611268402b3663f4ad1e0ff991049713fd78788f19":
        raise ValueError("Only frozen V2 results may be summarized")
    if len(records) != 1020 or any(record["status"] != "converged" for record in records):
        raise ValueError("Unexpected trial population or failure; inspect raw records first")
    result = {"scope": "Summary of canonical V2 table-only solver records; root differences are not source errors",
              "trials_sha256": checksum(OUTPUT / "trials.jsonl"),
              "jacobian_checks_sha256": checksum(OUTPUT / "jacobian_checks.json"),
              "manifest_sha256": baseline["manifest_sha256"],
              "status_counts": dict(Counter(record["status"] for record in records)),
              "families": {}, "paired": {}, "jacobian": {}}
    for family in ("common", "off_face_control", "manufactured", "restart_guess", "tolerance_sweep"):
        result["families"][family] = {}
        for name in ("candidate", "legacy", "smooth"):
            selected = [record for record in records if record["family"] == family and record["model"] == name]
            result["families"][family][name] = {
                "count": len(selected), "iterations": quantiles([record["iterations"] for record in selected]),
                "root_error_abs": quantiles([abs(record["root_error"]) for record in selected]),
                "rejected_total": sum(record["rejected"] for record in selected),
                "bisected_total": sum(record["bisected"] for record in selected),
                "selected_face_crossings": sum(record["crossed_selected_face"] for record in selected),
                "u_crossings": sum(record["initial_to_root_u_crossings"] for record in selected),
                "tau_crossings": sum(record["initial_to_root_tau_crossings"] for record in selected),
                "sampled_positive_jacobian_count": sum(record["reference_diagnostic"]["sampled_positive_jacobian"]
                                                       for record in selected),
                "sampled_one_root_count": sum(record["reference_diagnostic"]["sampled_sign_changes_or_zeros"] == 1
                                              for record in selected)}
    for family in ("common", "off_face_control"):
        paired = defaultdict(dict)
        for record in records:
            if record["family"] == family:
                paired[(record["case_id"], record["sign"], record["scale"])][record["model"]] = record
        result["paired"][family] = {}
        for category in ("ordinary", "tail"):
            groups = [group for group in paired.values() if group["candidate"]["category"] == category]
            no_native = [group for group in groups if group["candidate"]["initial_to_root_u_crossings"] == 0
                         and group["candidate"]["initial_to_root_tau_crossings"] == 0]
            differences = [group["candidate"]["iterations"] - group["legacy"]["iterations"] for group in groups]
            result["paired"][family][category] = {
                "count": len(groups), "candidate_minus_legacy_iterations": dict(sorted(Counter(differences).items())),
                "mean_iteration_difference": float(np.mean(differences)),
                "no_native_crossing_count": len(no_native),
                "no_native_iteration_difference": dict(sorted(Counter(
                    group["candidate"]["iterations"] - group["legacy"]["iterations"]
                    for group in no_native).items())),
                "model_root_separation_abs": quantiles([
                    abs(group["candidate"]["reference_root"] - group["legacy"]["reference_root"])
                    for group in groups])}
    result["manufactured_target_error_abs_max"] = max(abs(record["target_error"]) for record in records
                                                  if record["family"] == "manufactured")
    result["jacobian"]["within_cell_scaled_error"] = quantiles([
        abs(row["scaled_error"]) for check in checks for row in check["within_cell"]])
    result["jacobian"]["one_sided_abs_error_all_steps"] = quantiles([
        abs(row["finite_difference"] - row["exact_jacobian"])
        for check in checks for row in check["exact_face"]["one_sided"]])
    result["jacobian"]["one_sided_abs_error_small_step_max"] = max(
        abs(check["exact_face"]["one_sided"][index]["finite_difference"]
            - check["exact_face"]["one_sided"][index]["exact_jacobian"])
        for check in checks for index in (1, 3))
    result["tolerance_sweep"] = {}
    for tolerance in (1e-9, 1e-11, 1e-13):
        result["tolerance_sweep"][str(tolerance)] = {}
        for name in ("candidate", "legacy", "smooth"):
            selected = [record for record in records if record["family"] == "tolerance_sweep"
                        and record["model"] == name and record["tolerance_absolute_R"] == tolerance]
            result["tolerance_sweep"][str(tolerance)][name] = {
                "count": len(selected), "iterations": quantiles([record["iterations"] for record in selected]),
                "max_absolute_residual": max(abs(record["root_residual"]) for record in selected),
                "max_root_error": max(abs(record["root_error"]) for record in selected)}
    result["smooth_temporal_error_max_by_scale"] = {str(scale): max(
        abs(record["backward_euler_temporal_error"]) for record in records
        if record["family"] == "common" and record["model"] == "smooth" and record["scale"] == scale)
        for scale in (1e-4, 1e-3, 1e-2)}
    (OUTPUT / "analysis.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    plot_paired(paired_records=records)
    print(json.dumps({"trials": len(records), "status": result["status_counts"],
                      "paired": result["paired"], "jacobian": result["jacobian"]}, indent=2))


def plot_paired(paired_records):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    manifest = json.loads((WORKSPACE / "freeeos_biquintic_progress/experiments/038_solver_manifest/queries_v2.json").read_text())
    jumps = {case["case_id"]: abs(case["split_scaled_Cv_t_jump"]) for case in manifest["cases"]}
    figure, axis = plt.subplots(figsize=(8, 4.8))
    for family, marker, color, label in (("common", "o", "tab:blue", "Cross selected face"),
                                        ("off_face_control", "^", "tab:orange", "Stay on same side of selected face")):
        grouped = defaultdict(dict)
        for record in paired_records:
            if record["family"] == family:
                grouped[(record["case_id"], record["sign"], record["scale"])][record["model"]] = record
        axis.scatter([jumps[group["candidate"]["case_id"]] for group in grouped.values()],
                     [group["candidate"]["iterations"] - group["legacy"]["iterations"]
                      for group in grouped.values()], marker=marker, color=color, alpha=.5, label=label)
    axis.set_xscale("log")
    axis.set_yticks((0, 1))
    axis.set_xlabel("|Scaled jump in dCv/dlnT| at selected face")
    axis.set_ylabel("Biquintic minus legacy Newton iterations")
    axis.grid(alpha=.2)
    axis.legend()
    figure.tight_layout()
    figure.savefig(OUTPUT / "paired_iteration_difference.pdf")
    plt.close(figure)


if __name__ == "__main__":
    main()
