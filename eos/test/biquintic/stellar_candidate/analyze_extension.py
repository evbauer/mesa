"""Audit the authorized100-step pair without further stellar/EOS evaluation."""

from bisect import bisect_left
import json

from analyze import compare_profiles, summarize
from final_audit import audit_run
from harness import CODE, CONTROL, RUNS, SUMMARY, digest, inventory, placement, write_json


def common_age_profiles(baseline, candidate):
    baseline_ages = [float(header["star_age"]) for header, unused in baseline]
    candidate_ages = [float(header["star_age"]) for header, unused in candidate]
    if any(new <= old for old, new in zip(baseline_ages, baseline_ages[1:])) or any(
            new <= old for old, new in zip(candidate_ages, candidate_ages[1:])):
        raise ValueError("profile ages are not increasing")
    overlap = [max(baseline_ages[0], candidate_ages[0]), min(baseline_ages[-1], candidate_ages[-1])]
    comparisons = []
    for model_number, (baseline_profile, age) in enumerate(zip(baseline, baseline_ages), 1):
        if not overlap[0] <= age <= overlap[1]:
            continue
        upper = bisect_left(candidate_ages, age)
        if candidate_ages[upper] == age:
            interpolated = candidate[upper]
            method = {"kind": "exact_age", "candidate_model": upper + 1}
        else:
            lower = upper - 1
            span = candidate_ages[upper] - candidate_ages[lower]
            fraction = (age - candidate_ages[lower]) / span
            low_header, low_data = candidate[lower]
            unused_header, high_data = candidate[upper]
            if low_data["zone"] != high_data["zone"] or low_data["mass"] != high_data["mass"]:
                raise ValueError("temporal interpolation requires identical mass grids")
            interpolated_data = {name: values if name in ("zone", "mass") else
                                 [low + fraction * (high - low) for low, high in zip(values, high_data[name])]
                                 for name, values in low_data.items()}
            interpolated = (dict(low_header, star_age=str(age)), interpolated_data)
            method = {"kind": "approximate_linear_temporal_interpolation", "candidate_bracket_models": [lower + 1, upper + 1],
                      "bracket_years": [candidate_ages[lower], candidate_ages[upper]], "bracket_dt_years": span,
                      "interpolation_fraction": fraction}
        comparison = compare_profiles([baseline_profile], [interpolated])[0]
        comparison["model"] = model_number
        comparison["time_alignment"] = method
        comparisons.append(comparison)
    return {"overlap_years": overlap, "overlap_duration_years": max(0.0, overlap[1] - overlap[0]),
            "comparison_target": "baseline saved ages inside joint saved-profile interval",
            "no_extrapolation": True, "comparisons": comparisons,
            "exact_age_pairs": sum(row["time_alignment"]["kind"] == "exact_age" for row in comparisons),
            "interpolated_age_pairs": sum(row["time_alignment"]["kind"] != "exact_age" for row in comparisons)}


def main():
    execution = placement()
    baseline, baseline_profiles = summarize("baseline_100", 100)
    candidate, candidate_profiles = summarize("candidate_100", 100)
    baseline_manifest = json.loads((RUNS / "baseline_100/manifest.json").read_text())
    candidate_manifest = json.loads((RUNS / "candidate_100/manifest.json").read_text())
    same_inputs = baseline_manifest["inputs"] == candidate_manifest["inputs"]
    if not same_inputs:
        raise ValueError("extension input mismatch")
    for name in ("baseline_100", "candidate_100"):
        original = (CONTROL / "inlist").read_text()
        if (RUNS / name / "inlist").read_text() != original.replace("   max_model_number = 8\n", "   max_model_number = 100\n"):
            raise ValueError("unapproved control change")
    freeze = json.loads((RUNS / "driver_build/freeze.json").read_text())
    if inventory() != freeze["dependencies"] or digest(RUNS / "driver_build/star") != freeze["executable_sha256"]:
        raise ValueError("frozen dependencies/executable changed")
    for filename, expected in freeze["sources"].items():
        if digest(CODE / filename) != expected:
            raise ValueError("frozen source changed: " + filename)
    audits = {"baseline": audit_run("baseline_100", False), "candidate": audit_run("candidate_100", True)}
    same_executable = baseline["executable_sha256"] == candidate["executable_sha256"] == freeze["executable_sha256"]
    aligned = common_age_profiles(baseline_profiles, candidate_profiles)
    prefix = {}
    for backend in ("baseline", "candidate"):
        short, longer = RUNS / (backend + "_00"), RUNS / (backend + "_100")
        with (longer / "matrices.bin").open("rb") as stream:
            matrix_prefix_equal = stream.read((short / "matrices.bin").stat().st_size) == (short / "matrices.bin").read_bytes()
        prefix[backend] = {"first_eight_profile_bytes_identical": all(digest(short / f"LOGS/profile{number}.data") ==
                                                                    digest(longer / f"LOGS/profile{number}.data") for number in range(1, 9)),
                           "initial_matrix_stream_bytes_identical": matrix_prefix_equal}
    completed = same_inputs and same_executable and all(result["completed"] and result["callback_count_matches_final_routes"]
                   and result["component_errors"] == 0 and result["hybrid_failures"] == 0 and result["outcome"]["protected_unchanged"]
                   and result["outcome"].get("stop_reason") is None for result in (baseline, candidate))
    result = {"execution": execution, "extension": "TASK052_same_seed_100", "status": "COMPLETED" if completed else "INCOMPLETE_OR_FAILED",
              "scope": "exploratory frozen-composition stellar relaxation; not production or global tail-stress acceptance",
              "baseline": baseline, "candidate": candidate, "audits": audits, "identical_inputs": same_inputs,
              "identical_executable": same_executable, "only_input_change": "max_model_number:8to100 for both backends",
              "identical_accepted_timesteps": baseline["accepted_dt_seconds"] == candidate["accepted_dt_seconds"],
              "identical_accepted_iterations": baseline["accepted_iterations"] == candidate["accepted_iterations"],
              "maximum_equal_index_age_difference_years": max(abs(old - new) for old, new in zip(
                  baseline["accepted_age_years"], candidate["accepted_age_years"])),
              "first_pair_prefix_reproduction": prefix, "common_age_comparison": aligned,
              "profile_scaling": "absolute; scaled denominator=max(1,abs(base),abs(candidate)); log quantities in dex",
              "profile_maximum_absolute_over_common_ages": {quantity: max(row["differences"][quantity]["maximum_absolute"]
                  for row in aligned["comparisons"]) for quantity in aligned["comparisons"][0]["differences"]} if aligned["comparisons"] else {},
              "comparison_completed": completed, "additional_extension_started": False}
    write_json(SUMMARY / "extension_comparison.json", result)
    concise = {"extension": result["extension"], "comparison_completed": completed,
               "identical_inputs": same_inputs, "identical_executable": same_executable,
               "identical_accepted_timesteps": result["identical_accepted_timesteps"],
               "identical_accepted_iterations": result["identical_accepted_iterations"],
               "common_age_overlap_years": aligned["overlap_years"], "exact_age_pairs": aligned["exact_age_pairs"],
               "interpolated_age_pairs": aligned["interpolated_age_pairs"], "first_pair_prefix_reproduction": prefix,
               "baseline": {key: baseline[key] for key in ("accepted_steps", "attempts", "deliberate_retries", "numerical_retries",
                                                         "component_calls", "component_errors", "routes", "native_crossings", "outcome")},
               "candidate": {key: candidate[key] for key in ("accepted_steps", "attempts", "deliberate_retries", "numerical_retries",
                                                           "component_calls", "component_errors", "routes", "native_crossings", "outcome")},
               "final_profile": aligned["comparisons"][-1] if aligned["comparisons"] else None}
    write_json(SUMMARY / "extension_report_numbers.json", concise)
    print(json.dumps(concise, indent=2))
    if not completed:
        raise RuntimeError("extension did not complete its bounded comparison")


if __name__ == "__main__":
    main()
