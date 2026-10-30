"""Compare saved matched endpoints; no EOS query or stellar evolution."""

import json
import math
import time

import timestep as task
from analyze import compare_profiles, mesa_profile


QUANTITIES = ("logRho", "logT", "logP", "logR", "luminosity", "cp", "cv", "grada", "eps_grav", "eos_frac_FreeEOS")


def statistics(delta, scales):
    absolute = list(map(abs, delta))
    maximum_index = max(range(len(delta)), key=absolute.__getitem__)
    return {"minimum_signed": min(delta), "maximum_signed": max(delta), "maximum_absolute": max(absolute),
            "maximum_absolute_zone": maximum_index + 1, "rms_absolute": math.sqrt(sum(value * value for value in delta) / len(delta)),
            "maximum_scaled": max(abs(value) / scale for value, scale in zip(delta, scales)),
            "rms_scaled": math.sqrt(sum((value / scale)**2 for value, scale in zip(delta, scales)) / len(delta))}


def main():
    task.frozen.placement()
    if time.time() >= task.NUMERICAL_STOP:
        raise RuntimeError("TASK057numerical analysis cutoff20:30 reached")
    results = {name: json.loads((task.ANALYSIS / (name + ".json")).read_text()) for name in task.NAMES}
    continuation = task.continuation_authorized()
    if not continuation or not all(results[name]["audit_passed"] for name in ("coarse_baseline", "half_baseline", "half_candidate")):
        raise RuntimeError("specific exploratory continuation incomplete; original coarse gate remains failed")
    verified = task.verify_frozen(full=True)
    manifests = {name: json.loads((task.RUNS / name / "manifest.json").read_text()) for name in task.NAMES}
    commands = {name: json.loads((task.RUNS / name / "command.json").read_text()) for name in task.NAMES}
    pair_inputs = {resolution: manifests[resolution + "_baseline"]["inputs"] == manifests[resolution + "_candidate"]["inputs"]
                   for resolution in ("coarse", "half")}
    same_executable = len({command["executable_sha256"] for command in commands.values()}) == 1
    if not all(pair_inputs.values()) or not same_executable:
        raise ValueError("unmatched inputs or executable")
    for name in task.NAMES:
        if (task.RUNS / name / "inlist").read_text() != task.expected_inlist(name)[0]:
            raise ValueError("nonwhitelisted input change")
        for filename, expected in manifests[name]["inputs"].items():
            if task.frozen.digest(task.RUNS / name / filename) != expected:
                raise ValueError("post-run input identity mismatch")
    profiles = {name: mesa_profile(task.frozen.WORKSPACE / result["final_profile"])
                for name, result in results.items()}
    reference_header, reference = profiles["coarse_baseline"]
    for header, data in profiles.values():
        if data["zone"] != reference["zone"] or data["mass"] != reference["mass"]:
            raise ValueError("endpoint mass grids differ")
    ages = {name: float(header["star_age"]) for name, (header, unused) in profiles.items()}
    exact_ages = len(set(ages.values())) == 1
    comparisons, vectors, sensitivity = {}, {}, {}
    for quantity in QUANTITIES:
        coarse_base = profiles["coarse_baseline"][1][quantity]
        coarse_candidate = profiles["coarse_candidate"][1][quantity]
        half_base = profiles["half_baseline"][1][quantity]
        half_candidate = profiles["half_candidate"][1][quantity]
        delta = {
            "baseline_refinement": [half - coarse for half, coarse in zip(half_base, coarse_base)],
            "candidate_refinement": [half - coarse for half, coarse in zip(half_candidate, coarse_candidate)],
            "EOS_coarse": [candidate - baseline for candidate, baseline in zip(coarse_candidate, coarse_base)],
            "EOS_half": [candidate - baseline for candidate, baseline in zip(half_candidate, half_base)]}
        delta["EOS_change_under_refinement"] = [half - coarse for half, coarse in zip(delta["EOS_half"], delta["EOS_coarse"])]
        scales = [max(1.0, *(abs(value) for value in row)) for row in zip(coarse_base, coarse_candidate, half_base, half_candidate)]
        comparisons[quantity] = {key: statistics(values, scales) for key, values in delta.items()}
        eos_norm = comparisons[quantity]["EOS_half"]["maximum_absolute"]
        comparisons[quantity]["refinement_to_refined_EOS_max_norm_ratio"] = {
            backend: comparisons[quantity][backend + "_refinement"]["maximum_absolute"] / eos_norm if eos_norm else None
            for backend in ("baseline", "candidate")}
        floor_by_backend = {backend: results["coarse_" + backend]["coarse_reproduction_gate"]["final_profile_differences"][quantity]["maximum_absolute"]
                            for backend in ("baseline", "candidate")}
        measured_scale = sum(floor_by_backend.values())
        sensitivity[quantity] = {"coarse_clamp_maximum_absolute": floor_by_backend,
                                 "summed_empirical_endpoint_scale": measured_scale,
                                 "contrast_to_measured_scale": {key: values["maximum_absolute"] / measured_scale if measured_scale else None
                                     for key, values in comparisons[quantity].items() if key != "refinement_to_refined_EOS_max_norm_ratio"},
                                 "within_ten_times_measured_scale": [key for key, values in comparisons[quantity].items()
                                     if key != "refinement_to_refined_EOS_max_norm_ratio" and values["maximum_absolute"] <= 10 * measured_scale]}
        vectors[quantity] = delta
    first_step = {}
    for resolution in ("coarse", "half"):
        base = mesa_profile(task.RUNS / (resolution + "_baseline") / "LOGS/profile1.data")
        candidate = mesa_profile(task.RUNS / (resolution + "_candidate") / "LOGS/profile1.data")
        first_step[resolution] = compare_profiles([base], [candidate])[0]
        first_step[resolution]["warning"] = "different first-step ages across resolutions; not a matched-age temporal comparison"
    paired = {resolution: {"identical_accepted_timesteps": results[resolution + "_baseline"]["accepted_dt_seconds"] ==
                          results[resolution + "_candidate"]["accepted_dt_seconds"],
                          "identical_accepted_iterations": results[resolution + "_baseline"]["accepted_iterations"] ==
                          results[resolution + "_candidate"]["accepted_iterations"],
                          "identical_accepted_ages": results[resolution + "_baseline"]["accepted_age_years"] ==
                          results[resolution + "_candidate"]["accepted_age_years"]} for resolution in ("coarse", "half")}
    result = {"status": "TASK057_EXPLORATORY_CONTINUATION_COMPLETE_WITH_RETAINED_ORIGINAL_GATE_FAILURE", "target_age_years": task.TARGET_AGE,
              "scope": "two-resolution fixed-composition temporal sensitivity, not convergence order/source accuracy/production acceptance",
              "same_executable": same_executable, "identical_backend_inputs_at_each_resolution": pair_inputs,
              "verified052frozen_artifacts_after_runs": verified, "final_ages": ages, "exactly_identical_endpoint_ages": exact_ages,
              "endpoint_age_spread_years": max(ages.values()) - min(ages.values()), "endpoint_interpolation": False,
              "same_mass_grid": True, "paired_trajectories": paired, "runs": results, "endpoint_comparisons": comparisons,
              "original_gate_passed": False, "continuation_authorization": json.loads((task.EXPERIMENT / "CONTINUATION.json").read_text()),
              "decoder_authorization": json.loads((task.EXPERIMENT / "DECODER_CONTINUATION.json").read_text()),
              "endpoint_sensitivity": sensitivity,
              "sensitivity_scope": "measured coarse-clamp response only, not a rigorous half-run arithmetic bound; contrasts within10times this empirical scale are conservatively flagged unresolved",
              "scaling": "common per-zone denominator=max(1,absolute values of all four profiles); log quantities remain dex",
              "first_step_EOS_contrast": first_step,
              "limitations": ["two resolutions do not establish asymptotic order", "all composition/hybrid/rough-face scope limits from052remain",
                              "first-step EOS transients are sampled at different physical ages across resolutions",
                              "no time interpolation is used for endpoint comparison", "no further runs or retuning are implied"]}
    task.frozen.write_json(task.ANALYSIS / "comparison.json", result)
    task.frozen.write_json(task.ANALYSIS / "endpoint_vectors.json", {"zone": reference["zone"], "mass": reference["mass"],
                           "coarse_baseline_logT": reference["logT"], "coarse_baseline_cv": reference["cv"], "differences": vectors})
    concise = {"status": result["status"], "exact_endpoint_ages": exact_ages, "paired_trajectories": paired,
               "runs": {name: {key: result[key] for key in ("accepted_steps", "attempts", "deliberate_retries", "numerical_retries", "endpoint", "outcome")}
                        for name, result in results.items()},
               "original_gate_passed": False, "endpoint_sensitivity": sensitivity,
               "endpoint_maxima": {quantity: {key: values["maximum_absolute"] for key, values in data.items()
                                               if key != "refinement_to_refined_EOS_max_norm_ratio"}
                                   for quantity, data in comparisons.items()}}
    task.frozen.write_json(task.ANALYSIS / "report_numbers.json", concise)
    print(json.dumps(concise, indent=2))


if __name__ == "__main__":
    main()
