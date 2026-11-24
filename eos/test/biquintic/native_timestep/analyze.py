"""Audit TASK069 endpoints and compare mass-aligned profiles."""

import importlib.util
import json
import math
from pathlib import Path
import shlex
import sys

import numpy as np


CODE = Path(__file__).resolve().parent
WORKSPACE = CODE.parents[4]
EXPERIMENT = WORKSPACE / "freeeos_biquintic_progress/experiments/069_native_timestep"
FIELDS = ("logRho", "logT", "logP", "logR", "h1", "he4", "eos_frac_FreeEOS")
sys.dont_write_bytecode = True
specification = importlib.util.spec_from_file_location("task069_harness", CODE / "harness.py")
harness = importlib.util.module_from_spec(specification)
specification.loader.exec_module(harness)


def completion(path):
    return dict(line.split(maxsplit=1) for line in path.read_text().splitlines() if line.strip())


def table(path):
    with path.open() as stream:
        next(stream)
        names = next(stream).split()
        values = shlex.split(next(stream))
        next(stream)
        next(stream)
        columns = next(stream).split()
        data = np.loadtxt(stream, ndmin=2)
    return dict(zip(names, values)), {name: index for index, name in enumerate(columns)}, data


def final_profile(run_dir):
    entries = [list(map(int, line.split())) for line in (run_dir / "LOGS/profiles.index").read_text().splitlines()
               if len(line.split()) == 3 and all(part.isdigit() for part in line.split())]
    model, _, number = entries[-1]
    header, columns, data = table(run_dir / "LOGS" / f"profile{number}.data")
    if model != int(header["model_number"]):
        raise RuntimeError("final profile index/header mismatch")
    return {"retained": len(entries), "model": model, "number": number,
            "header": header, "columns": columns, "data": data}


def audit_run(backend, cap):
    run_dir = harness.directory(backend, cap)
    outcome = json.loads((run_dir / "outcome.json").read_text())
    state = completion(run_dir / "completion.txt")
    profile = final_profile(run_dir)
    history_header, history_columns, history = table(run_dir / "LOGS/history.data")
    events = [line.split() for line in (run_dir / "events.tsv").read_text().splitlines() if line and not line.startswith("#")]
    accepted = [row for row in events if row[0] == "accepted"]
    startup = next(row for row in events if row[0] == "startup")
    iterations = [line.split() for line in (run_dir / "iterations.tsv").read_text().splitlines() if line and not line.startswith("#")]
    ages = np.array([float(row[4]) for row in accepted])
    timesteps = np.array([float(row[5]) / 31557600 for row in accepted])
    cap_years = harness.CAPS[cap]
    total_years = float(np.sum(timesteps))
    cap_mask = timesteps >= cap_years * (1 - 1e-10)
    age = float(accepted[-1][4])
    endpoint = {
        "age_events_yr": age,
        "age_history_yr": float(history[-1, history_columns["star_age"]]),
        "age_profile_yr": float(profile["header"]["star_age"]),
        "age_target_ulp_error": abs(age - harness.MAX_AGE) / math.ulp(float(harness.MAX_AGE)),
        "model_number": int(state["model_number"]),
        "profile_model_number": profile["model"],
        "history_model_number": int(history[-1, history_columns["model_number"]]),
        "profile_number": profile["number"],
        "profiles_retained": profile["retained"],
        "zones": len(profile["data"]),
        "teff_k": float(profile["header"]["Teff"]),
        "luminosity_lsun": float(profile["header"]["photosphere_L"]),
        "radius_rsun": float(profile["header"]["photosphere_r"]),
        "center_h1": float(profile["header"]["center_h1"]),
        "center_he3": float(profile["header"]["center_he3"]),
        "center_he4": float(profile["header"]["center_he4"]),
        "center_h1_depletion": float(startup[10]) - float(accepted[-1][10]),
    }
    result = {"backend": backend, "cap": cap, "outcome": outcome,
              "completion": state, "endpoint": endpoint,
              "accepted_models": len(accepted),
              "attempt_rows": sum(row[0] == "attempt" for row in events),
              "retry_rows": sum(row[0] == "retry" for row in events),
              "iteration_rows": len(iterations),
              "maximum_newton_iteration": max((int(row[2]) for row in iterations), default=0),
              "maximum_retry_index": max((int(row[3]) for row in accepted), default=0),
              "dt_years_min": float(np.min(timesteps)),
              "dt_years_median": float(np.median(timesteps)),
              "dt_years_max": float(np.max(timesteps)),
              "cap_limited_steps": int(np.sum(cap_mask)),
              "cap_limited_fraction_steps": float(np.mean(cap_mask)),
              "cap_limited_fraction_time": float(np.sum(timesteps[cap_mask]) / total_years),
              "sum_accepted_dt_years": total_years,
              "first_ten_dt_years": timesteps[:10].tolist(),
              "last_ten_dt_years": timesteps[-10:].tolist(),
              "age_first_yr": float(ages[0]),
              "mass_grid_min_msun": float(np.min(profile["data"][:, profile["columns"]["mass"]])),
              "mass_grid_max_msun": float(np.max(profile["data"][:, profile["columns"]["mass"]])),
              "history_columns": list(history_columns),
              "profile_columns": list(profile["columns"]),
              "input_sha256": json.loads((run_dir / "manifest.json").read_text())["inputs"],
              "final_model_present": (run_dir / "final.mod").is_file(),
              "installed_unchanged": json.loads((run_dir / "installed_before.json").read_text()) == json.loads((run_dir / "installed_after.json").read_text())}
    result["endpoint_pass"] = (outcome["returncode"] == 0 and not outcome["timed_out"]
        and state["termination_code"] == "1" and endpoint["age_target_ulp_error"] <= 2
        and endpoint["age_history_yr"] == age and endpoint["age_profile_yr"] == age
        and endpoint["model_number"] == endpoint["profile_model_number"] == endpoint["history_model_number"]
        and result["final_model_present"] and endpoint["profiles_retained"] == len(accepted)
        and all(state[key] == "T" for key in ("saw_burn", "saw_mix", "saw_abundance_equations"))
        and state["native_failures"] == "0" and result["installed_unchanged"]
        and (int(state["native_successes"]) > 0 if backend == "native" else state["native_successes"] == "0"))
    return result, profile


def compare(left, right):
    left_mass = left["data"][:, left["columns"]["mass"]]
    right_mass = right["data"][:, right["columns"]["mass"]]
    left_order = np.argsort(left_mass)
    right_order = np.argsort(right_mass)
    left_mass = left_mass[left_order]
    right_mass = right_mass[right_order]
    if np.any(np.diff(left_mass) <= 0) or np.any(np.diff(right_mass) <= 0):
        raise RuntimeError("profile mass coordinates are not strictly monotone")
    low = max(left_mass[0], right_mass[0])
    high = min(left_mass[-1], right_mass[-1])
    grid = np.union1d(left_mass[(left_mass >= low) & (left_mass <= high)],
                      right_mass[(right_mass >= low) & (right_mass <= high)])
    if len(grid) < 2:
        raise RuntimeError("insufficient common mass support")
    matched_grid = len(left_mass) == len(right_mass) and np.array_equal(left_mass, right_mass)
    nearest = []
    for source, target in ((left_mass, right_mass), (right_mass, left_mass)):
        positions = np.searchsorted(target, source)
        lower = target[np.clip(positions - 1, 0, len(target) - 1)]
        upper = target[np.clip(positions, 0, len(target) - 1)]
        nearest.append(float(np.max(np.minimum(np.abs(source - lower), np.abs(source - upper)))))
    fields = {}
    for field in FIELDS:
        left_value = np.interp(grid, left_mass, left["data"][left_order, left["columns"][field]])
        right_value = np.interp(grid, right_mass, right["data"][right_order, right["columns"][field]])
        difference = right_value - left_value
        index = int(np.argmax(np.abs(difference)))
        fields[field] = {"max_abs": float(abs(difference[index])),
                         "signed_at_max": float(difference[index]),
                         "mass_at_max_msun": float(grid[index]),
                         "rms": float(np.sqrt(np.mean(difference**2)))}
    return {"left_zones": len(left_mass), "right_zones": len(right_mass),
            "exact_mass_grid_match": matched_grid,
            "max_nearest_mass_grid_distance_msun": max(nearest),
            "max_same_index_mass_difference_msun": float(np.max(np.abs(left_mass - right_mass))) if len(left_mass) == len(right_mass) else None,
            "common_mass_min_msun": float(low), "common_mass_max_msun": float(high),
            "union_grid_points": len(grid),
            "excluded_left_zones": int(np.sum((left_mass < low) | (left_mass > high))),
            "excluded_right_zones": int(np.sum((right_mass < low) | (right_mass > high))),
            "fields_right_minus_left": fields}


def contrast_change(profiles):
    masses = {}
    orders = {}
    for key, profile in profiles.items():
        values = profile["data"][:, profile["columns"]["mass"]]
        orders[key] = np.argsort(values)
        masses[key] = values[orders[key]]
    low = max(values[0] for values in masses.values())
    high = min(values[-1] for values in masses.values())
    grid = np.unique(np.concatenate([values[(values >= low) & (values <= high)] for values in masses.values()]))
    if len(grid) < 2:
        raise RuntimeError("insufficient four-way mass support")
    fields = {}
    for field in FIELDS:
        values = {key: np.interp(grid, masses[key], profile["data"][orders[key], profile["columns"][field]])
                  for key, profile in profiles.items()}
        difference = (values["native_half"] - values["legacy_half"]
                      - values["native_coarse"] + values["legacy_coarse"])
        index = int(np.argmax(np.abs(difference)))
        fields[field] = {"max_abs": float(abs(difference[index])),
                         "signed_at_max": float(difference[index]),
                         "mass_at_max_msun": float(grid[index]),
                         "rms": float(np.sqrt(np.mean(difference**2)))}
    return {"common_mass_min_msun": float(low), "common_mass_max_msun": float(high),
            "union_grid_points": len(grid),
            "excluded_zones": {key: int(np.sum((values < low) | (values > high))) for key, values in masses.items()},
            "fields": fields}


def main():
    results = {}
    profiles = {}
    for backend in harness.BACKENDS:
        for cap in harness.CAPS:
            key = f"{backend}_{cap}"
            results[key], profiles[key] = audit_run(backend, cap)
    contrasts = {}
    for cap in harness.CAPS:
        contrasts[f"eos_{cap}"] = compare(profiles[f"legacy_{cap}"], profiles[f"native_{cap}"])
    for backend in harness.BACKENDS:
        contrasts[f"cap_{backend}"] = compare(profiles[f"{backend}_coarse"], profiles[f"{backend}_half"])
    contrasts["change_in_eos_contrast"] = contrast_change(profiles)
    scalars = ("teff_k", "luminosity_lsun", "radius_rsun", "center_h1", "center_he3", "center_he4")
    scalar_contrasts = {cap: {field: results[f"native_{cap}"]["endpoint"][field] - results[f"legacy_{cap}"]["endpoint"][field]
                             for field in scalars} for cap in harness.CAPS}
    within_backend = {backend: {field: results[f"{backend}_half"]["endpoint"][field] - results[f"{backend}_coarse"]["endpoint"][field]
                                for field in scalars} for backend in harness.BACKENDS}
    change = {field: scalar_contrasts["half"][field] - scalar_contrasts["coarse"][field] for field in scalars}
    common_inputs = [{name: digest for name, digest in item["input_sha256"].items() if name != "inlist"}
                     for item in results.values()]
    inlists = {key: (harness.directory(item["backend"], item["cap"]) / "inlist").read_text()
               for key, item in results.items()}
    normalized = {key: inlists[key].replace("FreeEOS_use_biquintic_interpolation", "FreeEOS_use_HELM_interpolation").replace(f"max_years_for_timestep = {harness.CAPS[item['cap']]}", "max_years_for_timestep = @CAP@")
                        .replace(f"FreeEOS_use_HELM_interpolation = {harness.BACKENDS[item['backend']]}",
                                 "FreeEOS_use_HELM_interpolation = @SWITCH@")
                  for key, item in results.items()}
    final = {"runs": results, "profile_comparisons": contrasts,
             "endpoint_eos_native_minus_legacy": scalar_contrasts,
             "endpoint_half_minus_coarse": within_backend,
             "endpoint_change_in_eos_contrast": change,
             "common_non_inlist_input_hashes_equal": all(item == common_inputs[0] for item in common_inputs),
             "normalized_inlists_equal": len(set(normalized.values())) == 1,
             "endpoint_ages_exactly_equal": len({results[key]["endpoint"]["age_events_yr"] for key in results}) == 1}
    harness.frozen.record(EXPERIMENT / "analysis.json", final)
    print(json.dumps({"runs": {key: {"endpoint": item["endpoint"], "endpoint_pass": item["endpoint_pass"]}
                              for key, item in results.items()},
                      "endpoint_eos_native_minus_legacy": scalar_contrasts,
                      "endpoint_change_in_eos_contrast": change}, indent=2))


if __name__ == "__main__":
    main()
