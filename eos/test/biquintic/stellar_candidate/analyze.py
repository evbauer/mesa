"""Audit the bounded stellar A/B, including profiles and native cell transitions."""

from collections import Counter
import json
import math
from pathlib import Path
import shlex
import struct

from harness import CACHE, CONTROL, RUNS, SUMMARY, WORKSPACE, digest, placement, write_json


def records(path):
    if not path.exists():
        return
    with path.open() as stream:
        for line in stream:
            if line.strip() and not line.startswith("#"):
                yield line.split()


def mesa_profile(path):
    lines = path.read_text().splitlines()
    header = dict(zip(lines[1].split(), shlex.split(lines[2])))
    names = lines[5].split()
    rows = [list(map(float, line.split())) for line in lines[6:] if line.strip()]
    if not all(len(row) == len(names) and all(map(math.isfinite, row)) for row in rows):
        raise ValueError("malformed/nonfinite profile: " + str(path))
    return header, {name: [row[index] for row in rows] for index, name in enumerate(names)}


def matrix_inventory(path):
    count, scalar_count = 0, 0
    dimensions = set()
    with path.open("rb") as stream:
        while header := stream.read(16):
            if len(header) != 16:
                raise ValueError("truncated matrix header")
            model, attempt, hydro, zones = struct.unpack("=4i", header)
            if min(model, attempt, hydro, zones) < 1 or hydro != 4 or zones != 366:
                raise ValueError("unexpected structural matrix dimensions")
            values = hydro * zones + 3 * hydro * hydro * zones
            payload = stream.read(8 * values)
            if len(payload) != 8 * values or not all(math.isfinite(item[0]) for item in struct.iter_unpack("=d", payload)):
                raise ValueError("truncated/nonfinite matrix")
            count += 1
            scalar_count += values
            dimensions.add((hydro, zones))
    return {"records": count, "finite_scalars": scalar_count, "dimensions": sorted(dimensions)}


def cell_coordinates(logrho, logtemp, grid):
    umin, du, tmin, dtemp = grid
    coordinates = (((logrho - 2 * logtemp + 12) * math.log(10) - umin) / du,
                   (logtemp * math.log(10) - tmin) / dtemp)
    ambiguous = any(abs(value - round(value)) < 1e-10 for value in coordinates)
    return tuple(math.floor(value) for value in coordinates), ambiguous


def crossings(profiles):
    with CACHE.open() as stream:
        header = [next(stream) for unused in range(6)]
    grid = list(map(float, header[5].split()))
    counts, examples = Counter(), []
    for model_number in range(2, len(profiles) + 1):
        earlier, later = profiles[model_number - 2][1], profiles[model_number - 1][1]
        for zone_index, (old_rho, old_temp, new_rho, new_temp) in enumerate(zip(
                earlier["logRho"], earlier["logT"], later["logRho"], later["logT"])):
            old_cell, old_ambiguous = cell_coordinates(old_rho, old_temp, grid)
            new_cell, new_ambiguous = cell_coordinates(new_rho, new_temp, grid)
            if old_ambiguous or new_ambiguous:
                counts["near_face_endpoint_pairs"] += 1
            delta_u, delta_temp = abs(new_cell[0] - old_cell[0]), abs(new_cell[1] - old_cell[1])
            if delta_u or delta_temp:
                counts["changed_zone_steps"] += 1
                counts["minimum_u_faces"] += delta_u
                counts["minimum_temperature_faces"] += delta_temp
                pure = min(earlier["eos_frac_FreeEOS"][zone_index], later["eos_frac_FreeEOS"][zone_index]) >= 1 - 1e-12
                counts["pure_endpoint_crossings" if pure else "blended_endpoint_crossings"] += 1
                examples.append({"model_from": model_number - 1, "model_to": model_number, "zone": zone_index + 1,
                                 "old_cell": old_cell, "new_cell": new_cell, "pure_endpoints": pure,
                                 "old_logRho_logT": [old_rho, old_temp], "new_logRho_logT": [new_rho, new_temp]})
    return {"counts": dict(counts), "events": examples,
            "scope": f"accepted profiles1through{len(profiles)} only; observed endpoint cell changes give minimum crossings on a continuous path, not complete trial trajectories or first-step crossings"}


def summarize(name, expected_steps=8):
    directory = RUNS / name
    command = json.loads((directory / "command.json").read_text())
    outcome = json.loads((directory / "outcome.json").read_text())
    completion = dict(line.split() for line in (directory / "completion.txt").read_text().splitlines()) if outcome["completion_written"] else {}
    events = list(records(directory / "events.tsv"))
    accepted = [row for row in events if row[0] == "accepted"]
    attempts = [row for row in events if row[0] == "attempt"]
    forced = sum(row[0] == "forced_retry" for row in events)
    completed = outcome["returncode"] == 0 and len(accepted) == expected_steps and completion.get("model_number") == str(expected_steps)
    completed = completed and "termination code: max_model_number" in (directory / "run.log").read_text()
    routes, phases = Counter(), Counter()
    max_relative_error, min_fraction = 0.0, 1.0
    query_ranges = {"logRho": [math.inf, -math.inf], "logT": [math.inf, -math.inf]}
    faces = json.loads((WORKSPACE / "freeeos_biquintic_progress/experiments/038_solver_manifest/queries_v2.json").read_text())["cases"]
    nearest = {face["case_id"]: {"distance_native_cells": math.inf} for face in faces}
    for row in records(directory / "queries.tsv"):
        if len(row) != 15:
            raise ValueError("bad final-query schema")
        values = list(map(float, row[4:]))
        if not all(map(math.isfinite, values)):
            raise ValueError("nonfinite final query")
        logrho, logtemp = values[:2]
        fractions = values[2:10]
        freeeos = fractions[2]
        if fractions[7] != 0 or abs(sum(fractions) - 1) > 1e-12 or values[-1] > 5e-13:
            raise ValueError("route or relative-composition guard failed")
        min_fraction = min(min_fraction, freeeos)
        max_relative_error = max(max_relative_error, values[-1])
        pure = freeeos >= 1 - 1e-12 and max(fractions[:2] + fractions[3:]) <= 1e-12
        routes["pure_FreeEOS" if pure else "blended_FreeEOS" if freeeos > 0 else "other_EOS"] += 1
        phases[row[3]] += 1
        for key, value in (("logRho", logrho), ("logT", logtemp)):
            query_ranges[key] = [min(query_ranges[key][0], value), max(query_ranges[key][1], value)]
        if pure:
            logq = logrho - 2 * logtemp + 12
            for face in faces:
                distance = max(abs(logq - face["logQ"]) / .03, abs(logtemp - face["logT"]) / .02)
                if distance < nearest[face["case_id"]]["distance_native_cells"]:
                    nearest[face["case_id"]] = {"distance_native_cells": distance, "query": int(row[0]),
                                                "logRho": logrho, "logT": logtemp, "category": face["category"]}
    components, component_errors, hybrid_failures = 0, 0, 0
    for row in records(directory / "component_queries.tsv"):
        if len(row) != 32:
            raise ValueError("bad component-query schema")
        components += 1
        component_errors += int(row[10]) != 0
        hybrid_failures += row[11] != "T"
        if not all(math.isfinite(float(value)) for value in row[4:10] + row[12:]):
            raise ValueError("nonfinite component-query input")
    expected_components = routes["pure_FreeEOS"] + routes["blended_FreeEOS"] if command["backend"] == "CANDIDATE" else 0
    profiles = [mesa_profile(directory / f"LOGS/profile{number}.data") for number in range(1, len(accepted) + 1)]
    iterations = list(records(directory / "iterations.tsv"))
    total_retries = max((int(row[3]) for row in events), default=0)
    result = {"run": name, "backend": command["backend"], "scope": "exploratory frozen-composition stellar relaxation",
              "outcome": outcome, "completed": completed, "accepted_steps": len(accepted), "attempts": len(attempts),
              "total_retries": total_retries, "deliberate_retries": forced, "numerical_retries": total_retries - forced,
              "accepted_iterations": [int(row[6]) for row in accepted], "logged_iterations_including_rejected": len(iterations),
              "accepted_dt_seconds": [float(row[5]) for row in accepted],
              "accepted_age_years": [float(row[4]) for row in accepted], "final_age_years": float(accepted[-1][4]) if accepted else None,
              "routes": dict(routes), "phases": dict(phases), "min_FreeEOS_fraction": min_fraction,
              "query_ranges": query_ranges, "max_relative_composition_error": max_relative_error,
              "component_calls": components, "component_errors": component_errors, "hybrid_failures": hybrid_failures,
              "callback_count_matches_final_routes": components == expected_components,
              "completion_counters": completion, "native_crossings": crossings(profiles),
              "nearest_manifest_faces": nearest, "matrix_inventory": matrix_inventory(directory / "matrices.bin"),
              "matrices_sha256": digest(directory / "matrices.bin"), "executable_sha256": command["executable_sha256"]}
    write_json(SUMMARY / f"{name}.json", result)
    return result, profiles


def compare_profiles(baseline_profiles, candidate_profiles):
    result = []
    for model_number, (baseline, candidate) in enumerate(zip(baseline_profiles, candidate_profiles), 1):
        baseline_header, baseline_data = baseline
        candidate_header, candidate_data = candidate
        if baseline_data["zone"] != candidate_data["zone"] or baseline_data["mass"] != candidate_data["mass"]:
            raise ValueError("profiles do not share the same Lagrangian mesh")
        differences = {}
        for quantity in ("logRho", "logT", "logP", "logR", "luminosity", "cp", "cv", "grada", "eps_grav", "eos_frac_FreeEOS"):
            pairs = list(zip(baseline_data[quantity], candidate_data[quantity]))
            absolute = [abs(candidate_value - baseline_value) for baseline_value, candidate_value in pairs]
            scaled = [abs(candidate_value - baseline_value) / max(1.0, abs(candidate_value), abs(baseline_value))
                      for baseline_value, candidate_value in pairs]
            max_index = max(range(len(absolute)), key=absolute.__getitem__)
            differences[quantity] = {"maximum_absolute": max(absolute), "max_absolute_zone": max_index + 1,
                                     "maximum_scaled": max(scaled), "rms_absolute": math.sqrt(sum(value * value for value in absolute) / len(absolute))}
        result.append({"model": model_number, "baseline_age": float(baseline_header["star_age"]),
                       "candidate_age": float(candidate_header["star_age"]), "differences": differences})
    return result


def main():
    placement()
    baseline, baseline_profiles = summarize("baseline_00")
    candidate, candidate_profiles = summarize("candidate_00")
    baseline_inputs = json.loads((RUNS / "baseline_00/manifest.json").read_text())["inputs"]
    candidate_inputs = json.loads((RUNS / "candidate_00/manifest.json").read_text())["inputs"]
    same_inputs = baseline_inputs == candidate_inputs
    same_executable = baseline["executable_sha256"] == candidate["executable_sha256"]
    profiles = compare_profiles(baseline_profiles, candidate_profiles)
    old_profiles = [mesa_profile(CONTROL / f"LOGS/profile{number}.data") for number in range(1, 9)]
    old_comparison = compare_profiles(old_profiles, baseline_profiles)
    optout_exact = all(all(value["maximum_absolute"] == 0 for value in model["differences"].values()) for model in old_comparison)
    result = {"scope": "exploratory eight-step frozen-composition stellar relaxation; not production or global jump-stress acceptance",
              "baseline": baseline, "candidate": candidate, "identical_inputs": same_inputs, "identical_executable": same_executable,
              "identical_accepted_timesteps": baseline["accepted_dt_seconds"] == candidate["accepted_dt_seconds"],
              "profile_comparison": profiles, "profile_scaling": "absolute difference; scaled denominator=max(1,abs(base),abs(candidate)); log quantities are dex",
              "private_optout_matches049_profiles_exactly": optout_exact,
              "private_optout_matches049_matrices_exactly": digest(CONTROL / "matrices.bin") == baseline["matrices_sha256"],
              "comparison_completed": baseline["completed"] and candidate["completed"] and same_inputs and same_executable
                  and baseline["callback_count_matches_final_routes"] and candidate["callback_count_matches_final_routes"]
                  and candidate["component_errors"] == 0 and candidate["hybrid_failures"] == 0
                  and baseline["outcome"]["protected_unchanged"] and candidate["outcome"]["protected_unchanged"],
              "limitations": ["051tight1e-13stress failures remain retained; standard-tolerance consumer evidence motivated this exploratory run",
                              "no composition-changing evolution, no independent source-accuracy claim, no global rough-face stress",
                              "no additional extension started; any next bounded experiment requires lead clearance"]}
    write_json(SUMMARY / "comparison.json", result)
    concise = {"baseline": {key: baseline[key] for key in ("outcome", "routes", "phases", "native_crossings", "matrix_inventory")},
               "candidate": {key: candidate[key] for key in ("outcome", "routes", "phases", "native_crossings", "matrix_inventory")},
               "profile_maxima_over_eight_models": {quantity: max(model["differences"][quantity]["maximum_absolute"] for model in profiles)
                                                       for quantity in profiles[0]["differences"]} if profiles else {},
               "optout_matrix_exact": result["private_optout_matches049_matrices_exactly"],
               "nearest_candidate_tail": min((value | {"case_id": key} for key, value in candidate["nearest_manifest_faces"].items()
                                               if value.get("category") == "tail"), key=lambda value: value["distance_native_cells"])}
    write_json(SUMMARY / "report_numbers.json", concise)
    print(json.dumps({"comparison_completed": result["comparison_completed"], "identical_inputs": same_inputs,
                      "identical_executable": same_executable, "identical_accepted_timesteps": result["identical_accepted_timesteps"],
                      "optout_profiles_exact": optout_exact,
                      "baseline_iterations": baseline["accepted_iterations"], "candidate_iterations": candidate["accepted_iterations"],
                      "candidate_numerical_retries": candidate["numerical_retries"], "candidate_component_errors": candidate["component_errors"],
                      "candidate_component_calls": candidate["component_calls"], "candidate_max_relative_composition_error": candidate["max_relative_composition_error"],
                      "last_profile_differences": profiles[-1]["differences"] if profiles else {}}, indent=2))


if __name__ == "__main__":
    main()
