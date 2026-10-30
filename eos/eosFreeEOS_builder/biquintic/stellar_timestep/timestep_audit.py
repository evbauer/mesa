"""Read-only saved-state audits and coarse-prefix gate for TASK057."""

from collections import Counter
import json
import math
import struct

import timestep as task
from analyze import cell_coordinates, compare_profiles, matrix_inventory, mesa_profile, records


def retained_profiles(directory):
    index_lines = (directory / "LOGS/profiles.index").read_text().splitlines()
    indexed = []
    for line in index_lines[1:]:
        if not line.strip():
            continue
        model, priority, profile_number = map(int, line.split())
        path = directory / f"LOGS/profile{profile_number}.data"
        profile = mesa_profile(path)
        if int(profile[0]["model_number"]) != model:
            raise ValueError("profile index/model metadata mismatch")
        indexed.append((model, profile_number, profile))
    indexed.sort(key=lambda item: item[0])
    if len({item[0] for item in indexed}) != len(indexed):
        raise ValueError("duplicate retained model")
    return indexed


def retained_crossings(indexed):
    with task.frozen.CACHE.open() as stream:
        header = [next(stream) for unused in range(6)]
    grid = list(map(float, header[5].split()))
    counts, events, gaps = Counter(), [], []
    adjacent = 0
    for earlier, later in zip(indexed, indexed[1:]):
        if later[0] != earlier[0] + 1:
            gaps.append([earlier[0], later[0]])
            continue
        adjacent += 1
        old_data, new_data = earlier[2][1], later[2][1]
        for zone_index, (old_rho, old_temp, new_rho, new_temp) in enumerate(zip(
                old_data["logRho"], old_data["logT"], new_data["logRho"], new_data["logT"])):
            old_cell, old_near = cell_coordinates(old_rho, old_temp, grid)
            new_cell, new_near = cell_coordinates(new_rho, new_temp, grid)
            if old_near or new_near:
                counts["near_face_endpoint_pairs"] += 1
            delta_u, delta_temp = abs(new_cell[0] - old_cell[0]), abs(new_cell[1] - old_cell[1])
            if delta_u or delta_temp:
                counts["changed_zone_steps"] += 1
                counts["minimum_u_faces"] += delta_u
                counts["minimum_temperature_faces"] += delta_temp
                pure = min(old_data["eos_frac_FreeEOS"][zone_index], new_data["eos_frac_FreeEOS"][zone_index]) >= 1 - 1e-12
                counts["pure_endpoint_crossings" if pure else "blended_endpoint_crossings"] += 1
                events.append({"model_from": earlier[0], "model_to": later[0], "zone": zone_index + 1,
                               "old_cell": old_cell, "new_cell": new_cell, "pure_endpoints": pure})
    return {"counts": dict(counts), "events": events, "adjacent_retained_pairs": adjacent, "excluded_gaps": gaps,
            "scope": "only consecutive retained model pairs from profiles.index; missing models, seed-first, and within-step/rejected paths excluded"}


def coarse_gate(name, accepted, final_profile):
    backend = name.split("_", 1)[1]
    old_directory = task.frozen.RUNS / (backend + "_100")
    old_accepted = [row for row in records(old_directory / "events.tsv") if row[0] == "accepted"]
    same_count = len(accepted) == len(old_accepted) == 100
    same_prefix_events = same_count and all(new[1:7] == old[1:7] for new, old in zip(accepted[:-1], old_accepted[:-1]))
    same_prefix_profiles = same_count and all(task.frozen.digest(task.RUNS / name / f"LOGS/profile{number}.data") ==
                                              task.frozen.digest(old_directory / f"LOGS/profile{number}.data") for number in range(1, 100))
    same_matrix_prefix = True
    prefix_records = 0
    with (old_directory / "matrices.bin").open("rb") as old_stream, (task.RUNS / name / "matrices.bin").open("rb") as new_stream:
        while old_header := old_stream.read(16):
            model, attempt, hydro, zones = struct.unpack("=4i", old_header)
            if model >= 100:
                break
            payload_size = 8 * (hydro * zones + 3 * hydro * hydro * zones)
            same_matrix_prefix = same_matrix_prefix and new_stream.read(16) == old_header
            same_matrix_prefix = same_matrix_prefix and new_stream.read(payload_size) == old_stream.read(payload_size)
            prefix_records += 1
            if not same_matrix_prefix:
                break
    old_final = mesa_profile(old_directory / "LOGS/profile100.data")
    differences = compare_profiles([old_final], [final_profile])[0]["differences"]
    maximum_scaled = max(value["maximum_scaled"] for value in differences.values())
    final_dt_difference = abs(float(accepted[-1][5]) - float(old_accepted[-1][5]))
    final_dt_bound = 2 * math.ulp(task.TARGET_AGE) * 31557600
    final_roundoff_gate = maximum_scaled <= 5e-11 and final_dt_difference <= final_dt_bound
    final_roundoff_gate = final_roundoff_gate and accepted[-1][6] == old_accepted[-1][6]
    result = {"passed": same_count and same_prefix_events and same_prefix_profiles and same_matrix_prefix and final_roundoff_gate,
              "accepted_count_matches052": same_count, "pre_endpoint_events_exact": same_prefix_events,
              "pre_endpoint_profile_bytes_exact": same_prefix_profiles, "pre_endpoint_matrix_bytes_exact": same_matrix_prefix,
              "pre_endpoint_matrix_records": prefix_records,
              "final_profile_bytes_exact": task.frozen.digest(task.RUNS / name / "LOGS/profile100.data") ==
                                           task.frozen.digest(old_directory / "LOGS/profile100.data") if same_count else False,
              "final_maximum_scaled_difference": maximum_scaled, "final_scaled_roundoff_ceiling": 5e-11,
              "final_dt_difference_seconds": final_dt_difference, "final_dt_roundoff_ceiling_seconds": final_dt_bound,
              "final_profile_differences": differences,
              "limit_policy": "external coarse-reproduction gate only; never changes stellar tolerances or physics"}
    return result


def query_audit(directory, candidate):
    reference = next(list(map(float, line.split()[1:])) for line in (directory / "backend.txt").read_text().splitlines()
                     if line.split() and line.split()[0] == "pin_xa")
    if len(reference) != 20:
        raise ValueError("wrong frozen composition length")
    components = iter(records(directory / "component_queries.tsv"))
    routes, phases = Counter(), Counter()
    final_count, component_count = 0, 0
    max_relative_error, component_relative_error, min_fraction = 0.0, 0.0, 1.0
    query_ranges = {"logRho": [math.inf, -math.inf], "logT": [math.inf, -math.inf]}
    for sequence, row in enumerate(records(directory / "queries.tsv"), 1):
        if len(row) != 15 or int(row[0]) != sequence:
            raise ValueError("bad final-query schema or sequence")
        values = list(map(float, row[4:]))
        if not all(map(math.isfinite, values)):
            raise ValueError("nonfinite final-query record")
        fractions = values[2:10]
        if any(not 0 <= value <= 1 for value in fractions) or abs(sum(fractions) - 1) > 1e-12 or fractions[-1] != 0:
            raise ValueError("invalid final fractions/WDM")
        if not 0 <= values[-1] <= 5e-13:
            raise ValueError("relative composition violation")
        final_count = sequence
        freeeos = fractions[2]
        min_fraction = min(min_fraction, freeeos)
        max_relative_error = max(max_relative_error, values[-1])
        pure = freeeos >= 1 - 1e-12 and max(fractions[:2] + fractions[3:]) <= 1e-12
        routes["pure_FreeEOS" if pure else "blended_FreeEOS" if freeeos > 0 else "inactive_FreeEOS"] += 1
        phases[row[3]] += 1
        for key, value in zip(("logRho", "logT"), values[:2]):
            query_ranges[key] = [min(query_ranges[key][0], value), max(query_ranges[key][1], value)]
        if candidate and freeeos > 0:
            component = next(components, None)
            component_count += 1
            if component is None or len(component) != 32 or int(component[0]) != component_count:
                raise ValueError("missing or malformed candidate component query")
            if component[1:4] != row[1:4] or list(map(float, component[4:6])) != values[:2]:
                raise ValueError("component/final state or phase mismatch")
            if int(component[10]) != 0 or component[11] != "T":
                raise ValueError("component status or hybrid failure")
            if not all(math.isfinite(float(value)) for value in component[4:10] + component[12:]):
                raise ValueError("nonfinite component input")
            abundances = list(map(float, component[12:]))
            difference = max(abs(value - pinned) / pinned for value, pinned in zip(abundances, reference))
            component_relative_error = max(component_relative_error, difference)
            if difference > 5e-13 or min(abundances) <= 0 or abs(sum(abundances) - 1) > 5e-13:
                raise ValueError("candidate composition guard violation")
    if next(components, None) is not None:
        raise ValueError("extra candidate callback or nonzero baseline callback")
    return {"final_queries": final_count, "component_calls": component_count, "component_errors": 0,
            "hybrid_failures": 0, "component_final_sequence_state_phase_match": True,
            "routes": dict(routes), "phases": dict(phases), "min_FreeEOS_fraction": min_fraction,
            "max_relative_composition_error": max_relative_error,
            "max_component_relative_composition_error": component_relative_error if candidate else None,
            "query_ranges": query_ranges}


def audit_run(name):
    directory = task.RUNS / name
    outcome = json.loads((directory / "outcome.json").read_text())
    manifest = json.loads((directory / "manifest.json").read_text())
    for filename, expected in manifest["inputs"].items():
        if task.frozen.digest(directory / filename) != expected:
            raise ValueError("input mutated during run")
    events = list(records(directory / "events.tsv"))
    accepted = [row for row in events if row[0] == "accepted"]
    attempts = [row for row in events if row[0] == "attempt"]
    if not accepted or len(accepted) > 300:
        raise ValueError("no accepted model or exceeded300-model guard")
    if [int(row[1]) for row in accepted] != list(range(1, len(accepted) + 1)):
        raise ValueError("nonsequential accepted models")
    max_lifecycle_error = max(float(row[7]) for row in events)
    if not math.isfinite(max_lifecycle_error) or max_lifecycle_error > 5e-13:
        raise ValueError("lifecycle composition violation")
    indexed = retained_profiles(directory)
    if not indexed or indexed[-1][0] != int(accepted[-1][1]):
        raise ValueError("last accepted model is not retained")
    final_profile_number = indexed[-1][1]
    final_profile = indexed[-1][2]
    profile_age = float(final_profile[0]["star_age"])
    event_age = float(accepted[-1][4])
    offset = (profile_age - task.TARGET_AGE) / math.ulp(task.TARGET_AGE)
    termination = "max_age" if "termination code: max_age\n" in (directory / "run.log").read_text() else "other"
    endpoint_passed = termination == "max_age" and 0 <= offset <= 2 and event_age == profile_age
    completion = dict(line.split() for line in (directory / "completion.txt").read_text().splitlines())
    query_result = query_audit(directory, name.endswith("candidate"))
    matrix = matrix_inventory(directory / "matrices.bin")
    if int(completion["queries"]) != query_result["final_queries"] or int(completion["component_calls"]) != query_result["component_calls"]:
        raise ValueError("completion query/callback count mismatch")
    if int(completion["component_errors"]) != 0 or int(completion["hybrid_checks"]) != query_result["component_calls"]:
        raise ValueError("completion hybrid/error count mismatch")
    if int(completion["matrix_evaluations"]) != matrix["records"] or int(completion["attempts"]) != len(attempts):
        raise ValueError("completion matrix/attempt mismatch")
    deliberate = sum(row[0] == "forced_retry" for row in events)
    total_retries = max(int(row[3]) for row in events)
    result = {"name": name, "outcome": outcome, "requested_termination": "max_age", "actual_termination": termination,
              "accepted_steps": len(accepted), "attempts": len(attempts), "deliberate_retries": deliberate,
              "numerical_retries": total_retries - deliberate, "completion_counters": completion,
              "accepted_iterations": [int(row[6]) for row in accepted], "logged_iterations": len(list(records(directory / "iterations.tsv"))),
              "accepted_dt_seconds": [float(row[5]) for row in accepted], "accepted_age_years": [float(row[4]) for row in accepted],
              "endpoint": {"target_age_years": task.TARGET_AGE, "actual_profile_age_years": profile_age,
                           "actual_event_age_years": event_age, "offset_ulps": offset, "target_ulp_years": math.ulp(task.TARGET_AGE),
                           "at_or_above_requested_age": profile_age >= task.TARGET_AGE, "passed": endpoint_passed},
              "max_lifecycle_relative_composition_error": max_lifecycle_error, "matrix_inventory": matrix,
              "query_audit": query_result, "native_crossings": retained_crossings(indexed),
              "profile_retention": {"retained_models": [item[0] for item in indexed], "retained_count": len(indexed),
                                    "missing_models": sorted(set(range(1, len(accepted) + 1)) - {item[0] for item in indexed}),
                                    "final_model": indexed[-1][0], "final_file_number": final_profile_number},
              "final_profile": str((directory / f"LOGS/profile{final_profile_number}.data").relative_to(task.frozen.WORKSPACE)),
              "final_profile_sha256": task.frozen.digest(directory / f"LOGS/profile{final_profile_number}.data"),
              "final_model_sha256": task.frozen.digest(directory / "final.mod")}
    if name.startswith("coarse_"):
        result["coarse_reproduction_gate"] = coarse_gate(name, accepted, final_profile)
    result["audit_passed"] = endpoint_passed and outcome["returncode"] == 0 and outcome["stop_reason"] is None
    result["audit_passed"] = result["audit_passed"] and outcome["protected_unchanged"] and deliberate == 1
    result["audit_passed"] = result["audit_passed"] and result.get("coarse_reproduction_gate", {"passed": True})["passed"]
    task.frozen.write_json(task.ANALYSIS / (name + ".json"), result)
    print(json.dumps({"name": name, "audit_passed": result["audit_passed"], "accepted_steps": len(accepted),
                      "endpoint": result["endpoint"], "queries": query_result["final_queries"],
                      "component_calls": query_result["component_calls"],
                      "coarse_reproduction_gate": result.get("coarse_reproduction_gate", {}).get("passed")}), flush=True)
    if not result["audit_passed"]:
        raise RuntimeError("endpoint/coarse/run audit failed; stop before further runs")
    return result
