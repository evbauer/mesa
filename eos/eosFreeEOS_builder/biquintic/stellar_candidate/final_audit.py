"""Read-only first-pair audit and explicit, unstarted 100-step proposal."""

import json
import math
from itertools import zip_longest

from analyze import records
from harness import CODE, CONTROL, RUNS, SUMMARY, digest, inventory, placement, write_json


def byte_count(path):
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file()) if path.is_dir() else path.stat().st_size


def audit_run(name, candidate):
    directory = RUNS / name
    manifest = json.loads((directory / "manifest.json").read_text())
    for filename, expected in manifest["inputs"].items():
        if digest(directory / filename) != expected:
            raise ValueError("input identity mismatch: " + filename)
    backend_lines = [line.split() for line in (directory / "backend.txt").read_text().splitlines()]
    reference = [float(value) for line in backend_lines if line[0] == "pin_xa" for value in line[1:]]
    if len(reference) != 20:
        raise ValueError("wrong reference composition length")
    query_count = 0
    positive_count = 0
    component_count = 0
    relative_error = 0.0
    for sequence, row in enumerate(records(directory / "queries.tsv"), 1):
        query_count = sequence
        if int(row[0]) != sequence:
            raise ValueError("nonsequential final query record")
        fractions = [float(value) for value in row[6:14]]
        if any(not math.isfinite(value) or not 0 <= value <= 1 for value in fractions):
            raise ValueError("invalid fraction")
        if abs(sum(fractions) - 1) > 1e-12 or fractions[-1] != 0:
            raise ValueError("fraction normalization/WDM violation")
        if float(row[14]) > 5e-13:
            raise ValueError("final relative composition violation")
        if fractions[2] > 0:
            positive_count += 1
    positive_queries = (row for row in records(directory / "queries.tsv") if float(row[8]) > 0) if candidate else ()
    for sequence, (component, final) in enumerate(zip_longest(records(directory / "component_queries.tsv"), positive_queries), 1):
        if component is None or final is None:
            raise ValueError("component/final count mismatch")
        component_count = sequence
        if int(component[0]) != sequence or component[1:4] != final[1:4]:
            raise ValueError("component/final phase mismatch")
        if list(map(float, component[4:6])) != list(map(float, final[4:6])):
            raise ValueError("component/final state mismatch")
        if int(component[10]) != 0 or component[11] != "T":
            raise ValueError("candidate status/hybrid failure")
        abundances = list(map(float, component[12:]))
        if len(abundances) != 20 or any(not math.isfinite(value) or value <= 0 for value in abundances):
            raise ValueError("invalid component abundance")
        difference = max(abs(value - pinned) / pinned for value, pinned in zip(abundances, reference))
        relative_error = max(relative_error, difference)
        if difference > 5e-13 or abs(sum(abundances) - 1) > 5e-13:
            raise ValueError("component relative composition violation")
    event_error = max(float(row[7]) for row in records(directory / "events.tsv"))
    if not math.isfinite(event_error) or event_error > 5e-13:
        raise ValueError("lifecycle composition violation")
    raw_names = ("queries.tsv", "component_queries.tsv", "events.tsv", "iterations.tsv", "matrices.bin", "LOGS", "photos")
    raw_bytes = {name: byte_count(directory / name) for name in raw_names}
    return {"successful_final_queries": query_count, "component_calls": component_count,
            "positive_FreeEOS_queries": positive_count, "inactive_FreeEOS_queries": query_count - positive_count,
            "component_final_sequence_state_phase_match": True,
            "max_component_relative_composition_error": relative_error if candidate else None,
            "max_lifecycle_relative_composition_error": event_error,
            "raw_logging_bytes": raw_bytes, "total_raw_logging_bytes": sum(raw_bytes.values()),
            "private_cache_bytes": byte_count(directory / "runtime_cache"),
            "final_model_sha256": digest(directory / "final.mod")}


def main():
    execution = placement()
    comparison = json.loads((SUMMARY / "comparison.json").read_text())
    required = ("comparison_completed", "identical_inputs", "identical_executable", "identical_accepted_timesteps",
                "private_optout_matches049_profiles_exactly", "private_optout_matches049_matrices_exactly")
    if not all(comparison[key] for key in required):
        raise ValueError("first-pair comparison gate failed")
    freeze = json.loads((RUNS / "driver_build/freeze.json").read_text())
    if inventory() != freeze["dependencies"] or digest(RUNS / "driver_build/star") != freeze["executable_sha256"]:
        raise ValueError("frozen dependencies or executable changed")
    for filename, expected in freeze["sources"].items():
        if digest(CODE / filename) != expected:
            raise ValueError("frozen source changed: " + filename)
    baseline = audit_run("baseline_00", False)
    candidate = audit_run("candidate_00", True)
    profile_bytes_equal = all(digest(CONTROL / f"LOGS/profile{number}.data") ==
                              digest(RUNS / f"baseline_00/LOGS/profile{number}.data") for number in range(1, 9))
    for backend in ("baseline", "candidate"):
        result = comparison[backend]
        if (result["accepted_steps"], result["attempts"], result["deliberate_retries"], result["numerical_retries"]) != (8, 9, 1, 0):
            raise ValueError("unexpected step/retry accounting")
        if result["accepted_iterations"] != [2] * 8:
            raise ValueError("unexpected accepted iterations")
    projection_factor = 100 / 8
    extension = {"name": "TASK052_same_seed_100", "status": "PROPOSED_NOT_STARTED_REQUIRES_LEAD_CLEARANCE",
                 "run_names": ["baseline_100", "candidate_100"], "accepted_steps_each": 100,
                 "common_seed_sha256": "cd5b74e7ccd55d8253d045abc5bb95cba1506ba42fba70482bccf5f29cf9eb3e",
                 "executable_sha256": freeze["executable_sha256"],
                 "only_input_change_from_first_pair": "max_model_number=100 in both inlists",
                 "same_forced_retry": "one model2 retry in each run",
                 "threads_each": 1, "memory_bytes_each": 8 * 1024**3,
                 "timeout_seconds_each": 900, "sequential_runs": True,
                 "projected_stellar_seconds": {backend: comparison[backend]["outcome"]["seconds"] * projection_factor
                                               for backend in ("baseline", "candidate")},
                 "projected_raw_logging_bytes": {"baseline": baseline["total_raw_logging_bytes"] * projection_factor,
                                                 "candidate": candidate["total_raw_logging_bytes"] * projection_factor},
                 "private_cache_bytes_pair": baseline["private_cache_bytes"] + candidate["private_cache_bytes"],
                 "proposed_raw_logging_cap_bytes_each": 1024**3,
                 "proposed_total_incremental_disk_reserve_bytes": 12 * 1024**3,
                 "projection_limit": "linear eight-step extrapolation; extra iterations/retries may increase time and logs",
                 "age_comparison": "report overlap [max(first saved ages),min(last saved ages)]; exact-age comparisons first, otherwise bracketed per-zone temporal interpolation only within overlap, explicitly approximate; no extrapolation or model-index equivalence",
                 "failure_policy": "stop and preserve unsupported/composition/hybrid failure, timeout or logging-cap attempt; no fallback, control retune, EOS fix or automatic rerun"}
    result = {"execution": execution, "status": "FIRST_PAIR_AUDIT_PASS", "baseline": baseline, "candidate": candidate,
              "all_inputs_and_frozen_sources_verified": True, "private_installed_dependencies_unchanged": True,
              "private_optout_profiles_byte_exact_vs049": profile_bytes_equal,
              "first_pair_saved_profile_age_overlap_years": [comparison["baseline"]["accepted_age_years"][0],
                                                              comparison["baseline"]["final_age_years"]],
              "proposal": extension}
    write_json(SUMMARY / "HANDOFF.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
