"""Summarize real stellar logs and cache support without evaluating the candidate."""

import argparse
from bisect import bisect_right
from collections import Counter
import json
import math
from pathlib import Path
import struct

from harness import EXPERIMENT, WORKSPACE, digest, placement, write_json


CACHE = WORKSPACE / "freeeos_biquintic_progress/experiments/008_standalone_table/solar_split.cache"
FACES = WORKSPACE / "freeeos_biquintic_progress/experiments/038_solver_manifest/queries_v2.json"


def support_mask():
    with CACHE.open() as stream:
        magic = next(stream).strip()
        if magic != "FREEEOS_BIQUINTIC_V1":
            raise ValueError("unexpected archived cache")
        nq, nt = map(int, next(stream).split())
        next(stream)
        next(stream)
        next(stream)
        umin, du, tmin, dtemp = map(float, next(stream).split())
        valid = {}
        for record in stream:
            fields = record.split()
            valid[int(fields[0]) - 1, int(fields[1]) - 1] = int(fields[2]) == 1
    uaxis = [umin + index * du for index in range(nq)]
    taxis = [tmin + index * dtemp for index in range(nt)]

    def supported(logrho, logtemp):
        native_u = (logrho - 2 * logtemp + 12) * math.log(10)
        native_t = logtemp * math.log(10)
        if logtemp < 3.2 or not uaxis[0] <= native_u <= uaxis[-1] or not taxis[0] <= native_t <= taxis[-1]:
            return False
        cell_u = min(nq - 2, bisect_right(uaxis, native_u) - 1)
        cell_t = min(nt - 2, bisect_right(taxis, native_t) - 1)
        return all(valid[cell_u + offset_u, cell_t + offset_t] for offset_u in (0, 1) for offset_t in (0, 1))

    return supported


def records(path):
    with path.open() as stream:
        for line in stream:
            if line.strip() and not line.startswith("#"):
                yield line.split()


def matrix_inventory(path):
    count = 0
    dimensions = set()
    scalar_count = 0
    with path.open("rb") as stream:
        while header := stream.read(16):
            if len(header) != 16:
                raise ValueError("truncated matrix header")
            model, attempt, hydro, zones = struct.unpack("=4i", header)
            if min(model, attempt, hydro, zones) < 1 or hydro > 100 or zones > 100000:
                raise ValueError("invalid matrix dimensions")
            values = hydro * zones + 3 * hydro * hydro * zones
            payload = stream.read(8 * values)
            if len(payload) != 8 * values or not all(math.isfinite(item[0]) for item in struct.iter_unpack("=d", payload)):
                raise ValueError("truncated or nonfinite matrix payload")
            count += 1
            scalar_count += values
            dimensions.add((hydro, zones))
    return {"records": count, "dimensions": sorted(dimensions), "finite_scalars": scalar_count}


def summarize(name):
    placement()
    run = EXPERIMENT / name
    outcome = json.loads((run / "run_outcome.json").read_text())
    supported = support_mask()
    faces = json.loads(FACES.read_text())["cases"]
    nearest = {face["case_id"]: {"distance_native_cells": math.inf, "query": None} for face in faces}
    categories, phases = Counter(), Counter()
    maximum_error = 0.0
    ranges = {"logRho": [math.inf, -math.inf], "logT": [math.inf, -math.inf]}
    uncovered_examples = []
    for row in records(run / "queries.tsv"):
        if len(row) != 15:
            raise ValueError("query schema changed: requires all eight final EOS fractions")
        logrho, logtemp = map(float, row[4:6])
        fractions = list(map(float, row[6:14]))
        freeeos = fractions[2]
        if fractions[7] != 0 or abs(sum(fractions) - 1) > 1e-12:
            raise ValueError("invalid final EOS route")
        pure = freeeos >= 1 - 1e-12 and max(fractions[:2] + fractions[3:]) <= 1e-12
        category = "pure_FreeEOS" if pure else "blend_FreeEOS" if freeeos > 1e-12 else "other_EOS"
        categories[category] += 1
        phases[row[3]] += 1
        maximum_error = max(maximum_error, float(row[14]))
        for key, value in (("logRho", logrho), ("logT", logtemp)):
            ranges[key] = [min(ranges[key][0], value), max(ranges[key][1], value)]
        if freeeos > 0 and not supported(logrho, logtemp):
            categories["FreeEOS_queries_without_candidate_support"] += 1
            if len(uncovered_examples) < 8:
                uncovered_examples.append({"query": int(row[0]), "phase": row[3], "logRho": logrho,
                                           "logT": logtemp, "FreeEOS_fraction": freeeos})
        if pure:
            logq = logrho - 2 * logtemp + 12
            for face in faces:
                distance = max(abs(logq - face["logQ"]) / 0.03, abs(logtemp - face["logT"]) / 0.02)
                if distance < nearest[face["case_id"]]["distance_native_cells"]:
                    nearest[face["case_id"]] = {"distance_native_cells": distance, "query": int(row[0]),
                                                "logRho": logrho, "logT": logtemp, "category": face["category"]}
    events = list(records(run / "events.tsv"))
    accepted = [row for row in events if row[0] == "accepted"]
    attempts = [row for row in events if row[0] == "attempt"]
    termination = (run / "completion.txt").read_text() if (run / "completion.txt").exists() else ""
    completion_fields = dict(line.split() for line in termination.splitlines())
    completed = (outcome["returncode"] == 0 and len(accepted) == 8 and completion_fields.get("model_number") == "8"
                 and "termination code: max_model_number" in (run / "run.log").read_text())
    result = {"run": name, "scope": "installed legacy stellar control; no candidate evaluation", "outcome": outcome,
              "bounded_evolution_completed": completed, "accepted_steps": len(accepted), "attempts": len(attempts),
              "retry_count": max((int(row[3]) for row in events), default=0),
              "deliberate_retry_events": sum(row[0] == "forced_retry" for row in events),
              "accepted_iterations": [int(row[6]) for row in accepted],
              "accepted_dt_seconds": [float(row[5]) for row in accepted],
              "final_age_years": float(accepted[-1][4]) if accepted else None,
              "query_route_counts": dict(categories), "query_phase_counts": dict(phases),
              "query_ranges": ranges, "max_composition_error": maximum_error,
              "candidate_support_check": "read-only four-corner mask plus floor 3.2; no thermodynamic candidate evaluation",
              "candidate_cache_sha256": digest(CACHE), "unsupported_examples": uncovered_examples,
              "nearest_manifest_faces_in_native_cell_Linf_metric": nearest,
              "matrix_sha256": digest(run / "matrices.bin"), "matrix_inventory": matrix_inventory(run / "matrices.bin"),
              "completion": termination}
    write_json(run / "summary.json", result)
    print(json.dumps({key: result[key] for key in ("run", "bounded_evolution_completed", "accepted_steps", "retry_count",
                                                   "accepted_iterations", "query_route_counts", "max_composition_error")}, indent=2))
    return result


def compare(names):
    summaries = [summarize(name) for name in names]
    runs = [EXPERIMENT / name for name in names]
    commands = [json.loads((run / "run_command.json").read_text()) for run in runs]
    required = ["seed.mod", "inlist", "uniform_xa.data", "frozen20.net", "history_columns.list", "profile_columns.list"]
    same_inputs = all(len({digest(run / filename) for run in runs}) == 1 for filename in required)
    same_binary = len({command["executable_sha256"] for command in commands}) == 1
    same_retry = len({command["force_retry"] for command in commands}) == 1
    files = ["matrices.bin", "iterations.tsv", "events.tsv", "final.mod", "LOGS/history.data"]
    files += [f"LOGS/profile{number}.data" for number in range(1, 9)]
    byte_equal = {filename: len({digest(run / filename) for run in runs}) == 1 for filename in files}
    result = {"runs": names, "abundance_modes": [command["abundance_mode"] for command in commands],
              "all_completed": all(summary["bounded_evolution_completed"] for summary in summaries),
              "identical_input_files": same_inputs, "identical_executable": same_binary, "identical_forced_retry_policy": same_retry,
              "byte_equal": byte_equal,
              "passed": all(summary["bounded_evolution_completed"] for summary in summaries) and same_inputs and same_binary
                        and same_retry and all(byte_equal.values()) and {command["abundance_mode"] for command in commands} == {-1, 0, 1},
              "meaning": "dependency invariance for legacy frozen-composition structural residual/Jacobian; candidate remains untested"}
    write_json(EXPERIMENT / "abundance_invariance.json", result)
    print(json.dumps(result, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", nargs="+")
    parser.add_argument("--compare", action="store_true")
    arguments = parser.parse_args()
    if arguments.compare:
        compare(arguments.runs)
    else:
        for name in arguments.runs:
            summarize(name)


if __name__ == "__main__":
    main()
