"""Independently account for archived coupled trials and sampled MESA routing."""

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from build_manifest import WORKSPACE, checksum


def common_pairs(records, expected_count=288):
    pairs = defaultdict(dict)
    for record in records:
        if record["family"] != "common":
            continue
        key = (record["case_id"], record["sign"], record["scale"],
               record["density_over_temperature_increment"])
        if record["model"] in pairs[key]:
            raise ValueError("Duplicate common-forcing model")
        pairs[key][record["model"]] = record
    if len(pairs) != expected_count or any(set(pair) != {"candidate", "legacy"} for pair in pairs.values()):
        raise ValueError("Incomplete common-forcing matrix")
    for pair in pairs.values():
        for field in ("initial", "nominal_target", "pressure_target", "heating", "energy_scale", "box"):
            if pair["candidate"][field] != pair["legacy"][field]:
                raise ValueError("Common models differ in " + field)
    return pairs


def classify_rows(path, records, expected_ids):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if {row["trial_id"] for row in rows} != expected_ids:
        raise ValueError("Routing and trial inventories disagree")
    groups = defaultdict(Counter)
    nonpure_trials = set()
    for row in rows:
        record = records[row["trial_id"]]
        if row["case_id"] != record["case_id"] or row["model"] != record["model"]:
            raise ValueError("Routing row identifiers disagree")
        groups[record["family"]][row["classification"]] += 1
        if row["classification"] != "pure_FreeEOS":
            nonpure_trials.add(row["trial_id"])
    return {"row_count": len(rows), "by_family": dict(groups),
            "nonpure_trial_ids": sorted(nonpure_trials)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trials", type=Path, default=WORKSPACE / "freeeos_biquintic_progress/experiments/042_coupled_consumer/trials.jsonl")
    parser.add_argument("--routing", type=Path, default=WORKSPACE / "freeeos_biquintic_progress/experiments/043_coupled_path_selection")
    arguments = parser.parse_args()
    records = [json.loads(line) for line in arguments.trials.read_text().splitlines() if line.strip()]
    by_id = {record["trial_id"]: record for record in records}
    if len(by_id) != len(records):
        raise ValueError("Duplicate trial identifiers")
    pairs = common_pairs(records)
    paired = []
    for key, pair in pairs.items():
        candidate, legacy = pair["candidate"], pair["legacy"]
        accepted = candidate["status"] == legacy["status"] == "converged"
        paired.append({"case_id": key[0], "sign": key[1], "scale": key[2], "ratio": key[3],
                       "both_converged": accepted,
                       "iteration_difference": candidate["newton"]["iterations"] - legacy["newton"]["iterations"] if accepted else None,
                       "reference_root_difference_inf": float(np.max(np.abs(np.array(candidate["reference"]["root"]) - legacy["reference"]["root"]))) if accepted else None})
    groups = {}
    for family, model in sorted({(record["family"], record["model"]) for record in records}):
        selected = [record for record in records if (record["family"], record["model"]) == (family, model)]
        successful = [record for record in selected if record["status"] == "converged"]
        errors = [float(np.max(np.abs(np.array(record["newton"]["root"]) - record["reference"]["root"]))) for record in successful]
        for record, error in zip(successful, errors):
            if error != record["root_error_inf"]:
                raise ValueError("Saved root error disagrees with saved roots")
        groups[family + "/" + model] = {
            "count": len(selected), "statuses": dict(Counter(record["status"] for record in selected)),
            "reference_statuses": dict(Counter(record["reference"]["status"] for record in selected)),
            "iterations": dict(Counter(record["newton"]["iterations"] for record in successful)),
            "selected_face_crossings": sum(record["crossings"].get("selected_face_crossed_accepted", False) for record in selected),
            "rejected_steps": sum(len(record["newton"]["rejected_history"]) for record in selected),
            "maximum_root_error_inf": max(errors, default=None)}
    states = arguments.routing / "states/classified_paths.csv"
    spans = arguments.routing / "spans/classified_paths.csv"
    result = {"scope": "Lead accounting of saved coupled results, not a new EOS/reference calculation or continuum routing proof",
              "provenance": {"script_sha256": checksum(__file__), "trials_sha256": checksum(arguments.trials),
                             "states_sha256": checksum(states), "spans_sha256": checksum(spans)},
              "trial_count": len(records), "groups": groups, "common_pairs": paired,
              "paired_iteration_differences": dict(Counter(pair["iteration_difference"] for pair in paired)),
              "maximum_intermodel_root_difference_inf": max((pair["reference_root_difference_inf"] for pair in paired if pair["both_converged"]), default=None),
              "states": classify_rows(states, by_id, set(by_id)),
              "spans": classify_rows(spans, by_id, {record["trial_id"] for record in records if record["newton"].get("root") is not None})}
    (arguments.routing / "lead_review.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in ("trial_count", "groups", "paired_iteration_differences", "maximum_intermodel_root_difference_inf", "states", "spans")}, indent=2))


if __name__ == "__main__":
    main()
