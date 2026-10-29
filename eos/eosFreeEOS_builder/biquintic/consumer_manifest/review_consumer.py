"""Independent accounting of saved thermal-consumer trials and path routing."""

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from build_manifest import WORKSPACE, checksum


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trials", type=Path, default=WORKSPACE / "freeeos_biquintic_progress/experiments/039_thermal_consumer/trials.jsonl")
    parser.add_argument("--paths", type=Path, default=WORKSPACE / "freeeos_biquintic_progress/experiments/041_consumer_path_selection/states/classified_paths.csv")
    parser.add_argument("--output", type=Path, default=WORKSPACE / "freeeos_biquintic_progress/experiments/041_consumer_path_selection/lead_review.json")
    arguments = parser.parse_args()
    records = [json.loads(line) for line in arguments.trials.read_text().splitlines() if line.strip()]
    identifiers = [record["trial_id"] for record in records]
    if len(set(identifiers)) != len(records):
        raise ValueError("Duplicate trial identifiers")
    rows_by_trial = defaultdict(list)
    with arguments.paths.open(newline="") as stream:
        for row in csv.DictReader(stream):
            rows_by_trial[row["trial_id"]].append(row)
    if set(rows_by_trial) != set(identifiers):
        raise ValueError("Probe/trial inventories disagree")
    groups = defaultdict(list)
    for record in records:
        groups[(record["family"], record["model"], record["category"], record["direction"])].append(record)
    summaries = {}
    for key, selected in groups.items():
        successful = [record for record in selected if record["status"] == "converged"]
        pure = []
        for record in selected:
            checked = [row for row in rows_by_trial[record["trial_id"]]
                       if row["stage"] in ("initial", "reference_root", "accepted_iterate")]
            if not checked:
                raise ValueError("Trial has no probed accepted/initial state")
            if all(row["classification"] == "pure_FreeEOS" for row in checked):
                pure.append(record["trial_id"])
        summaries["/".join(key)] = {
            "count": len(selected), "statuses": dict(Counter(record["status"] for record in selected)),
            "all_recorded_accepted_states_pure_FreeEOS": len(pure),
            "selected_face_crossing_count": sum(bool(record.get("crossed_selected_face")) for record in successful),
            "iterations_median": float(np.median([record["iterations"] for record in successful])) if successful else None,
            "iterations_max": max((record["iterations"] for record in successful), default=None),
            "rejections_total": sum(record.get("rejected", 0) for record in selected),
            "bisections_total": sum(record.get("bisected", 0) for record in selected),
            "max_abs_root_error_lnT": max((abs(record["root_error"]) for record in successful), default=None),
            "max_abs_root_residual": max((abs(record["root_residual"]) for record in successful), default=None),
        }
    paired = defaultdict(dict)
    for record in records:
        if record["family"] == "common":
            key = (record["case_id"], record["sign"], record["scale"])
            if record["model"] in paired[key]:
                raise ValueError("Duplicate common-forcing model comparison")
            paired[key][record["model"]] = record
    if len(paired) != 144 or any(set(pair) != {"candidate", "legacy", "smooth"} for pair in paired.values()):
        raise ValueError("Expected 144 complete common-forcing comparisons")
    differences = []
    for key, pair in paired.items():
        first, second = pair["candidate"], pair["legacy"]
        if first["initial_y"] != second["initial_y"] or first["H_erg_per_g_s"] != second["H_erg_per_g_s"]:
            raise ValueError("Common-forcing models do not use identical initial state/forcing")
        if first["status"] == second["status"] == "converged":
            differences.append({"case_id": key[0], "sign": key[1], "scale": key[2],
                "iteration_difference_candidate_minus_legacy": first["iterations"] - second["iterations"],
                "physical_model_root_difference_lnT": first["reference_root"] - second["reference_root"]})
    result = {"scope": "Lead accounting of saved solver outputs; not an independent recomputation of EOS roots or proof of continuous-path routing",
              "provenance": {"trials_sha256": checksum(arguments.trials), "classified_paths_sha256": checksum(arguments.paths),
                             "script_sha256": checksum(__file__)},
              "trial_count": len(records), "groups": summaries, "paired_common_forcing": differences,
              "paired_iteration_difference_counts": dict(Counter(record["iteration_difference_candidate_minus_legacy"] for record in differences)),
              "maximum_abs_common_model_root_difference_lnT": max(abs(record["physical_model_root_difference_lnT"]) for record in differences)}
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in ("trial_count", "paired_iteration_difference_counts", "maximum_abs_common_model_root_difference_lnT")}, indent=2))


if __name__ == "__main__":
    main()
