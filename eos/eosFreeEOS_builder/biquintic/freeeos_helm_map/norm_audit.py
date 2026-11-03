"""Audit archived 048 FreeEOS seeds and HELM face denominator floors."""

import json

import numpy as np

from common import OUTPUT, WORKSPACE, digest, save_json, scaled_jump


OLD048 = WORKSPACE / "freeeos_biquintic_progress/experiments/048_helm_control/results.json"


def traces(item, left_key, right_key):
    left, right = item[left_key], item[right_key]
    response_left = left["Cv"] if isinstance(left, dict) else left[0]
    response_right = right["Cv"] if isinstance(right, dict) else right[0]
    partial_left = left["Cv_t"] if isinstance(left, dict) else left[2]
    partial_right = right["Cv_t"] if isinstance(right, dict) else right[2]
    jump, denominator = scaled_jump(response_left, response_right, partial_left, partial_right)
    return float(jump), float(denominator)


def population(items, left_key, right_key, jump_key):
    checked = []
    for item in items:
        jump, denominator = traces(item, left_key, right_key)
        old = item[jump_key]["Cv_t"]["signed_scaled"]
        checked.append((denominator, abs(jump - old)))
    return {"count": len(checked), "minimum_denominator": min(value[0] for value in checked),
            "denominator_below_one": sum(value[0] < 1 for value in checked),
            "max_old_vs_no_floor_jump_difference": max(value[1] for value in checked)}


def main():
    frozen = json.loads(OLD048.read_text())
    seeds = [item for item in frozen["seeds"] if item["group"] == "frozen_FreeEOS_24"]
    faces = frozen["faces"]
    report = {"archived048_FreeEOS_seeds": population(seeds, "split_left", "split_right", "FreeEOS_jump"),
              "archived048_HELM_faces": population(faces, "total_left", "total_right", "total_jump"),
              "old048_sha256": digest(OLD048),
              "metric": "no denominator floor; Cv and Cv_t terms retained"}
    output = OUTPUT / "normalization_audit.json"
    if output.exists():
        raise FileExistsError("normalization audit already saved")
    save_json(output, report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
