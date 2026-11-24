"""Compare retained TASK067 profiles at common accepted model numbers."""

import argparse
import importlib.util
import json
import math
from pathlib import Path
import sys


CODE = Path(__file__).resolve().parent
WORKSPACE = CODE.parents[4]
EXPERIMENT = WORKSPACE / "freeeos_biquintic_progress/experiments/067_stellar_burning"
ORIGINAL = CODE.parent / "stellar_experiment/compare_profiles.py"
sys.dont_write_bytecode = True
specification = importlib.util.spec_from_file_location("frozen_profile_reader", ORIGINAL)
reader = importlib.util.module_from_spec(specification)
specification.loader.exec_module(reader)
FIELDS = reader.FIELDS


def index(run_dir):
    entries = {}
    for line in (run_dir / "LOGS/profiles.index").read_text().splitlines():
        values = line.split()
        if len(values) == 3 and all(value.isdigit() for value in values):
            entries[int(values[0])] = int(values[2])
    return entries


def compare(pair):
    baseline = EXPERIMENT / f"baseline_{pair}"
    candidate = EXPERIMENT / f"candidate_{pair}"
    baseline_index = index(baseline)
    candidate_index = index(candidate)
    common = sorted(set(baseline_index) & set(candidate_index))
    maxima = {field: {"absolute_difference": 0.0, "model": None, "zone": None}
              for field in FIELDS}
    result = {
        "pair": pair, "baseline_retained": len(baseline_index),
        "candidate_retained": len(candidate_index),
        "common_model_numbers": common, "profiles_compared": 0,
        "mass_grid_mismatch_models": [], "max_abs_age_difference_yr": 0.0,
        "max_abs_mass_difference_msun": 0.0,
        "max_abs_profile_differences": maxima,
    }
    for model_number in common:
        left_path = baseline / "LOGS" / f"profile{baseline_index[model_number]}.data"
        right_path = candidate / "LOGS" / f"profile{candidate_index[model_number]}.data"
        left_stream, left_header, left_columns = reader.profile_lines(left_path)
        right_stream, right_header, right_columns = reader.profile_lines(right_path)
        try:
            if int(left_header["model_number"]) != model_number or int(right_header["model_number"]) != model_number:
                raise ValueError(f"profile index/header mismatch at model {model_number}")
            if not set(FIELDS + ("mass",)).issubset(left_columns) or not set(FIELDS + ("mass",)).issubset(right_columns):
                raise ValueError(f"required profile column absent at model {model_number}")
            result["max_abs_age_difference_yr"] = max(result["max_abs_age_difference_yr"],
                abs(float(left_header["star_age"]) - float(right_header["star_age"])))
            left_rows = [list(map(float, line.split())) for line in left_stream]
            right_rows = [list(map(float, line.split())) for line in right_stream]
            if len(left_rows) != len(right_rows):
                result["mass_grid_mismatch_models"].append(model_number)
                continue
            model_mass_difference = max((abs(left[left_columns["mass"]] - right[right_columns["mass"]])
                                         for left, right in zip(left_rows, right_rows)), default=0.0)
            result["max_abs_mass_difference_msun"] = max(result["max_abs_mass_difference_msun"], model_mass_difference)
            if model_mass_difference > 1e-11:
                result["mass_grid_mismatch_models"].append(model_number)
                continue
            for zone, (left, right) in enumerate(zip(left_rows, right_rows), start=1):
                if not all(math.isfinite(value) for value in left + right):
                    raise ValueError(f"nonfinite profile value at model {model_number}, zone {zone}")
                for field in FIELDS:
                    difference = abs(left[left_columns[field]] - right[right_columns[field]])
                    if difference > maxima[field]["absolute_difference"]:
                        maxima[field] = {"absolute_difference": difference,
                                         "model": model_number, "zone": zone}
            result["profiles_compared"] += 1
            if model_number == common[-1]:
                result["final_headers"] = {backend: {key: float(header[key]) for key in
                    ("star_age", "Teff", "photosphere_L", "photosphere_r", "center_h1", "center_he3")}
                    for backend, header in (("baseline", left_header), ("candidate", right_header))}
        finally:
            left_stream.close()
            right_stream.close()
    result["max_abs_profile_differences"] = maxima
    output = EXPERIMENT / f"profile_comparison_{pair}.json"
    if output.exists():
        raise FileExistsError(output)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: result[key] for key in
          ("baseline_retained", "candidate_retained", "profiles_compared",
           "max_abs_age_difference_yr", "max_abs_mass_difference_msun",
           "mass_grid_mismatch_models")}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair", required=True)
    compare(parser.parse_args().pair)


if __name__ == "__main__":
    main()
