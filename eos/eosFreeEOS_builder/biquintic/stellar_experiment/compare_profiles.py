"""Compare accepted MESA profiles at matched ages and mass coordinates."""

import argparse
import json
import math
from pathlib import Path
import shlex
from itertools import zip_longest


WORKSPACE = Path(__file__).resolve().parents[5]
EXPERIMENT = WORKSPACE / "freeeos_biquintic_progress/experiments/064_stellar_native"
FIELDS = ("logRho", "logT", "logP", "logR", "h1", "he4", "eos_frac_FreeEOS")


def profile_lines(path):
    stream = path.open()
    next(stream)
    header_names = next(stream).split()
    header_values = shlex.split(next(stream))
    header = dict(zip(header_names, header_values))
    next(stream)
    next(stream)
    columns = next(stream).split()
    return stream, header, {name: index for index, name in enumerate(columns)}


def compare(pair):
    baseline = EXPERIMENT / f"baseline_{pair}"
    candidate = EXPERIMENT / f"candidate_{pair}"
    expected = json.loads((baseline / "manifest.json").read_text()).get("max_model_number", 8)
    maxima = {name: {"absolute_difference": 0.0, "model": None, "zone": None} for name in FIELDS}
    max_mass_difference = 0.0
    max_age_difference = 0.0
    zone_counts = []
    final_headers = {}
    for model_number in range(1, expected + 1):
        left_path = baseline / "LOGS" / f"profile{model_number}.data"
        right_path = candidate / "LOGS" / f"profile{model_number}.data"
        left_stream, left_header, left_columns = profile_lines(left_path)
        right_stream, right_header, right_columns = profile_lines(right_path)
        try:
            if int(left_header["model_number"]) != model_number or int(right_header["model_number"]) != model_number:
                raise ValueError(f"model header mismatch at {model_number}")
            age_difference = abs(float(left_header["star_age"]) - float(right_header["star_age"]))
            max_age_difference = max(max_age_difference, age_difference)
            if not set(FIELDS + ("mass",)).issubset(left_columns) or not set(FIELDS + ("mass",)).issubset(right_columns):
                raise ValueError(f"missing comparison column at {model_number}")
            zones = 0
            for zones, (left_line, right_line) in enumerate(zip_longest(left_stream, right_stream), start=1):
                if left_line is None or right_line is None:
                    raise ValueError(f"zone count differs at {model_number}")
                left_values = list(map(float, left_line.split()))
                right_values = list(map(float, right_line.split()))
                if not all(map(math.isfinite, left_values + right_values)):
                    raise ValueError(f"nonfinite profile value at {model_number}:{zones}")
                mass_difference = abs(left_values[left_columns["mass"]] - right_values[right_columns["mass"]])
                max_mass_difference = max(max_mass_difference, mass_difference)
                for field in FIELDS:
                    difference = abs(left_values[left_columns[field]] - right_values[right_columns[field]])
                    if difference > maxima[field]["absolute_difference"]:
                        maxima[field] = {"absolute_difference": difference, "model": model_number, "zone": zones}
            zone_counts.append(zones)
            if model_number == expected:
                final_headers = {backend: {key: float(header[key]) for key in
                                 ("star_age", "Teff", "photosphere_L", "photosphere_r", "center_h1", "center_he3")}
                                 for backend, header in (("baseline", left_header), ("candidate", right_header))}
        finally:
            left_stream.close()
            right_stream.close()
    result = {"pair": pair, "profiles_compared": expected, "zone_counts": sorted(set(zone_counts)),
              "max_abs_age_difference_yr": max_age_difference,
              "max_abs_mass_difference_msun": max_mass_difference,
              "zone_aligned": max_mass_difference <= 1e-11,
              "max_abs_profile_differences": maxima, "final_headers": final_headers}
    output = EXPERIMENT / f"profile_comparison_{pair}.json"
    if output.exists():
        raise FileExistsError(output)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: result[key] for key in ("profiles_compared", "zone_counts", "max_abs_age_difference_yr",
                                                   "max_abs_mass_difference_msun", "zone_aligned")}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pair", choices=("01", "100"))
    compare(parser.parse_args().pair)


if __name__ == "__main__":
    main()
