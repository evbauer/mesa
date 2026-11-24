"""Classify recorded solver states using the unchanged actual-MESA patchwork."""

import argparse
import csv
import io
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from build_manifest import CODE, WORKSPACE, checksum


def classify(status, fraction):
    if status != 0 or not np.isfinite(fraction) or fraction < -1e-12 or fraction > 1 + 1e-12:
        return "failed"
    if fraction >= 1 - 1e-12:
        return "pure_FreeEOS"
    if fraction > 1e-12:
        return "blend"
    return "inactive"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=WORKSPACE / "freeeos_biquintic_progress/experiments/041_consumer_path_selection")
    arguments = parser.parse_args()
    source = arguments.input.resolve()
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    with source.open(newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"case_id", "trial_id", "model", "stage", "logRho", "logT"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Path CSV lacks required identifiers or coordinates")
        rows = list(reader)
    if not rows:
        raise ValueError("No recorded solver query states")
    coordinates = np.array([[float(row["logRho"]), float(row["logT"])] for row in rows])
    if not np.all(np.isfinite(coordinates)):
        raise ValueError("Nonfinite solver path coordinates")
    unique, inverse = np.unique(coordinates, axis=0, return_inverse=True)
    input_text = "".join(f"{density:.17g} {temperature:.17g}\n" for density, temperature in unique)
    input_path = output / "probe_input.txt"
    table_path = output / "probe.tsv"
    source_identity = checksum(source)
    results_path = output / "results.json"
    if results_path.exists():
        prior = json.loads(results_path.read_text())
        if prior["input_sha256"] != source_identity or checksum(table_path) != prior["probe_output_sha256"]:
            raise ValueError("Saved path audit differs from input or payload; choose a new output directory")
        print(json.dumps(prior["summary"], indent=2))
        return
    if table_path.exists():
        raise ValueError("Partial prior audit exists; inspect it before a new probe")
    input_path.write_text(input_text)
    command = [sys.executable, "-B", str(CODE / "mesa_probe/run_probe.py"), "--output", str(output)]
    with table_path.open("w") as stdout, (output / "probe.stderr").open("w") as stderr:
        completed = subprocess.run(command, input=input_text, text=True, stdout=stdout, stderr=stderr, timeout=180)
    if completed.returncode:
        raise RuntimeError("Actual-MESA path probe failed; inspect its metadata/stderr")
    with table_path.open() as stream:
        header = stream.readline().split()
        data = np.loadtxt(io.StringIO(stream.read()), ndmin=2)
    if data.shape != (len(unique), len(header)) or not np.allclose(data[:, :2], unique, rtol=0, atol=1e-12):
        raise ValueError("Probe count, order or coordinates differ from queries")
    status_index, fraction_index = header.index("ierr"), header.index("frac_FreeEOS")
    successful = data[:, status_index] == 0
    fraction_columns = [index for index, name in enumerate(header) if name.startswith("frac_")]
    if len(fraction_columns) != 8 or not np.allclose(data[successful][:, fraction_columns].sum(axis=1), 1., rtol=0, atol=1e-10):
        raise ValueError("Successful probe fraction inventory is not normalized")
    categories = [classify(record[status_index], record[fraction_index]) for record in data]
    summary = {"recorded_rows": len(rows), "unique_queries": len(unique),
               "unique_categories": {category: categories.count(category) for category in
                                     ("pure_FreeEOS", "blend", "inactive", "failed")}}
    grouped = {}
    fieldnames = list(rows[0]) + ["probe_index", "ierr", "fraction_FreeEOS", "classification"]
    with (output / "classified_paths.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row, probe_index in zip(rows, inverse):
            category = categories[probe_index]
            key = "/".join(row[field] for field in ("case_id", "model", "stage"))
            counts = grouped.setdefault(key, {})
            counts[category] = counts.get(category, 0) + 1
            writer.writerow({**row, "probe_index": int(probe_index), "ierr": int(data[probe_index, status_index]),
                             "fraction_FreeEOS": float(data[probe_index, fraction_index]), "classification": category})
    result = {"scope": "Unchanged actual-MESA selection at recorded solver query states, not candidate runtime outputs or proof of whole continuous-path routing",
              "input": str(source.relative_to(WORKSPACE)), "input_sha256": source_identity,
              "command": command, "script_sha256": checksum(__file__),
              "probe_wrapper_sha256": checksum(CODE / "mesa_probe/run_probe.py"),
              "probe_output_sha256": checksum(table_path), "summary": summary, "by_case_model_stage": grouped}
    results_path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
