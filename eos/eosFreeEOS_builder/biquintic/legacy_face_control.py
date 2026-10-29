"""Matched near-face baseline gaps: installed MESA versus the potential prototype."""

import argparse
import io
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from active_faces import classifications
from freeeos_potential import LN10, checksum, construct_twists, evaluate, read_table, source_constraints, thermodynamic_partials
from run_audit import distribution


def gap(first, second, first_value, second_value):
    scale = np.maximum.reduce([abs(first), abs(second), abs(first_value), abs(second_value), np.ones_like(first)])
    return (first - second) / scale


def main():
    workspace = Path(__file__).resolve().parents[4]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=workspace / "freeeos_biquintic_progress/experiments/014_active_domain/active_faces.npz")
    parser.add_argument("--output", type=Path, default=workspace / "freeeos_biquintic_progress/experiments/018_legacy_face_control")
    parser.add_argument("--random-faces", type=int, default=400)
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    saved = np.load(arguments.input)
    table = workspace / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
    _, logq, logt, rows = read_table(table)
    low, valid, _ = source_constraints(rows, logq, logt)
    nodes, complete, _ = construct_twists(low, valid, (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10)
    generator = np.random.default_rng(2026102818)
    steps = (1e-4, 1e-5, 1e-6, 1e-7)
    batches = []
    queries = []
    for axis, direction in ((0, "u"), (1, "tau")):
        center_probe = saved[f"{direction}_probe"]
        active = classifications(center_probe[:, 2], center_probe[:, 5])["active"]
        indices = np.flatnonzero(active)
        random = generator.choice(indices, min(arguments.random_faces, len(indices)), replace=False)
        extreme = indices[np.argsort(abs(saved[f"{direction}_Cv_t_jump"][indices]))[-25:]]
        selected = np.unique(np.concatenate((random, extreme)))
        sample_q, sample_t = saved[f"{direction}_q"][selected], saved[f"{direction}_t"][selected]
        batch = {"axis": axis, "direction": direction, "q": sample_q, "t": sample_t,
                 "fractions": center_probe[selected, 5], "exact_Cv_t_jump": saved[f"{direction}_Cv_t_jump"][selected], "states": []}
        for step in steps:
            for sign in (-1, 1):
                nearby_q = sample_q + sign * step / LN10 if axis == 0 else sample_q
                nearby_t = sample_t + sign * step / LN10 if axis == 1 else sample_t
                queries.append(np.column_stack((nearby_q + 2 * nearby_t - 12, nearby_t)))
                state = thermodynamic_partials(evaluate(nodes, complete, logq, logt, nearby_q, nearby_t), nearby_q, nearby_t)
                batch["states"].append(state)
        batches.append(batch)
    query_coordinates = np.concatenate(queries)
    query_text = "".join(f"{density:.17g} {temperature:.17g}\n" for density, temperature in query_coordinates)
    (output / "queries.txt").write_text(query_text)
    command = [sys.executable, "-B", str(Path(__file__).parent / "mesa_probe/run_probe.py")]
    completed = subprocess.run(command, input=query_text, capture_output=True, text=True, timeout=600)
    (output / "probe.tsv").write_text(completed.stdout)
    (output / "probe.stderr").write_text(completed.stderr)
    if completed.returncode:
        raise RuntimeError("MESA probe failed")
    probe = np.loadtxt(io.StringIO(completed.stdout), skiprows=1, ndmin=2)
    if probe.shape != (len(query_coordinates), 14) or not np.allclose(probe[:, :2], query_coordinates, rtol=0, atol=1e-12):
        raise ValueError("Mismatched MESA probe output")
    result = {"purpose": "Matched shrinking two-sided finite-distance gaps, not exact legacy polynomial face limits",
              "sampling": "400 random active faces per direction plus25 largest active Cv_t jumps; union, same centers at all steps, separately report pure/blended center status",
              "steps": "natural native normal coordinate displacement; u face moves q only, tau face moves T and rho together at fixedq; returned partials remain physicallnrho/lnT",
              "normalization": "each model uses its own max(1,abs(bothCv),abs(bothresponsepartials)); do not identify the nearfacegap with an error",
              "interpretation": "A continuous baseline partial has shrinking gap; a discontinuous candidate approaches its nonzero exact face jump. This is smoothness, not accuracy or stellar performance.",
              "provenance": {"source_sha256": checksum(table), "script_sha256": checksum(__file__), "preprocessing_sha256": checksum(Path(__file__).parent / "freeeos_potential.py"), "command": command}, "faces": {}}
    archive = {}
    offset = 0
    for batch in batches:
        direction, count = batch["direction"], len(batch["q"])
        local = probe[offset:offset + len(steps) * 2 * count].reshape(len(steps), 2, count, 14)
        offset += len(steps) * 2 * count
        if np.any(local[..., 2] != 0):
            raise ValueError("Failed near-face MESA queries")
        categories = {"all_active_centers": np.ones(count, dtype=bool), "pure_FreeEOS_centers": batch["fractions"] >= 1 - 1e-12,
                      "blended_centers": batch["fractions"] < 1 - 1e-12}
        face_result = {"count": count, "categories": {}}
        exact_jump = batch["exact_Cv_t_jump"]
        gaps = {"potential": {}, "installed_MESA": {}}
        for step_index, step in enumerate(steps):
            left, right = batch["states"][2 * step_index:2 * step_index + 2]
            for field, column in (("Cv_x", 12), ("Cv_t", 13)):
                gaps["potential"].setdefault(field, []).append(gap(left[field], right[field], left["Cv"], right["Cv"]))
                gaps["installed_MESA"].setdefault(field, []).append(gap(local[step_index, 0, :, column], local[step_index, 1, :, column], local[step_index, 0, :, 11], local[step_index, 1, :, 11]))
        for category, mask in categories.items():
            category_result = {"count": int(np.count_nonzero(mask))}
            if np.any(mask):
                category_result["exact_potential_Cv_t_jump"] = distribution(exact_jump[mask])
                category_result["near_face_gaps"] = {model: {field: {str(step): distribution(values[mask]) for step, values in zip(steps, field_gaps)} for field, field_gaps in model_gaps.items()} for model, model_gaps in gaps.items()}
            face_result["categories"][category] = category_result
        archive.update({f"{direction}_q": batch["q"], f"{direction}_t": batch["t"], f"{direction}_probe": local,
                        f"{direction}_center_frac": batch["fractions"], f"{direction}_exact_Cv_t_jump": exact_jump})
        for model, model_gaps in gaps.items():
            archive.update({f"{direction}_{model}_{field}_gaps": np.asarray(values) for field, values in model_gaps.items()})
        result["faces"][direction] = face_result
    np.savez_compressed(output / "gaps.npz", **archive)
    (output / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({direction: data["categories"]["pure_FreeEOS_centers"] for direction, data in result["faces"].items()}, indent=2))


if __name__ == "__main__":
    main()
