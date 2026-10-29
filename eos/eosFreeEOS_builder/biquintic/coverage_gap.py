"""Audit actual default-MESA contribution where native potential stencils reject."""

import argparse
import io
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

from active_faces import classifications
from freeeos_potential import LN10, checksum, construct_twists, query_cells, read_table, source_constraints, valid_cells
from run_audit import own_use


def main():
    workspace = Path(__file__).resolve().parents[4]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=workspace / "freeeos_report/plots/offgrid_sample_data.npz")
    parser.add_argument("--output", type=Path, default=workspace / "freeeos_biquintic_progress/experiments/019_coverage_gap")
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    table = workspace / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
    _, logq, logt, rows = read_table(table)
    low, valid, _ = source_constraints(rows, logq, logt)
    _, complete, _ = construct_twists(low, valid, (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10)
    with np.load(arguments.input) as saved:
        sample_q, sample_t = saved["logQ"].ravel(), saved["logT"].ravel()
    iq, it, inside = query_cells(logq, logt, sample_q, sample_t)
    accepted = inside & valid_cells(complete)[iq, it]
    rejected_indices = np.flatnonzero(~accepted)
    rejected_q, rejected_t = sample_q[rejected_indices], sample_t[rejected_indices]
    coordinates = np.column_stack((rejected_q + 2 * rejected_t - 12, rejected_t))
    query_text = "".join(f"{density:.17g} {temperature:.17g}\n" for density, temperature in coordinates)
    command = [sys.executable, "-B", str(Path(__file__).parent / "mesa_probe/run_probe.py")]
    completed = subprocess.run(command, input=query_text, capture_output=True, text=True, timeout=900)
    (output / "probe.tsv").write_text(completed.stdout)
    (output / "probe.stderr").write_text(completed.stderr)
    if completed.returncode:
        raise RuntimeError("MESA probe failed during coverage audit")
    probe = np.loadtxt(io.StringIO(completed.stdout), skiprows=1, ndmin=2)
    if probe.shape != (len(rejected_indices), 14) or not np.allclose(probe[:, :2], coordinates, rtol=0, atol=1e-12):
        raise ValueError("Mismatched coverage probe output")
    masks = classifications(probe[:, 2], probe[:, 5])
    subsets = {"all_rejected": np.ones(len(rejected_q), dtype=bool),
               "own_use_rejected": own_use(rejected_q, rejected_t),
               "below_logT3.1": rejected_t < 3.1, "logT3.1_to3.2": (rejected_t >= 3.1) & (rejected_t < 3.2),
               "logT3.2_to3.28": (rejected_t >= 3.2) & (rejected_t < 3.28),
               "logT_atleast3.28": rejected_t >= 3.28}
    result = {"purpose": "Coverage audit only: actual blend fractions at every saved query rejected by radius2 native potential stencils",
              "limitation": "Archive is stratified rather than physical population weighting; no table bug repaired, no candidate state returned, no installed runtime integration",
              "provenance": {"source_sha256": checksum(table), "archive_sha256": checksum(arguments.input), "script_sha256": checksum(__file__),
                             "preprocessing_sha256": checksum(Path(__file__).parent / "freeeos_potential.py"), "probe_command": command},
              "counts": {"saved_queries": len(sample_q), "accepted_potential_queries": int(np.count_nonzero(accepted)), "rejected_potential_queries": len(rejected_indices)},
              "clean_complete_node_logT_min": float(logt[np.flatnonzero(np.any(complete, axis=0))[0]]),
              "subsets": {}}
    for subset, selected in subsets.items():
        subset_result = {"count": int(np.count_nonzero(selected)), "categories": {}}
        for category, category_mask in masks.items():
            selection = selected & category_mask
            category_result = {"count": int(np.count_nonzero(selection))}
            if np.any(selection):
                category_result["logQ_range"] = [float(np.min(rejected_q[selection])), float(np.max(rejected_q[selection]))]
                category_result["logT_range"] = [float(np.min(rejected_t[selection])), float(np.max(rejected_t[selection]))]
            subset_result["categories"][category] = category_result
        result["subsets"][subset] = subset_result
    np.savez_compressed(output / "rejected_queries.npz", archive_flat_indices=rejected_indices, logQ=rejected_q, logT=rejected_t, probe=probe)
    (output / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    os.environ.setdefault("MPLCONFIGDIR", str(output / ".mplconfig"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure, axis = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    inactive = masks["inactive"]
    axis.scatter(coordinates[inactive, 0], rejected_t[inactive], c=".85", s=1, rasterized=True, label="FreeEOS inactive")
    active = masks["active"]
    image = axis.scatter(coordinates[active, 0], rejected_t[active], c=probe[active, 5], vmin=0, vmax=1, cmap="viridis", s=2, rasterized=True)
    axis.set_xlabel("log10 rho [g cm^-3]")
    axis.set_ylabel("log10 T [K]")
    axis.set_title("Saved off-node queries rejected by clean centered potential stencils\nColors = actual default MESA FreeEOS fraction, not a consistency metric", fontsize=10)
    axis.legend(loc="upper left", fontsize=8)
    figure.colorbar(image, ax=axis, label="FreeEOS contribution fraction")
    for extension in ("pdf", "png"):
        figure.savefig(output / f"coverage_gap.{extension}", dpi=180)
    plt.close(figure)
    print(json.dumps({subset: {category: values["count"] for category, values in data["categories"].items()} for subset, data in result["subsets"].items()}, indent=2))


if __name__ == "__main__":
    main()
