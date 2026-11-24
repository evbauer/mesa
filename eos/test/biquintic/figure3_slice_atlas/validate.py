"""Focused archive, face, coordinate, and within-cell derivative checks."""

import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
NEW = ROOT / "freeeos_biquintic_progress/experiments/071_figure3_slice_atlas"
sys.path.insert(0, str(HERE))
from atlas import CACHE, SOURCE, load_construction
sys.path.insert(0, str(HERE.parent))
from freeeos_potential import LN10, checksum
from active_slice_extended.run_extended import potential_state


def main():
    selection = json.loads((NEW / "selection.json").read_text())
    routing = json.loads((NEW / "routing_summary.json").read_text())
    metadata = json.loads((NEW / "curve_metadata.json").read_text())
    probe_metadata = json.loads((NEW / "probe_metadata.json").read_text())
    if not probe_metadata["installed_cache_inventory_unchanged"] or probe_metadata["private_replaced_cache_links"]:
        raise ValueError("MESA probe cache safety condition failed")
    if checksum(SOURCE) != selection["source_sha256"] or checksum(CACHE) != selection["cache_sha256"]:
        raise ValueError("Source or installed cache changed")
    logq, logt, nodes, complete = load_construction()
    checks = {}
    for case in selection["cases"][1:]:
        name = case["id"]
        path = NEW / f"{name}_curves.npz"
        if checksum(path) != metadata[name]["curves_sha256"]:
            raise ValueError(f"Saved {name} curves changed")
        data = np.load(path)
        if not np.allclose(data["logRho"], data["logQ"] + 2 * data["logT"] - 12, rtol=0, atol=3e-14):
            raise ValueError(f"{name} physical density coordinate mismatch")
        anchor = int(np.argmin(abs(data["logQ"] - case["logQ"])))
        if abs(data["logQ"][anchor] - case["logQ"]) > 1e-12:
            raise ValueError(f"{name} query misses anchor")
        if case["orientation"] == "tau":
            direct = data["lower_Cv_t"][anchor] - data["upper_Cv_t"][anchor]
            error = abs(direct - (case["Cv_t_left"] - case["Cv_t_right"]))
            if error > 1e-7:
                raise ValueError(f"{name} exact tau trace mismatch: {error}")
        else:
            duplicates = np.flatnonzero(abs(data["trace_logQ"] - case["logQ"]) < 1e-12)
            if len(duplicates) != 2:
                raise ValueError(f"{name} u face must have two exact endpoint traces")
            measured = data["trace_single_Cv_t"][duplicates]
            error = float(np.max(abs(measured - [case["Cv_t_left"], case["Cv_t_right"]])))
            if error > 1e-7:
                raise ValueError(f"{name} exact u trace mismatch: {error}")
        if routing[name]["supported_count"] != metadata[name]["source_support_count"]:
            raise ValueError(f"{name} source support count mismatch")
        center_q = case["logQ"] + .005
        center_t = case["logT"] + (.003 if case["orientation"] == "tau" else 0)
        center = potential_state(nodes, complete, logq, logt, np.array([center_q]), np.array([center_t]))
        partial_errors = {}
        for axis in ("x", "t"):
            step = 1e-5
            shift_q = step / LN10 if axis == "x" else -2 * step / LN10
            shift_t = 0 if axis == "x" else step / LN10
            plus = potential_state(nodes, complete, logq, logt, np.array([center_q + shift_q]),
                                   np.array([center_t + shift_t]))
            minus = potential_state(nodes, complete, logq, logt, np.array([center_q - shift_q]),
                                    np.array([center_t - shift_t]))
            numerical = (plus["Cv"][0] - minus["Cv"][0]) / (2 * step)
            analytical = center["Cv_" + axis][0]
            scale = max(abs(center["Cv"][0]), abs(analytical))
            partial_errors[axis] = abs(numerical - analytical) / scale
            if partial_errors[axis] > 2e-5:
                raise ValueError(f"{name} {axis} within-cell derivative mismatch")
        checks[name] = {"face_Cv_t_error_cgs": float(error), "within_cell_scaled_Cv_partial_error": partial_errors,
                        "routing_pure": routing[name]["pure_FreeEOS_count"],
                        "source_supported": routing[name]["supported_count"]}
    result = {"all_passed": True, "construction_source_sha256": checksum(SOURCE),
              "installed_cache_sha256": checksum(CACHE), "routing_probe_sha256": checksum(NEW / "routing_probe.tsv"),
              "routing_queries_sha256": checksum(NEW / "routing_queries.txt"),
              "probe_installed_cache_inventory_unchanged": True, "cases": checks}
    (NEW / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
