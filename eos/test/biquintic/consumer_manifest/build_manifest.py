"""Freeze active-face cases before controlled-solver outcomes are measured."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np


CODE = Path(__file__).resolve().parents[1]
WORKSPACE = CODE.parents[3]
sys.path.insert(0, str(CODE))

from freeeos_potential import LN10, checksum, construct_twists, evaluate, read_table, source_constraints, thermodynamic_partials, valid_cells


def select_indices(jumps, eligible):
    available = np.flatnonzero(eligible)
    if len(available) < 24:
        raise ValueError("Insufficient eligible faces for a balanced fixed sample")
    ordered = available[np.argsort(np.abs(jumps[available]), kind="stable")]
    ordinary = ordered[np.rint(np.linspace(.2, .7, 6) * (len(ordered) - 1)).astype(int)]
    tail = ordered[-6:][::-1]
    if len(np.unique(np.concatenate((ordinary, tail)))) != 12:
        raise ValueError("Selected ordinary and tail cases must be distinct")
    return {"ordinary": ordinary, "tail": tail}


def composition_metadata(header):
    return {"X": float(header[1]), "Z": float(header[2])}


def face_indices(logq, logt, sample_q, sample_t, direction):
    density_cell = int(np.searchsorted(logq, sample_q, side="right") - 1)
    temperature_cell = int(np.searchsorted(logt, sample_t, side="right") - 1)
    if direction == "u":
        knot = int(np.argmin(abs(logq - sample_q)))
        if abs(logq[knot] - sample_q) > 1e-12:
            raise ValueError("Archive u-face coordinate is not on a source knot")
        return [knot - 1, temperature_cell], [knot, temperature_cell]
    if direction == "tau":
        knot = int(np.argmin(abs(logt - sample_t)))
        if abs(logt[knot] - sample_t) > 1e-12:
            raise ValueError("Archive tau-face coordinate is not on a source knot")
        return [density_cell, knot - 1], [density_cell, knot]
    raise ValueError("Unknown face direction")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=WORKSPACE / "freeeos_biquintic_progress/experiments/038_solver_manifest")
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    archive = WORKSPACE / "freeeos_biquintic_progress/experiments/014_active_domain/active_faces.npz"
    table = WORKSPACE / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
    header, logq, logt, rows = read_table(table)
    low, valid, _ = source_constraints(rows, logq, logt, split_radiation=True, minimum_logT=3.2)
    nodes, complete, _ = construct_twists(low, valid, (logq[1] - logq[0]) * LN10,
                                         (logt[1] - logt[0]) * LN10, radius=2, route="mean")
    cells = valid_cells(complete)
    cases = []
    populations = {}
    with np.load(archive) as saved:
        for direction in ("u", "tau"):
            probe = saved[f"{direction}_probe"]
            eligible = (probe[:, 2] == 0) & (probe[:, 5] >= 1 - 1e-12)
            jumps = saved[f"{direction}_Cv_t_jump"]
            selections = select_indices(jumps, eligible)
            populations[direction] = int(np.count_nonzero(eligible))
            for category, indices in selections.items():
                for number, archive_index in enumerate(indices, start=1):
                    sample_q = float(saved[f"{direction}_q"][archive_index])
                    sample_t = float(saved[f"{direction}_t"][archive_index])
                    left, right = face_indices(logq, logt, sample_q, sample_t, direction)
                    if not cells[tuple(left)] or not cells[tuple(right)]:
                        raise ValueError("Selected face lacks complete baseline support")
                    states = []
                    for indices_side in (left, right):
                        cell = tuple(np.array([index]) for index in indices_side)
                        state = thermodynamic_partials(evaluate(nodes, complete, logq, logt,
                            np.array([sample_q]), np.array([sample_t]), add_radiation=True,
                            cell_indices=cell), np.array([sample_q]), np.array([sample_t]))
                        states.append({field: float(state[field][0]) for field in ("Cv", "Cv_t", "Cv_x")})
                    scale = max(1., *(abs(state[field]) for state in states for field in ("Cv", "Cv_t")))
                    cases.append({"case_id": f"{direction}_{category}_{number:02d}",
                        "direction": direction, "category": category, "archive_index": int(archive_index),
                        "logQ": sample_q, "logT": sample_t, "logRho": sample_q + 2 * sample_t - 12,
                        "left_cell": left, "right_cell": right,
                        "fraction_FreeEOS_at_seed": float(probe[archive_index, 5]),
                        "archived_scaled_Cv_t_jump": float(jumps[archive_index]),
                        "split_scaled_Cv_t_jump": (states[0]["Cv_t"] - states[1]["Cv_t"]) / scale,
                        "split_left": states[0], "split_right": states[1]})
    result = {"schema": "FREEEOS_CONSUMER_MANIFEST_V2", "status": "case selection frozen before solver outcomes; composition labels corrected without changing cases",
        "supersedes": {"file": "queries.json", "sha256": "74f5079a7b5830abeae2947c4063d78e214a9dde2be43dac59f507d6b6007eb8",
                       "reason": "V1 metadata reversed X and Z header labels; correct source table and all 24 selected cases are unchanged"},
        "selection": "six ordinary per direction at absolute archived Cv_t-jump rank quantiles .2,.3,.4,.5,.6,.7; six largest absolute jumps per direction; stable ranking; pure-FreeEOS seed fraction >=1-1e-12",
        "domain_limit": "Actual MESA fraction is known at seed only. Solver paths must be independently classified; support and seed classification do not imply whole-path pure FreeEOS.",
        "normalization": "signed left-minus-right Cv_t jump / max(1,abs(left/right Cv),abs(left/right Cv_t)); physical natural-log T at fixed rho",
        "candidate": {**composition_metadata(header), "composition": "audited normalized mass_frac_solar.txt metal pattern",
            "radius": 2, "route": "mean", "split_radiation": True, "minimum_logT": 3.2, "boundary_stencil": "centered"},
        "provenance": {"source_table": str(table.relative_to(WORKSPACE)), "source_sha256": checksum(table),
            "archive": str(archive.relative_to(WORKSPACE)), "archive_sha256": checksum(archive),
            "core_sha256": checksum(CODE / "freeeos_potential.py"), "script_sha256": checksum(__file__)},
        "eligible_counts": populations, "case_count": len(cases), "cases": cases}
    previous = output / "queries.json"
    if previous.exists() and json.loads(previous.read_text())["cases"] != cases:
        raise ValueError("Metadata correction must not change any frozen case")
    destination = output / "queries_v2.json"
    if destination.exists() and json.loads(destination.read_text()) != result:
        raise ValueError("Refusing to overwrite a frozen manifest with a changed selection or provenance")
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"manifest": str(destination), "sha256": checksum(destination), "case_count": len(cases), "eligible_counts": populations}, indent=2))


if __name__ == "__main__":
    main()
