"""Separate twist-dependent interior-face jumps from fixed nodal C3 obstacles."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from freeeos_potential import LN10, construct_twists, evaluate, read_table, source_constraints, thermodynamic_partials, valid_cells
from run_audit import distribution, own_use
from validate_table import material_nodes


def main():
    workspace = Path(__file__).resolve().parents[4]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesa", type=Path, default=workspace / "mesa")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    _, logq, logt, rows = read_table(arguments.mesa / "data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data")
    low, valid, _ = source_constraints(rows, logq, logt)
    du, dtau = (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10
    schemes = {f"radius{radius}_{route}": construct_twists(low, valid, du, dtau, radius, route)
               for radius, route in ((1, "mean"), (2, "mean"), (3, "mean"), (2, "first"), (2, "second"))}
    common_cells = valid_cells(material_nodes(rows, logq, logt))
    for _, complete, _ in schemes.values():
        common_cells &= valid_cells(complete)
    result = {"purpose": "normal third potential derivatives at nodal faces are fixed by f/first/second constraints; no twist choice can enforce their C3 continuity",
              "metric": "jump divided by max(1,response magnitude,both one-sided response derivatives)", "faces": {}}
    generator = np.random.default_rng(20261030)
    for axis, direction in ((0, "u"), (1, "tau")):
        adjacent = common_cells[:-1] & common_cells[1:] if axis == 0 else common_cells[:, :-1] & common_cells[:, 1:]
        density_index, temperature_index = np.nonzero(adjacent)
        choice = generator.choice(len(density_index), min(10000, len(density_index)), replace=False)
        density_index, temperature_index = density_index[choice], temperature_index[choice]
        fractions = generator.uniform(.1, .9, len(choice))
        direction_result = {}
        for mode, fraction in (("nodal_tangent", np.zeros(len(choice))), ("offnode_tangent", fractions)):
            if axis == 0:
                sample_q = logq[density_index + 1]
                sample_t = logt[temperature_index] + fraction * (logt[1] - logt[0])
                right_indices = density_index + 1, temperature_index
            else:
                sample_q = logq[density_index] + fraction * (logq[1] - logq[0])
                sample_t = logt[temperature_index + 1]
                right_indices = density_index, temperature_index + 1
            selected = own_use(sample_q, sample_t)
            sample_q, sample_t = sample_q[selected], sample_t[selected]
            left_indices = density_index[selected], temperature_index[selected]
            right_indices = right_indices[0][selected], right_indices[1][selected]
            scheme_results = {}
            third_jumps = []
            for name, (nodes, complete, _) in schemes.items():
                left = evaluate(nodes, complete, logq, logt, sample_q, sample_t, cell_indices=left_indices)
                right = evaluate(nodes, complete, logq, logt, sample_q, sample_t, cell_indices=right_indices)
                left_state, right_state = thermodynamic_partials(left, sample_q, sample_t), thermodynamic_partials(right, sample_q, sample_t)
                normal_order = (3, 0) if axis == 0 else (0, 3)
                third_jumps.append(left[normal_order] - right[normal_order])
                metrics = {}
                for field in ("Cv", "Cp", "chiRho", "chiT", "gamma1", "grad_ad"):
                    for suffix in ("_x", "_t"):
                        first, second = left_state[field + suffix], right_state[field + suffix]
                        scale = np.maximum.reduce([abs(first), abs(second), abs(left_state[field]), abs(right_state[field]), np.ones(len(sample_q))])
                        jump = (first - second) / scale
                        worst = int(np.argmax(abs(jump)))
                        metrics[field + suffix] = {**distribution(jump), "worst": {"logQ": float(sample_q[worst]), "logT": float(sample_t[worst]), "signed_jump": float(jump[worst])}}
                scheme_results[name] = metrics
            reference_jump = third_jumps[0]
            third_scale = np.maximum(abs(reference_jump), 1)
            direction_result[mode] = {"count": len(sample_q), "schemes": scheme_results,
                                      "normal_third_jump_change_between_schemes": [distribution((jump - reference_jump) / third_scale) for jump in third_jumps[1:]]}
        result["faces"][direction] = direction_result
    (output / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({direction: {mode: {"count": data["count"], "max_normal_jump_scheme_variation": max(value["max"] for value in data["normal_third_jump_change_between_schemes"])}
                                    for mode, data in modes.items()} for direction, modes in result["faces"].items()}, indent=2))


if __name__ == "__main__":
    main()
