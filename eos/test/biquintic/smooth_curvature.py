"""One-dimensional C3 quintic control: quantify required curvature reconciliation."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import solve_banded

from freeeos_potential import LN10, read_table, source_constraints
from run_audit import distribution, own_use
from validate_table import material_nodes


def reconcile_curvatures(values, slopes, curvatures, spacing):
    count = len(values)
    if count < 3:
        return curvatures.copy()
    right = 20 * (values[:-2] - 2 * values[1:-1] + values[2:]) / spacing**2
    right += 8 * (slopes[:-2] - slopes[2:]) / spacing
    right[0] += curvatures[0]
    right[-1] += curvatures[-1]
    band = np.zeros((3, count - 2))
    band[0, 1:] = -1
    band[1] = 6
    band[2, :-1] = -1
    result = curvatures.copy()
    result[1:-1] = solve_banded((1, 1), band, right)
    return result


def normal_third_jump(values, slopes, curvatures, spacing):
    return ((-60 * values[:-2] + 120 * values[1:-1] - 60 * values[2:]) / spacing**3
            + (-24 * slopes[:-2] + 24 * slopes[2:]) / spacing**2
            + (-3 * curvatures[:-2] + 18 * curvatures[1:-1] - 3 * curvatures[2:]) / spacing)


def main():
    workspace = Path(__file__).resolve().parents[4]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesa", type=Path, default=workspace / "mesa")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    _, logq, logt, rows = read_table(arguments.mesa / "data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data")
    nodes, valid, _ = source_constraints(rows, logq, logt)
    query_q, query_t = np.broadcast_arrays(logq[:, None], logt[None, :])
    selected = valid & material_nodes(rows, logq, logt) & own_use(query_q, query_t)
    result = {"purpose": "1D feasibility control, not implemented 2D EOS: preserve f and first derivative, reconcile second derivative to make quintic C3 along each native line",
              "boundary": "retain supplied second derivatives at ends of each clean contiguous source run",
              "limitation": "separate line solves are not a globally consistent 2D potential or production recommendation", "axes": {}}
    for axis, direction, first_orders, second_orders in ((0, "u", (1, 0), (2, 0)), (1, "tau", (0, 1), (0, 2))):
        spacing = ((logq[1] - logq[0]) if axis == 0 else (logt[1] - logt[0])) * LN10
        free = np.moveaxis(nodes[0, 0], axis, 0)
        slopes = np.moveaxis(nodes[first_orders], axis, 0)
        supplied = np.moveaxis(nodes[second_orders], axis, 0)
        accepted = np.moveaxis(valid, axis, 0)
        correction = np.full_like(supplied, np.nan)
        old_jumps, new_jumps, jump_scales = [], [], []
        for line in range(free.shape[1]):
            changes = np.flatnonzero(np.diff(np.r_[False, accepted[:, line], False]))
            for start, end in zip(changes[::2], changes[1::2]):
                if end - start < 5:
                    continue
                values, gradients, curvatures = free[start:end, line], slopes[start:end, line], supplied[start:end, line]
                revised = reconcile_curvatures(values, gradients, curvatures, spacing)
                correction[start:end, line] = revised - curvatures
                old_jumps.append(normal_third_jump(values, gradients, curvatures, spacing))
                new_jumps.append(normal_third_jump(values, gradients, revised, spacing))
                jump_scales.append(np.maximum(abs(values[1:-1]) / spacing**3, 1))
        correction = np.moveaxis(correction, 0, axis)
        usable = selected & np.isfinite(correction)
        characteristic = np.maximum.reduce([abs(nodes[second_orders]), abs(nodes[first_orders]), np.ones_like(correction)])
        normalized = correction[usable] / characteristic[usable]
        worst_flat = np.flatnonzero(usable)[int(np.argmax(abs(normalized)))]
        worst_q, worst_t = np.unravel_index(worst_flat, usable.shape)
        axis_result = {"selected_nodes": int(usable.sum()), "curvature_correction_over_max_curvature_slope_floor1": distribution(normalized),
                       "old_normal_third_jump_over_F_h_cubed": distribution(np.concatenate(old_jumps) / np.concatenate(jump_scales)),
                       "reconciled_normal_third_jump_over_F_h_cubed": distribution(np.concatenate(new_jumps) / np.concatenate(jump_scales)),
                       "worst": {"logQ": float(logq[worst_q]), "logT": float(logt[worst_t]), "supplied_curvature": float(nodes[second_orders][worst_q, worst_t]),
                                 "correction": float(correction[worst_q, worst_t])}}
        if axis == 0:
            physical_change = correction / nodes[1, 0]
            axis_result["induced_nodal_chiRho_absolute_change"] = distribution(physical_change[usable])
            axis_result["induced_nodal_chiRho_relative_change"] = distribution(physical_change[usable] / rows[:, :, 4][usable])
        else:
            physical_cv_change = -correction / (10**logt[None, :])
            axis_result["induced_nodal_Cv_relative_change_if_other_constraints_fixed"] = distribution(physical_cv_change[usable] / rows[:, :, 7][usable])
        result["axes"][direction] = axis_result
    (output / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["axes"], indent=2))


if __name__ == "__main__":
    main()
