"""Local source-line compatibility and physical face-jump diagnostics."""

import sys
from pathlib import Path

import numpy as np


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from freeeos_potential import evaluate, thermodynamic_partials


def quintic_third_jump(values, first, second, spacing):
    """Right-minus-left third derivative at the middle of three knots."""
    values = np.asarray(values, dtype=np.longdouble)
    first = np.asarray(first, dtype=np.longdouble)
    second = np.asarray(second, dtype=np.longdouble)
    step = np.longdouble(spacing)
    return float(60 * (values[2] - 2 * values[1] + values[0]) / step**3
                 - 24 * (first[2] - first[0]) / step**2
                 + 3 * (second[2] - 6 * second[1] + second[0]) / step)


def line_compatibility(low, direction, normal_index, tangent_index, spacing):
    offsets = np.arange(normal_index - 1, normal_index + 2)
    if direction == "u":
        triples = (low[order, 0, offsets, tangent_index] for order in range(3))
    elif direction == "tau":
        triples = (low[0, order, tangent_index, offsets] for order in range(3))
    else:
        raise ValueError("Direction must be u or tau")
    return quintic_third_jump(*triples, spacing)


def face_state(nodes, complete, logq, logt, case, split):
    sample_q = np.array([case["logQ"]])
    sample_t = np.array([case["logT"]])
    sides = {}
    for side in ("left", "right"):
        cell = case[f"{side}_cell"]
        native = evaluate(nodes, complete, logq, logt, sample_q, sample_t,
                          cell_indices=(np.array([cell[0]]), np.array([cell[1]])),
                          add_radiation=split)
        state = thermodynamic_partials(native, sample_q, sample_t)
        normal = native[3, 0, 0] if case["direction"] == "u" else native[0, 3, 0]
        sides[side] = {field: float(state[field][0]) for field in ("Cv", "Cv_t", "Cv_x")}
        sides[side]["normal_third"] = float(normal)
    scale = max(1.0, *(abs(sides[side][field]) for side in sides for field in ("Cv", "Cv_t")))
    jump = sides["left"]["Cv_t"] - sides["right"]["Cv_t"]
    normal_jump = sides["left"]["normal_third"] - sides["right"]["normal_third"]
    temperature = 10.0**case["logT"]
    expected_t = (8.0 if case["direction"] == "u" else -1.0) * normal_jump / temperature
    expected_x = (-4.0 if case["direction"] == "u" else 0.0) * normal_jump / temperature
    return {
        "left": sides["left"], "right": sides["right"], "Cv_t_jump": jump,
        "Cv_x_jump": sides["left"]["Cv_x"] - sides["right"]["Cv_x"],
        "scale": scale, "scaled_Cv_t_jump": jump / scale,
        "normal_third_jump": normal_jump, "shear_expected_Cv_t_jump": expected_t,
        "shear_expected_Cv_x_jump": expected_x,
    }


def knot_line(case, low, logq, logt, seed_scale):
    direction = case["direction"]
    left = case["left_cell"]
    normal_index = case["right_cell"][0 if direction == "u" else 1]
    tangent_index = left[1 if direction == "u" else 0]
    coordinate = logt if direction == "u" else logq
    spacing = (logq[1] - logq[0] if direction == "u" else logt[1] - logt[0]) * np.log(10)
    result = []
    for endpoint in (tangent_index, tangent_index + 1):
        third_right_minus_left = line_compatibility(low, direction, normal_index, endpoint, spacing)
        temperature = 10.0**(coordinate[endpoint] if direction == "u" else logt[normal_index])
        physical_jump = (-8.0 if direction == "u" else 1.0) * third_right_minus_left / temperature
        result.append({"tangent_index": int(endpoint), "tangent_coordinate": float(coordinate[endpoint]),
                       "normal_third_right_minus_left": third_right_minus_left,
                       "Cv_t_left_minus_right": physical_jump,
                       "Cv_t_jump_over_seed_scale": physical_jump / seed_scale})
    return result


def nodal_closure(rows, source, indices):
    density = source["rho"][indices]
    temperature = source["T"][indices]
    pressure_scale = source["P"][indices] / density
    mixed = pressure_scale * rows[..., 5][indices]
    pairs = {
        "mixed_energy": (mixed, pressure_scale - density * rows[..., 8][indices]),
        "mixed_entropy": (mixed, -density * temperature * rows[..., 10][indices]),
        "Cv_entropy": (rows[..., 7][indices], temperature * rows[..., 9][indices]),
    }
    return {name: float(np.max(np.abs(first - second) /
                            np.maximum(np.maximum(np.abs(first), np.abs(second)), np.finfo(float).tiny)))
            for name, (first, second) in pairs.items()}


def face_source_nodes(case):
    direction = case["direction"]
    normal = case["right_cell"][0 if direction == "u" else 1]
    tangent = case["left_cell"][1 if direction == "u" else 0]
    normal_indices = np.repeat(np.arange(normal - 1, normal + 2), 2)
    tangent_indices = np.tile(np.arange(tangent, tangent + 2), 3)
    return (normal_indices, tangent_indices) if direction == "u" else (tangent_indices, normal_indices)
