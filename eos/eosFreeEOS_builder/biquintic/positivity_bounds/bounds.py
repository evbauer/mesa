"""Numerical Bernstein sign bounds for fixed total-potential Hermite polynomials."""

import math

import numpy as np


DTYPE = np.longdouble
FIELDS = ("Cv_numerator", "P_rho")
SAFETY = 4096
HERMITE = np.array([
    [1, 0, 0, -10, 15, -6], [0, 1, 0, -6, 8, -3],
    [0, 0, .5, -1.5, 1.5, -.5], [0, 0, 0, 10, -15, 6],
    [0, 0, 0, -4, 7, -3], [0, 0, 0, .5, -1, .5],
], dtype=DTYPE)
POWER_TO_BERNSTEIN = np.array([
    [DTYPE(math.comb(row, column)) / math.comb(5, column) if column <= row else 0
     for column in range(6)] for row in range(6)
], dtype=DTYPE)
LEFT = np.array([
    [DTYPE(math.comb(row, column)) / 2**row if column <= row else 0
     for column in range(6)] for row in range(6)
], dtype=DTYPE)
RIGHT = LEFT[::-1, ::-1].copy()
STATUS = {"unresolved": 0, "positive": 1, "negative_witness": 2}


def cell_power(jets, du, dtau):
    """Endpoint-u, endpoint-tau, derivative-u, derivative-tau -> powers x,y."""
    jets = np.asarray(jets, dtype=DTYPE)
    if jets.shape[-4:] != (2, 2, 3, 3) or min(du, dtau) <= 0:
        raise ValueError("Require (...,2,2,3,3) jets and positive native spacings")
    scale = DTYPE(du)**np.arange(3)[:, None] * DTYPE(dtau)**np.arange(3)[None, :]
    scaled = (jets * scale).swapaxes(-3, -2).reshape(*jets.shape[:-4], 6, 6)
    power = HERMITE.T @ scaled @ HERMITE
    magnitude = abs(HERMITE.T) @ abs(scaled) @ abs(HERMITE)
    return power, magnitude


def derivative(power, density_order, temperature_order, du, dtau):
    result = np.zeros_like(power, dtype=DTYPE)
    density_factors = np.array([
        math.factorial(index) // math.factorial(index - density_order)
        for index in range(density_order, 6)
    ], dtype=DTYPE) / DTYPE(du)**density_order
    temperature_factors = np.array([
        math.factorial(index) // math.factorial(index - temperature_order)
        for index in range(temperature_order, 6)
    ], dtype=DTYPE) / DTYPE(dtau)**temperature_order
    result[..., :6-density_order, :6-temperature_order] = (
        power[..., density_order:, temperature_order:]
        * density_factors[:, None] * temperature_factors[None, :])
    return result


def sign_powers(power, magnitude, du, dtau):
    terms = {
        "Cv_numerator": ((1, 0, 1), (-2, 1, 0), (-1, 0, 2), (4, 1, 1), (-4, 2, 0)),
        "P_rho": ((1, 1, 0), (1, 2, 0)),
    }
    result, scales = {}, {}
    for field, orders in terms.items():
        result[field] = sum(weight * derivative(power, density_order, temperature_order, du, dtau)
                            for weight, density_order, temperature_order in orders)
        scales[field] = sum(abs(weight) * derivative(magnitude, density_order, temperature_order, du, dtau)
                            for weight, density_order, temperature_order in orders)
    return result, scales


def bernstein_bounds(power, magnitude, safety=SAFETY):
    coefficients = POWER_TO_BERNSTEIN @ power @ POWER_TO_BERNSTEIN.T
    scale = np.max(POWER_TO_BERNSTEIN @ magnitude @ POWER_TO_BERNSTEIN.T, axis=(-2, -1))
    assembly_margin = safety * DTYPE(np.finfo(DTYPE).eps) * scale
    core_margin = safety * DTYPE(np.finfo(np.float64).eps) * scale
    margin = assembly_margin + core_margin
    return coefficients, margin, assembly_margin, core_margin, scale


def normalize_closed_cell(power, density_stretch, temperature_stretch):
    density_stretch = np.asarray(density_stretch, dtype=DTYPE)
    temperature_stretch = np.asarray(temperature_stretch, dtype=DTYPE)
    scale = (density_stretch[..., None, None]**np.arange(6)[:, None]
             * temperature_stretch[..., None, None]**np.arange(6)[None, :])
    return power * scale


def polynomial_value(power, density_phase, temperature_phase):
    density_phase, temperature_phase = np.broadcast_arrays(
        np.asarray(density_phase, dtype=DTYPE), np.asarray(temperature_phase, dtype=DTYPE))
    result = np.zeros(np.broadcast_shapes(power.shape[:-2], density_phase.shape), dtype=DTYPE)
    for density_index in range(5, -1, -1):
        row = np.zeros_like(result)
        for temperature_index in range(5, -1, -1):
            row = row * temperature_phase + power[..., density_index, temperature_index]
        result = result * density_phase + row
    return result


def bernstein_value(coefficients, density_phase, temperature_phase):
    density_basis = np.array([math.comb(5, index) * density_phase**index * (1-density_phase)**(5-index)
                              for index in range(6)], dtype=DTYPE)
    temperature_basis = np.array([math.comb(5, index) * temperature_phase**index * (1-temperature_phase)**(5-index)
                                  for index in range(6)], dtype=DTYPE)
    return density_basis @ coefficients @ temperature_basis


def subdivide(coefficients):
    return tuple(density_matrix @ coefficients @ temperature_matrix.T
                 for density_matrix in (LEFT, RIGHT) for temperature_matrix in (LEFT, RIGHT))


def adaptive_bound(coefficients, margin, witness=None, max_depth=4, max_nodes=341):
    """Depth-first covering; only an actual evaluated negative point is a failure."""
    if not 0 <= max_depth <= 4 or not 1 <= max_nodes <= 341:
        raise ValueError("Bounded depth 0..4 and node budget 1..341 required")
    stack = [(coefficients, DTYPE(margin), 0, (0., 1., 0., 1.))]
    visited, positive_leaves, unresolved_leaves, deepest = 0, 0, 0, 0
    minimum_lower, minimum_upper = DTYPE(np.inf), DTYPE(np.inf)
    witness_queries = 0
    while stack:
        current, current_margin, depth, box = stack.pop()
        visited += 1
        deepest = max(deepest, depth)
        lower, upper = current.min() - current_margin, current.max() + current_margin
        if lower > 0:
            positive_leaves += 1
            minimum_lower = min(minimum_lower, lower)
            minimum_upper = min(minimum_upper, upper)
            continue
        density_min, density_max, temperature_min, temperature_max = box
        if witness is not None:
            control_density, control_temperature = np.unravel_index(np.argmin(current), (6, 6))
            phases = [(control_density / 5, control_temperature / 5), (.5, .5),
                      (0., 0.), (0., 1.), (1., 0.), (1., 1.)]
            for local_density, local_temperature in dict.fromkeys(phases):
                if bernstein_value(current, DTYPE(local_density), DTYPE(local_temperature)) < -current_margin:
                    density_phase = density_min + (density_max-density_min) * local_density
                    temperature_phase = temperature_min + (temperature_max-temperature_min) * local_temperature
                    evidence = witness(density_phase, temperature_phase)
                    witness_queries += 1
                    if evidence is not None:
                        return {"status": "negative_witness", "visited_nodes": visited,
                                "positive_leaves": positive_leaves, "unresolved_leaves": 0,
                                "maximum_depth": deepest, "lower_bound": float(coefficients.min()-margin),
                                "minimum_leaf_upper_bound": float(upper), "witness": evidence,
                                "witness_core_queries": witness_queries}
        if depth == max_depth or visited + len(stack) + 4 > max_nodes:
            unresolved_leaves += 1
            minimum_lower = min(minimum_lower, lower)
            minimum_upper = min(minimum_upper, upper)
            continue
        density_middle = (density_min + density_max) / 2
        temperature_middle = (temperature_min + temperature_max) / 2
        child_margin = current_margin + SAFETY * DTYPE(np.finfo(DTYPE).eps) * np.max(abs(current))
        boxes = [(density_min, density_middle, temperature_min, temperature_middle),
                 (density_min, density_middle, temperature_middle, temperature_max),
                 (density_middle, density_max, temperature_min, temperature_middle),
                 (density_middle, density_max, temperature_middle, temperature_max)]
        stack.extend((child, child_margin, depth + 1, child_box)
                     for child, child_box in zip(subdivide(current), boxes))
    return {"status": "positive" if unresolved_leaves == 0 else "unresolved",
            "visited_nodes": visited, "positive_leaves": positive_leaves,
            "unresolved_leaves": unresolved_leaves, "maximum_depth": deepest,
            "lower_bound": float(minimum_lower), "minimum_leaf_upper_bound": float(minimum_upper),
            "witness": None, "witness_core_queries": witness_queries}


def domain_masks(logq, logt, corner_ratio):
    density_min, density_max = logq[:-1, None], logq[1:, None]
    temperature_min, temperature_max = logt[None, :-1], logt[None, 1:]
    shape = (len(logq)-1, len(logt)-1)
    branches = {
        "high_temperature": np.broadcast_to(temperature_min >= 6.5, shape),
        "low_q": np.broadcast_to(density_max <= -.5, shape),
        "low_density": density_max.astype(DTYPE) + 2 * temperature_max.astype(DTYPE) - 12 <= DTYPE(-6.5),
    }
    full_own = (temperature_min >= 3.2) & np.logical_or.reduce(list(branches.values()))
    ratio_max = np.maximum.reduce([corner_ratio[:-1, :-1], corner_ratio[1:, :-1],
                                   corner_ratio[:-1, 1:], corner_ratio[1:, 1:]])
    material = ((density_min >= -8) & (density_max <= 5) & (temperature_min >= 3.2)
                & (temperature_max <= 7.8) & np.isfinite(ratio_max) & (ratio_max < 1e-3))
    return full_own, material, branches, ratio_max
