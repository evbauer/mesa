"""Offline C3 tensor Hermite control using only six native source constraints."""

import math
import sys
from fractions import Fraction
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from freeeos_potential import (
    EVAL_ORDERS, LN10, LOW_ORDERS, local_derivative, query_cells,
    radiation_potential, source_constraints, valid_cells,
)


JET_ORDERS = tuple((first, second) for first in range(4) for second in range(4))
ESTIMATED_ORDERS = tuple(orders for total in range(3, 7) for orders in JET_ORDERS
                         if sum(orders) == total)
ROUTE_GRAPH = {
    orders: tuple((parent, axis) for parent, axis in (
        ((orders[0] - 1, orders[1]), 0), ((orders[0], orders[1] - 1), 1))
        if min(parent) >= 0)
    for orders in ESTIMATED_ORDERS
}


def endpoint_polynomial_matrix():
    return [[Fraction(math.factorial(power), math.factorial(power - derivative))
             * Fraction(endpoint)**(power - derivative) if power >= derivative else Fraction(0)
             for power in range(8)] for endpoint in range(2) for derivative in range(4)]


def analytic_coefficients():
    matrix = endpoint_polynomial_matrix()
    augmented = [row + [Fraction(int(index == column)) for column in range(8)]
                 for index, row in enumerate(matrix)]
    for column in range(8):
        pivot = next(index for index in range(column, 8) if augmented[index][column])
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [entry / divisor for entry in augmented[column]]
        for index in range(8):
            if index != column:
                multiplier = augmented[index][column]
                augmented[index] = [entry - multiplier * pivot_entry
                                    for entry, pivot_entry in zip(augmented[index], augmented[column])]
    return tuple(tuple(augmented[power][8 + constraint] for power in range(8))
                 for constraint in range(8))


EXACT_COEFFICIENTS = analytic_coefficients()
COEFFICIENTS = np.asarray(EXACT_COEFFICIENTS, dtype=float)


def septic_basis(position, spacing, order):
    position = np.atleast_1d(np.asarray(position, dtype=float))
    if position.ndim != 1 or order not in range(4) or not np.isfinite(spacing) or spacing <= 0:
        raise ValueError("Require vector positions, positive spacing, and derivative order 0:3")
    result = np.empty((position.size, 2, 4))
    for endpoint in range(2):
        for derivative in range(4):
            polynomial = np.polynomial.polynomial.polyder(COEFFICIENTS[derivative], order)
            coordinate = position if endpoint == 0 else 1 - position
            sign = 1 if endpoint == 0 else (-1)**(derivative + order)
            result[:, endpoint, derivative] = (sign * spacing**(derivative - order)
                * np.polynomial.polynomial.polyval(coordinate, polynomial))
    for endpoint in range(2):
        exact = position == endpoint
        result[exact] = 0
        result[exact, endpoint, order] = 1
    return result


def construct_jets(low, valid, du, dtau, radius=2):
    if radius not in (2, 3):
        raise ValueError("C3 control requires centered radius 2 or 3")
    if low.shape != (3, 3, *valid.shape) or valid.ndim != 2:
        raise ValueError("Require six source constraints on a two-dimensional native grid")
    nodes = np.full((4, 4, *valid.shape), np.nan)
    masks, routes = {}, {}
    source_valid = valid.copy()
    for orders in LOW_ORDERS:
        source_valid &= np.isfinite(low[orders])
    for orders in LOW_ORDERS:
        nodes[orders] = np.where(source_valid, low[orders], np.nan)
        masks[orders] = source_valid.copy()
    for orders in ESTIMATED_ORDERS:
        estimates, accepted = [], []
        for parent, axis in ROUTE_GRAPH[orders]:
            estimate, mask = local_derivative(nodes[parent], masks[parent],
                (du, dtau)[axis], axis, radius)
            estimates.append(estimate)
            accepted.append(mask & np.isfinite(estimate))
        mask = np.logical_and.reduce(accepted)
        with np.errstate(over="ignore", invalid="ignore"):
            combined = np.mean(estimates, axis=0)
        mask &= np.isfinite(combined)
        nodes[orders] = np.where(mask, combined, np.nan)
        masks[orders] = mask
        routes[orders] = {"parents": ROUTE_GRAPH[orders], "values": tuple(estimates), "valid": mask}
    complete = np.logical_and.reduce([masks[orders] for orders in JET_ORDERS])
    nodes[:, :, ~complete] = np.nan
    return nodes, complete, {"jet_masks": masks, "routes": routes}


def evaluate(nodes, complete, logq, logt, sample_q, sample_t, orders=EVAL_ORDERS,
             cell_indices=None, add_radiation=False):
    sample_q, sample_t = np.broadcast_arrays(np.atleast_1d(sample_q), np.atleast_1d(sample_t))
    if sample_q.ndim != 1 or nodes.shape != (4, 4, len(logq), len(logt)):
        raise ValueError("Require query vectors and a 4x4 nodal jet grid")
    if complete.shape != nodes.shape[2:]:
        raise ValueError("Nodal validity must match grid")
    for grid in (logq, logt):
        if len(grid) < 2 or not np.all(np.isfinite(grid)) or grid[1] <= grid[0] or not np.allclose(
                np.diff(grid), grid[1] - grid[0], rtol=1e-11, atol=1e-13):
            raise ValueError("Require uniform increasing native-coordinate grids")
    if any(orders_pair not in JET_ORDERS for orders_pair in orders):
        raise ValueError("Evaluate only tensor partials 0:3 in each coordinate")
    density_index, temperature_index, inside = query_cells(logq, logt, sample_q, sample_t)
    if cell_indices is not None:
        density_index, temperature_index = map(np.asarray, cell_indices)
        if density_index.shape != sample_q.shape or temperature_index.shape != sample_t.shape:
            raise ValueError("Forced cell indices must match queries")
        if (not np.issubdtype(density_index.dtype, np.integer)
                or not np.issubdtype(temperature_index.dtype, np.integer)
                or np.any(density_index < 0) or np.any(density_index >= len(logq) - 1)
                or np.any(temperature_index < 0) or np.any(temperature_index >= len(logt) - 1)):
            raise ValueError("Forced cell indices outside grid")
        inside &= (sample_q >= logq[density_index]) & (sample_q <= logq[density_index + 1])
        inside &= (sample_t >= logt[temperature_index]) & (sample_t <= logt[temperature_index + 1])
    if not np.all(inside & valid_cells(complete)[density_index, temperature_index]):
        raise ValueError("Query outside complete native-only C3 stencil domain")
    if not np.all(np.isfinite(nodes[:, :, complete])):
        raise ValueError("Complete nodes contain nonfinite jets")
    dq, dt = logq[1] - logq[0], logt[1] - logt[0]
    position_q = (sample_q - logq[density_index]) / dq
    position_t = (sample_t - logt[temperature_index]) / dt
    position_q[sample_q == logq[density_index + 1]] = 1
    position_t[sample_t == logt[temperature_index + 1]] = 1
    density_basis = {order: septic_basis(position_q, dq * LN10, order).reshape(-1, 8)
                     for order in {pair[0] for pair in orders}}
    temperature_basis = {order: septic_basis(position_t, dt * LN10, order).reshape(-1, 8)
                         for order in {pair[1] for pair in orders}}
    corners = np.empty((sample_q.size, 2, 4, 2, 4))
    for density_end in range(2):
        for temperature_end in range(2):
            corners[:, density_end, :, temperature_end, :] = nodes[:, :,
                density_index + density_end, temperature_index + temperature_end].transpose(2, 0, 1)
    corners = corners.reshape(-1, 8, 8)
    result = np.full((4, 4, sample_q.size), np.nan)
    for first, second in orders:
        result[first, second] = np.einsum("ni,nij,nj->n", density_basis[first], corners,
                                         temperature_basis[second], optimize=False)
    if add_radiation:
        radiation = radiation_potential(sample_q, sample_t, orders)
        for pair in orders:
            result[pair] += radiation[pair]
    return result
