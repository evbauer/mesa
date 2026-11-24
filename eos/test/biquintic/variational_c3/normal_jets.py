"""Experimental C3 shared jets with linewise minimum-fourth-derivative normal jets."""

import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.linalg import solve_banded

HERE = Path(__file__).resolve().parent
BIQUINTIC = HERE.parent
sys.path.insert(0, str(BIQUINTIC / 'c3_control'))
sys.path.insert(0, str(BIQUINTIC))

from c3_potential import EXACT_COEFFICIENTS, ESTIMATED_ORDERS, JET_ORDERS, ROUTE_GRAPH
from freeeos_potential import LOW_ORDERS, local_derivative


def fourth_coefficients():
    return tuple(tuple(coefficient * Fraction(math.factorial(power), math.factorial(power - 4))
                       for power, coefficient in enumerate(row) if power >= 4)
                 for row in EXACT_COEFFICIENTS)


FOURTH_EXACT = fourth_coefficients()
FOURTH = np.asarray(FOURTH_EXACT, dtype=float)
STIFFNESS = np.asarray([[float(sum((left[first] * right[second] / Fraction(first + second + 1)
                                  for first in range(4) for second in range(4)), Fraction(0)))
                        for right in FOURTH_EXACT] for left in FOURTH_EXACT])
GAUSS_POSITION, GAUSS_WEIGHT = np.polynomial.legendre.leggauss(4)
GAUSS_POSITION, GAUSS_WEIGHT = (GAUSS_POSITION + 1) / 2, GAUSS_WEIGHT / 2
FOURTH_AT_GAUSS = np.stack([np.polynomial.polynomial.polyval(GAUSS_POSITION, row)
                           for row in FOURTH])


def coarse_bridge_guard(full_valid, coarse_valid):
    if coarse_valid.shape != full_valid[::2, ::2].shape:
        raise ValueError('Require stride-2 coarse validity on the original native grid')
    midpoint_q = full_valid[1:2 * (coarse_valid.shape[0] - 1):2, ::2]
    midpoint_t = full_valid[::2, 1:2 * (coarse_valid.shape[1] - 1):2]
    bridges_q = coarse_valid[:-1] & coarse_valid[1:] & ~midpoint_q
    bridges_t = coarse_valid[:, :-1] & coarse_valid[:, 1:] & ~midpoint_t
    guard = np.ones_like(coarse_valid)
    guard[:-1] &= ~bridges_q
    guard[1:] &= ~bridges_q
    guard[:, :-1] &= ~bridges_t
    guard[:, 1:] &= ~bridges_t
    return guard, {'axis0_invalid_original_midpoints': int(bridges_q.sum()),
                   'axis1_invalid_original_midpoints': int(bridges_t.sum()),
                   'previously_valid_coarse_nodes_excluded': int((coarse_valid & ~guard).sum())}


def known_segments(value, first, second, spacing):
    known = np.zeros((len(value) - 1, 8))
    known[:, 4] = np.diff(value) - spacing * first[:-1] - .5 * spacing**2 * second[:-1]
    known[:, 5] = spacing * np.diff(first) - spacing**2 * second[:-1]
    known[:, 6] = spacing**2 * np.diff(second)
    return known


def line_system(value, first, second, spacing):
    value, first, second = (np.asarray(array, dtype=float) for array in (value, first, second))
    if value.ndim != 1 or first.shape != value.shape or second.shape != value.shape or len(value) < 2:
        raise ValueError('Require matching one-dimensional endpoint data with at least two nodes')
    if not np.isfinite(spacing) or spacing <= 0 or not np.all(np.isfinite([value, first, second])):
        raise ValueError('Require positive finite spacing and finite endpoint data')
    known = known_segments(value, first, second, spacing)
    diagonal = np.zeros(len(value))
    diagonal[:-1] += STIFFNESS[3, 3]
    diagonal[1:] += STIFFNESS[7, 7]
    off_diagonal = np.full(len(value) - 1, STIFFNESS[3, 7])
    rhs = np.zeros(len(value))
    rhs[:-1] -= known @ STIFFNESS[3]
    rhs[1:] -= known @ STIFFNESS[7]
    return diagonal, off_diagonal, rhs


def solve_line(value, first, second, spacing):
    diagonal, off_diagonal, rhs = line_system(value, first, second, spacing)
    bands = np.zeros((3, len(diagonal)))
    bands[0, 1:] = off_diagonal
    bands[1] = diagonal
    bands[2, :-1] = off_diagonal
    scaled_third = solve_banded((1, 1), bands, rhs, check_finite=False)
    residual = diagonal * scaled_third - rhs
    residual[:-1] += off_diagonal * scaled_third[1:]
    residual[1:] += off_diagonal * scaled_third[:-1]
    stationarity = float(np.max(abs(residual)) / max(1., np.max(abs(rhs))))
    return scaled_third / spacing**3, stationarity


def line_energy(value, first, second, third, spacing):
    data = known_segments(np.asarray(value), np.asarray(first), np.asarray(second), spacing)
    data[:, 3] = spacing**3 * np.asarray(third)[:-1]
    data[:, 7] = spacing**3 * np.asarray(third)[1:]
    fourth = data @ FOURTH_AT_GAUSS
    return float(np.sum(fourth**2 * GAUSS_WEIGHT) / spacing**7)


def solve_normal(value, first, second, valid, spacing, axis, radius):
    fields = [np.moveaxis(array, axis, 0) for array in (value, first, second)]
    accepted = np.moveaxis(valid, axis, 0)
    result = np.full_like(fields[0], np.nan)
    mask = np.zeros_like(accepted)
    fd, fd_mask = local_derivative(second, valid, spacing, axis, radius)
    fd, fd_mask = np.moveaxis(fd, axis, 0), np.moveaxis(fd_mask, axis, 0)
    segments = []
    for line_index in range(accepted.shape[1]):
        transitions = np.diff(np.r_[False, accepted[:, line_index], False].astype(int))
        for start, stop in zip(np.flatnonzero(transitions == 1), np.flatnonzero(transitions == -1)):
            if stop - start < 2:
                continue
            data = [field[start:stop, line_index] for field in fields]
            third, stationarity = solve_line(*data, spacing)
            result[start:stop, line_index] = third
            mask[start:stop, line_index] = True
            fd_candidate = np.where(fd_mask[start:stop, line_index], fd[start:stop, line_index], third)
            minimum = line_energy(*data, third, spacing)
            competitor = line_energy(*data, fd_candidate, spacing)
            segments.append({'line': int(line_index), 'start': int(start), 'stop_exclusive': int(stop),
                             'nodes': int(stop - start), 'stationarity_relative': stationarity,
                             'minimum_energy': minimum, 'fd_interior_minimum_elsewhere_energy': competitor,
                             'energy_ratio': minimum / competitor if competitor > 0 else None})
    return np.moveaxis(result, 0, axis), np.moveaxis(mask, 0, axis), segments


def construct_variational_jets(low, valid, du, dtau, radius=2):
    if radius != 2 or low.shape != (3, 3, *valid.shape) or valid.ndim != 2:
        raise ValueError('Require radius-2 six-constraint native grid')
    source_valid = valid.copy()
    for orders in LOW_ORDERS:
        source_valid &= np.isfinite(low[orders])
    nodes = np.full((4, 4, *valid.shape), np.nan)
    masks = {}
    for orders in LOW_ORDERS:
        nodes[orders] = np.where(source_valid, low[orders], np.nan)
        masks[orders] = source_valid.copy()
    normal_segments = {}
    for orders in ESTIMATED_ORDERS:
        if orders in ((3, 0), (0, 3)):
            axis = int(orders == (0, 3))
            first, second = ((1, 0), (2, 0)) if axis == 0 else ((0, 1), (0, 2))
            nodes[orders], masks[orders], segments = solve_normal(
                nodes[0, 0], nodes[first], nodes[second], source_valid, (du, dtau)[axis], axis, radius)
            normal_segments[str(orders)] = segments
            continue
        estimates, accepted = [], []
        for parent, axis in ROUTE_GRAPH[orders]:
            estimate, mask = local_derivative(nodes[parent], masks[parent], (du, dtau)[axis], axis, radius)
            estimates.append(estimate)
            accepted.append(mask & np.isfinite(estimate))
        mask = np.logical_and.reduce(accepted)
        combined = np.mean(estimates, axis=0)
        mask &= np.isfinite(combined)
        nodes[orders] = np.where(mask, combined, np.nan)
        masks[orders] = mask
    complete = np.logical_and.reduce([masks[orders] for orders in JET_ORDERS])
    nodes[:, :, ~complete] = np.nan
    return nodes, complete, {'jet_masks': masks, 'normal_segments': normal_segments}
