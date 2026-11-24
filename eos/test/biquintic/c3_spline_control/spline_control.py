"""Patch-local degree-five double-interior-knot C3 potential experiment."""

import sys
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import BSpline
from scipy.linalg import lstsq, null_space
from scipy.sparse import csr_matrix, kron, vstack

HERE = Path(__file__).resolve().parent
BIQUINTIC = HERE.parent
sys.path.insert(0, str(BIQUINTIC))

from freeeos_potential import EVAL_ORDERS, LOW_ORDERS


class RankDeficientFit(RuntimeError):
    def __init__(self, diagnostics):
        self.diagnostics = diagnostics
        super().__init__(f"Rank deficient QR: rank={diagnostics['rank']}, unknowns={diagnostics['unknowns']}")


def knots_for(grid, multiplicity=2):
    grid = np.asarray(grid)
    if grid.ndim != 1 or len(grid) < 3 or not np.all(np.isfinite(grid)) or not np.all(np.diff(grid) > 0):
        raise ValueError('Require strictly increasing finite one-dimensional grid with at least three nodes')
    if multiplicity not in (1, 2):
        raise ValueError('Only single C4 or double C3 interior knots are supported')
    return np.r_[np.repeat(grid[0], 6), np.repeat(grid[1:-1], multiplicity), np.repeat(grid[-1], 6)]


def spline_basis(grid, query, order=0, multiplicity=2):
    if order not in range(6):
        raise ValueError('Derivative order must be zero through five')
    knots = knots_for(grid, multiplicity)
    count = len(knots) - 6
    return BSpline(knots, np.eye(count), 5, extrapolate=False)(query, nu=order)


def evaluate_fit(coefficients, grid_u, grid_tau, sample_u, sample_tau,
                 orders=EVAL_ORDERS, multiplicity=2, chunk_size=512):
    sample_u, sample_tau = np.broadcast_arrays(np.atleast_1d(sample_u), np.atleast_1d(sample_tau))
    if sample_u.ndim != 1 or chunk_size < 1:
        raise ValueError('Require vector queries and positive chunk size')
    expected = (len(knots_for(grid_u, multiplicity)) - 6, len(knots_for(grid_tau, multiplicity)) - 6)
    if coefficients.shape != expected or not np.all(np.isfinite(coefficients)):
        raise ValueError('Invalid coefficient dimensions or nonfinite fit')
    result = np.full((4, 4, len(sample_u)), np.nan)
    for start in range(0, len(sample_u), chunk_size):
        stop = min(len(sample_u), start + chunk_size)
        basis_u = {order: spline_basis(grid_u, sample_u[start:stop], order, multiplicity)
                   for order in {pair[0] for pair in orders}}
        basis_tau = {order: spline_basis(grid_tau, sample_tau[start:stop], order, multiplicity)
                     for order in {pair[1] for pair in orders}}
        for first, second in orders:
            result[first, second, start:stop] = np.einsum('qi,ij,qj->q', basis_u[first], coefficients,
                                                        basis_tau[second], optimize=True)
    return result


def thermodynamic_combinations(low):
    ft = low[0, 1] - 2 * low[1, 0]
    ftt = low[0, 2] - 4 * low[1, 1] + 4 * low[2, 0]
    return {'E': low[0, 0] - ft, 'TS': -ft, 'P_over_rho': low[1, 0],
            'P_rho': low[1, 0] + low[2, 0], 'P_over_rho_chiT': low[1, 1] - 2 * low[2, 0],
            'TCv': ft - ftt}


def weighted_system(low, grid_u, grid_tau, response_weight, multiplicity=2):
    if response_weight not in (.01, .1):
        raise ValueError('Bounded control permits only response weights .01 and .1')
    if low.shape != (3, 3, len(grid_u), len(grid_tau)) or not all(np.all(np.isfinite(low[pair])) for pair in LOW_ORDERS):
        raise ValueError('Require all six finite native supplied constraints')
    matrices = {orders: kron(csr_matrix(spline_basis(grid_tau, grid_tau, orders[1], multiplicity)),
                             csr_matrix(spline_basis(grid_u, grid_u, orders[0], multiplicity)), format='csr')
                for orders in LOW_ORDERS}
    combinations = thermodynamic_combinations(matrices)
    values = thermodynamic_combinations(low)
    value_scale = float(np.median(abs(low[0, 0])))
    if value_scale <= 0 or not np.isfinite(value_scale):
        raise ValueError('Invalid potential scale')
    matrices_weighted, targets, relative_floors = [], [], {}
    for field, matrix in combinations.items():
        value = values[field].ravel(order='F')
        floor = .01 * np.median(abs(value))
        relative_scale = np.maximum(abs(value), floor)
        if np.any(relative_scale <= 0):
            raise ValueError(f'Zero relative scale for {field}')
        importance = response_weight if field in ('P_rho', 'P_over_rho_chiT', 'TCv') else 1
        weights = importance * value_scale / relative_scale
        matrices_weighted.append(matrix.multiply(weights[:, None]))
        targets.append(value / value_scale * weights)
        relative_floors[field] = float(floor)
    return vstack(matrices_weighted, format='csr'), np.concatenate(targets), value_scale, relative_floors


def fit_thermodynamic_constraints(low, grid_u, grid_tau, response_weight, max_unknowns=4500, multiplicity=2):
    started = time.perf_counter()
    unknowns = (len(knots_for(grid_u, multiplicity)) - 6) * (len(knots_for(grid_tau, multiplicity)) - 6)
    if unknowns > max_unknowns:
        raise ValueError(f'QR bound: {unknowns} unknowns exceed {max_unknowns}')
    matrix, target, value_scale, floors = weighted_system(low, grid_u, grid_tau, response_weight, multiplicity)
    original_target = target.copy()
    column_scale = np.maximum(np.sqrt(np.asarray(matrix.power(2).sum(axis=0)).ravel()), 1)
    equilibrated = np.asfortranarray(matrix.multiply(1 / column_scale).toarray())
    dense_bytes = int(equilibrated.nbytes)
    solve_started = time.perf_counter()
    solution, _, rank, _ = lstsq(equilibrated, target, cond=1e-13, lapack_driver='gelsy',
                                overwrite_a=True, overwrite_b=True, check_finite=False)
    solve_seconds = time.perf_counter() - solve_started
    diagnostics = {'method': 'dense rank-revealing gelsy QR, cond=1e-13, column equilibration; overwrite Fortran-contiguous dense work array',
        'rank': int(rank), 'unknowns': unknowns, 'constraints': int(matrix.shape[0]), 'dense_matrix_bytes': dense_bytes,
        'nullity': int(unknowns - rank),
        'max_unknowns': max_unknowns, 'response_weight': response_weight, 'relative_scale_floors': floors,
        'solve_seconds': solve_seconds, 'build_and_solve_seconds': time.perf_counter() - started,
        'weighting': 'Same experiment-015 linear relative E,TS,P/rho,P_rho,P/rho*chiT,T*Cv; response weight only on last three; no extra regularizer'}
    unscaled = solution / column_scale
    diagnostics['weighted_residual_norm'] = float(np.linalg.norm(matrix @ unscaled - original_target))
    if rank != unknowns:
        raise RankDeficientFit(diagnostics)
    shape = (len(knots_for(grid_u, multiplicity)) - 6, len(knots_for(grid_tau, multiplicity)) - 6)
    return unscaled.reshape(shape, order='F') * value_scale, diagnostics


def polynomial_coefficients(polynomial, grid_u, grid_tau):
    transforms = []
    for grid, powers in [(grid_u, polynomial.shape[0]), (grid_tau, polynomial.shape[1])]:
        matrix = np.concatenate([spline_basis(grid, grid, order) for order in range(3)])
        columns = []
        for power in range(powers):
            monomial = np.zeros(power + 1)
            monomial[power] = 1
            targets = np.concatenate([np.polynomial.polynomial.polyval(grid,
                np.polynomial.polynomial.polyder(monomial, order)) for order in range(3)])
            solution, _, rank, _ = lstsq(matrix, targets, cond=1e-13, lapack_driver='gelsy')
            if rank != matrix.shape[1]:
                raise ValueError('One-dimensional polynomial-representation system deficient')
            columns.append(solution)
        transforms.append(np.stack(columns, axis=1))
    return transforms[0] @ polynomial @ transforms[1].T


def null_witness(grid_u, grid_tau):
    directions, dimensions = [], []
    for grid in [grid_u, grid_tau]:
        spacing = np.median(np.diff(grid))
        matrix = np.concatenate([spline_basis(grid, grid), spacing * spline_basis(grid, grid, 1)])
        kernel = null_space(matrix, rcond=1e-13)
        if kernel.shape[1] < 1:
            raise ValueError('No value/first-derivative null direction')
        directions.append(kernel[:, 0])
        dimensions.append(int(kernel.shape[1]))
    return np.outer(*directions), dimensions
