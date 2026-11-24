"""One-patch exact-first thermodynamic constraints with explicit C3 roughness selector."""

import os
import sys
import time
from pathlib import Path

for variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[variable] = '16'

import numpy as np
from scipy.linalg import lstsq, qr, solve_triangular
from scipy.sparse import vstack
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
BIQUINTIC = HERE.parent
sys.path.insert(0, str(BIQUINTIC / 'c3_spline_control'))
from regularized import roughness_operator
from spline_control import knots_for, weighted_system


class HardConstraintFailure(RuntimeError):
    def __init__(self, diagnostics):
        self.diagnostics = diagnostics
        super().__init__(diagnostics['reason'])


def solve_hard_system(hard, target, soft, rhs):
    columns = hard.shape[1]
    if hard.shape[0] > columns or columns > 4500:
        raise ValueError('Hard-row/unknown resource bound violated')
    scale = np.maximum(np.sqrt(np.sum(hard**2, axis=0) + np.asarray(soft.power(2).sum(axis=0)).ravel()), 1e-30)
    equilibrated_hard = hard / scale
    orthogonal, upper, pivot = qr(equilibrated_hard.T, mode='full', pivoting=True)
    diagonal = abs(np.diag(upper[:hard.shape[0]]))
    rank = int((diagonal > 1e-13 * diagonal.max()).sum())
    diagnostics = {'hard_rows': int(hard.shape[0]), 'unknowns': columns, 'hard_rank': rank,
                   'hard_row_rank_policy': 'Require full independent hard-row rank; no rank-deficient minimum-norm hard solution accepted'}
    if rank != hard.shape[0]:
        diagnostics['reason'] = f'Hard equality matrix deficient: {rank}/{hard.shape[0]}'
        raise HardConstraintFailure(diagnostics)
    independent = upper[:hard.shape[0], :hard.shape[0]]
    coordinates = solve_triangular(independent.T, target[pivot], lower=True)
    particular = orthogonal[:, :hard.shape[0]] @ coordinates
    kernel = orthogonal[:, hard.shape[0]:]
    scaled_soft = soft.multiply(1 / scale).tocsr()
    reduced = scaled_soft @ kernel
    reduced_rhs = rhs - scaled_soft @ particular
    work = np.asfortranarray(reduced)
    solution, _, reduced_rank, _ = lstsq(work, reduced_rhs, cond=1e-13, lapack_driver='gelsy',
                                         overwrite_a=True, overwrite_b=False, check_finite=False)
    diagnostics.update({'nullspace_columns': int(kernel.shape[1]), 'reduced_rank': int(reduced_rank),
        'reduced_dense_bytes': int(reduced.nbytes), 'orthogonal_dense_bytes': int(orthogonal.nbytes),
        'particular_hard_relative_residual_max': float(np.max(abs(equilibrated_hard @ particular - target)))})
    if reduced_rank != kernel.shape[1]:
        diagnostics['reason'] = f'Selector reduced system deficient: {reduced_rank}/{kernel.shape[1]}'
        raise HardConstraintFailure(diagnostics)
    scaled_solution = particular + kernel @ solution
    hard_residual = equilibrated_hard @ scaled_solution - target
    reduced_gradient = reduced.T @ (reduced @ solution - reduced_rhs)
    diagnostics.update({'hard_relative_residual_max': float(np.max(abs(hard_residual))),
        'reduced_stationarity_relative': float(np.linalg.norm(reduced_gradient) /
                                              max(1., np.linalg.norm(reduced.T @ reduced_rhs)))})
    if diagnostics['hard_relative_residual_max'] > 1e-9 or diagnostics['reduced_stationarity_relative'] > 1e-10:
        diagnostics['reason'] = 'Hard-constraint residual or reduced-stationarity guard failed'
        raise HardConstraintFailure(diagnostics)
    return scaled_solution / scale, diagnostics


def solve_constrained(low, grid_u, grid_tau, regularization, response_weight=.1):
    started = time.perf_counter()
    data, target, value_scale, floors = weighted_system(low, grid_u, grid_tau, response_weight)
    node_count = len(grid_u) * len(grid_tau)
    hard_count = 3 * node_count
    hard = data[:hard_count].toarray()
    hard_target = target[:hard_count]
    source_soft = data[hard_count:] / np.sqrt(data.shape[0])
    soft_target = target[hard_count:] / np.sqrt(data.shape[0])
    rough, rough_metadata = roughness_operator(grid_u, grid_tau)
    soft = vstack([source_soft, np.sqrt(regularization) * rough], format='csr')
    rhs = np.r_[soft_target, np.zeros(rough.shape[0])]
    with threadpool_limits(limits=16):
        normalized, diagnostics = solve_hard_system(hard, hard_target, soft, rhs)
    diagnostics.update({'lambda': regularization, 'response_weight': response_weight, 'elapsed_seconds': time.perf_counter() - started,
        'method': 'Column-scaled pivoted QR of hard equality transpose, orthogonal nullspace elimination, full-rank reduced gelsy QR',
        'source_objective_normalization': 'Sum soft response residual squared divided by six*source_nodes, matching031; hard first residuals zero.',
        'roughness': rough_metadata, 'potential_scale': value_scale, 'relative_source_floors': floors,
        'soft_data_objective': float(np.sum((source_soft @ normalized - soft_target)**2)),
        'roughness_objective': float(np.sum((rough @ normalized)**2))})
    diagnostics['total_objective'] = diagnostics['soft_data_objective'] + regularization * diagnostics['roughness_objective']
    shape = (len(knots_for(grid_u)) - 6, len(knots_for(grid_tau)) - 6)
    return normalized.reshape(shape, order='F') * value_scale, diagnostics
