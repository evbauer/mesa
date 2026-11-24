"""Explicit integrated-normal-fourth-derivative selector for one archived C3 patch."""

import argparse
import hashlib
import json
import os
import resource
import signal
import sys
import time
from datetime import datetime
from importlib.metadata import version
from pathlib import Path

for variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[variable] = '16'

import numpy as np
from scipy.linalg import lstsq
from scipy.interpolate import BSpline
from scipy.sparse import csr_matrix, kron, vstack
from scipy.sparse.linalg import lsmr
from threadpoolctl import threadpool_info, threadpool_limits

from spline_control import (
    BIQUINTIC, HERE, evaluate_fit, knots_for, null_witness, polynomial_coefficients,
    spline_basis, thermodynamic_combinations, weighted_system,
)
from freeeos_potential import CRAD, LN10, LOW_ORDERS, read_table, source_constraints, thermodynamic_partials

WORKSPACE = BIQUINTIC.parents[3]
OUTPUT = WORKSPACE / 'freeeos_biquintic_progress/experiments/031_regularized_c3'
ARCHIVED = WORKSPACE / 'freeeos_biquintic_progress/experiments/015_smooth_validation'
TABLE = WORKSPACE / 'mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data'
TABLE_SHA = 'e6ca7ce27bee8197bb257a8c4886a06366e6951e58958fbed3735c877ee68502'
LAMBDAS = (.001, .01, .1)
RESPONSE_WEIGHT = .1
FIELDS = ('P', 'Pgas', 'E', 'S', 'Cv', 'Cp', 'chiRho', 'chiT')
PARTIALS = tuple(field + suffix for field in FIELDS for suffix in ('_x', '_t'))
MAX_DENSE_BYTES = 128 * 1024**2
MAX_FIT_SECONDS = 90


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def gauss_grid(grid, count):
    positions, weights = np.polynomial.legendre.leggauss(count)
    half_width = np.diff(grid) / 2
    centers = (grid[:-1] + grid[1:]) / 2
    return (centers[:, None] + half_width[:, None] * positions).ravel(), (half_width[:, None] * weights).ravel()


def roughness_operator(grid_u, grid_tau):
    factors = []
    for grid in [grid_u, grid_tau]:
        fourth_query, fourth_weight = gauss_grid(grid, 2)
        mass_query, mass_weight = gauss_grid(grid, 6)
        fourth = csr_matrix(spline_basis(grid, fourth_query, 4) * np.sqrt(fourth_weight)[:, None])
        mass = csr_matrix(spline_basis(grid, mass_query) * np.sqrt(mass_weight)[:, None])
        factors.append((fourth, mass))
    spacing_u, spacing_tau = np.median(np.diff(grid_u)), np.median(np.diff(grid_tau))
    area = (grid_u[-1] - grid_u[0]) * (grid_tau[-1] - grid_tau[0])
    density = kron(factors[1][1], factors[0][0], format='csr') * spacing_u**4 / np.sqrt(area)
    temperature = kron(factors[1][0], factors[0][1], format='csr') * spacing_tau**4 / np.sqrt(area)
    return vstack([density, temperature], format='csr'), {
        'formula': '(h_u^8*integral(f_uuuu^2)+h_tau^8*integral(f_tautautautau^2))/(patch_area*median_abs_f^2)',
        'quadrature': '2 Gauss points per interval in fourth-derivative axis; 6 in value/mass axis. Exact piecewise polynomial integration (degrees 2 and 10).',
        'native_spacings': [float(spacing_u), float(spacing_tau)], 'native_patch_area': float(area),
        'rows': int(density.shape[0] + temperature.shape[0]), 'columns': int(density.shape[1]),
        'sparse_bytes': int(sum(array.nbytes for array in (density.data, density.indices, density.indptr,
                                                          temperature.data, temperature.indices, temperature.indptr)))}


def summary(values):
    values = np.asarray(values)
    finite = np.isfinite(values)
    result = {'count': int(values.size), 'nonfinite': int((~finite).sum())}
    if finite.any():
        result.update(zip(('median', 'p99', 'max'), map(float, np.quantile(abs(values[finite]), [.5, .99, 1]))))
    return result


def positive_status(state):
    return {field: {'nonfinite': int((~np.isfinite(state[field])).sum()),
                    'nonpositive_finite': int((np.isfinite(state[field]) & (state[field] <= 0)).sum())}
            for field in FIELDS}


def timeout_handler(signum, frame):
    raise TimeoutError('Per-fit 90-second bound exceeded')


def insert_double_knots(coefficients, grid_u, grid_tau):
    first = BSpline(knots_for(grid_u, 1), coefficients, 5, axis=0)
    for knot in grid_u[1:-1]:
        first = first.insert_knot(knot)
    second = BSpline(knots_for(grid_tau, 1), first.c.T, 5, axis=0)
    for knot in grid_tau[1:-1]:
        second = second.insert_knot(knot)
    return second.c.T


def solve_regularized(low, grid_u, grid_tau, regularization, response_weight=RESPONSE_WEIGHT,
                      force_dense=False, verify_sparse=False, initial_coefficients=None):
    started = time.perf_counter()
    data, target, value_scale, floors = weighted_system(low, grid_u, grid_tau, response_weight)
    source_rows = data.shape[0]
    data = data / np.sqrt(source_rows)
    target = target / np.sqrt(source_rows)
    rough, rough_metadata = roughness_operator(grid_u, grid_tau)
    matrix = vstack([data, np.sqrt(regularization) * rough], format='csr')
    rhs = np.r_[target, np.zeros(rough.shape[0])]
    scale = np.maximum(np.sqrt(np.asarray(matrix.power(2).sum(axis=0)).ravel()), 1e-30)
    equilibrated = matrix.multiply(1 / scale).tocsr()
    estimated_dense = int(np.prod(equilibrated.shape)) * 8
    if matrix.shape[1] > 4500:
        raise ValueError('Regularized unknown bound >4500')
    diagnostics = {'lambda': regularization, 'response_weight': response_weight, 'unknowns': int(matrix.shape[1]),
        'source_rows': source_rows, 'augmented_rows': int(matrix.shape[0]), 'estimated_dense_bytes': estimated_dense,
        'dense_byte_bound': MAX_DENSE_BYTES, 'roughness': rough_metadata, 'potential_scale': value_scale,
        'source_objective': 'Mean (over six combinations/all source nodes) squared weighted relative residual; response importance .1 enters squared.',
        'selector': 'Explicit J = mean weighted-relative source residual^2 + lambda*normalized integrated normal-fourth roughness; sqrt(lambda) multiplies operator rows.',
        'relative_source_floors': floors}
    old_handler = signal.signal(signal.SIGALRM, timeout_handler)
    signal.setitimer(signal.ITIMER_REAL, MAX_FIT_SECONDS)
    try:
        with threadpool_limits(limits=16):
            if estimated_dense <= MAX_DENSE_BYTES or force_dense:
                dense = np.asfortranarray(equilibrated.toarray())
                solution, _, rank, _ = lstsq(dense, rhs.copy(), cond=1e-13, lapack_driver='gelsy',
                                            overwrite_a=True, overwrite_b=True, check_finite=False)
                diagnostics.update({'method': 'Dense rank-revealing gelsy QR after column equilibration', 'rank': int(rank)})
                if rank != matrix.shape[1]:
                    raise RuntimeError(f'Regularized QR still deficient: {rank}/{matrix.shape[1]}')
                if verify_sparse:
                    sparse = lsmr(equilibrated, rhs, atol=1e-11, btol=1e-11, conlim=1e12, maxiter=20000)
                    diagnostics['sparse_crosscheck'] = {'stop': int(sparse[1]), 'iterations': int(sparse[2]),
                        'relative_scaled_solution_difference': float(np.linalg.norm(sparse[0] - solution) / np.linalg.norm(solution))}
            else:
                if initial_coefficients is None:
                    raise ValueError('Full sparse fit requires an explicit initial coefficient representation')
                initial = initial_coefficients.ravel(order='F') / value_scale * scale
                shifted_rhs = rhs - equilibrated @ initial
                sparse = lsmr(equilibrated, shifted_rhs, atol=1e-12, btol=1e-12, conlim=1e12, maxiter=40000)
                solution = initial + sparse[0]
                diagnostics.update({'method': 'Sparse LSMR after column equilibration; full dense QR exceeds 128MiB work-array bound',
                    'rank': None, 'rank_evidence': 'Analytic: roughness kernel consists of tensor bicubic polynomials, eliminated by nodal value constraints on >=4 nodes per axis; coarse augmented QR checks full rank numerically.',
                    'lsmr_stop': int(sparse[1]), 'lsmr_iterations': int(sparse[2]), 'lsmr_residual_norm': float(sparse[3]),
                    'lsmr_normal_residual_norm': float(sparse[4]), 'lsmr_condition_estimate': float(sparse[6]),
                    'initialization': 'Exact C4-to-C3 knot-insertion representation; solve for coefficient correction using shifted residual. Initialization is not a selector: augmented objective is full rank.',
                    'shifted_rhs_norm': float(np.linalg.norm(shifted_rhs))})
                if sparse[1] not in (2, 5):
                    raise RuntimeError(f'LSMR did not converge, stop={sparse[1]}, iterations={sparse[2]}')
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old_handler)
    normalized = solution / scale
    residual = equilibrated @ solution - rhs
    normal_residual = equilibrated.T @ residual
    diagnostics.update({'elapsed_seconds': time.perf_counter() - started,
        'data_objective': float(np.sum((data @ normalized - target)**2)),
        'roughness_objective': float(np.sum((rough @ normalized)**2)),
        'augmented_stationarity_infinity': float(np.max(abs(normal_residual))),
        'augmented_stationarity_relative': float(np.linalg.norm(normal_residual) /
                                               max(1., np.linalg.norm(equilibrated.T @ rhs)))})
    diagnostics['total_objective'] = diagnostics['data_objective'] + regularization * diagnostics['roughness_objective']
    if diagnostics['augmented_stationarity_relative'] > 1e-10:
        raise RuntimeError(f'Augmented stationarity guard failed: {diagnostics["augmented_stationarity_relative"]}')
    shape = (len(knots_for(grid_u)) - 6, len(knots_for(grid_tau)) - 6)
    return normalized.reshape(shape, order='F') * value_scale, diagnostics


def state_at(coefficients, grid_q, grid_t, sample_q, sample_t):
    native = evaluate_fit(coefficients, grid_q * LN10, grid_t * LN10, sample_q * LN10, sample_t * LN10)
    return thermodynamic_partials(native, sample_q, sample_t), native


def errors(state, reference):
    return {field: summary((state[field] - reference[field]) / abs(reference[field])) for field in FIELDS}


def continuity(coefficients, grid_q, grid_t):
    result = {}
    for axis, grid in enumerate([grid_q * LN10, grid_t * LN10]):
        face = grid[len(grid) // 2]
        sample_u = np.array([(grid_q[0] + grid_q[-1]) * LN10 / 2])
        sample_tau = np.array([(grid_t[0] + grid_t[-1]) * LN10 / 2])
        traces = []
        for direction in [-np.inf, np.inf]:
            query_u, query_tau = sample_u.copy(), sample_tau.copy()
            (query_u if axis == 0 else query_tau)[0] = np.nextafter(face, direction)
            traces.append(state_at(coefficients, grid_q, grid_t, query_u / LN10, query_tau / LN10)[0])
        result[f'axis{axis}'] = {'method': 'floating nextafter one-sided limits',
            'scaled_jumps': {field: summary((traces[1][field] - traces[0][field]) /
                np.maximum(1., np.maximum(abs(traces[0][field]), abs(traces[1][field])))) for field in FIELDS + PARTIALS},
            'left_positivity': positive_status(traces[0]), 'right_positivity': positive_status(traces[1])}
    return result


def run(args):
    started = time.perf_counter()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    if (OUTPUT / 'results.json').exists() and not args.overwrite:
        raise ValueError('Products exist; use --overwrite')
    previous_attempts = {}
    if (OUTPUT / 'results.json').exists():
        previous = json.loads((OUTPUT / 'results.json').read_text())
        previous_attempts = {resolution: {label: {key: candidate[key] for key in ['status', 'reason', 'solve'] if key in candidate}
            for label, candidate in fits.items()} for resolution, fits in previous.get('fits', {}).items()}
    if digest(TABLE) != TABLE_SHA:
        raise ValueError('Pinned source table changed')
    archived_result = json.loads((ARCHIVED / 'results.json').read_text())
    patch = archived_result['patches']['active_chiT']
    archive_path = ARCHIVED / 'active_chiT/fits_and_queries.npz'
    archive = np.load(archive_path)
    _, logq, logt, rows = read_table(TABLE)
    low, valid, _ = source_constraints(rows, logq, logt)
    start_q, start_t, count_q, count_t = patch['selected_indices']
    patch_q, patch_t = logq[start_q:start_q + count_q], logt[start_t:start_t + count_t]
    np.testing.assert_array_equal(patch_q, archive['patch_logQ'])
    np.testing.assert_array_equal(patch_t, archive['patch_logT'])
    if not valid[start_q:start_q + count_q, start_t:start_t + count_t].all():
        raise ValueError('Active patch no longer all native valid')
    patch_low = low[:, :, start_q:start_q + count_q, start_t:start_t + count_t]
    withheld_q, withheld_t = archive['withheld_logQ'], archive['withheld_logT']
    density, temperature = np.meshgrid(np.arange(1, count_q, 2), np.arange(1, count_t, 2), indexing='ij')
    reference_rows = rows[start_q + density.ravel(), start_t + temperature.ravel()]
    products = {'patch_logQ': patch_q, 'patch_logT': patch_t, 'withheld_logQ': withheld_q, 'withheld_logT': withheld_t,
        'withheld_iq': start_q + density.ravel(), 'withheld_it': start_t + temperature.ravel(),
        'offnode_logQ': archive['logQ'], 'offnode_logT': archive['logT'], 'metrics': archive['metrics'],
        'archived_biquintic_values': archive['biquintic_values']}
    reference = {'Pgas': 10.0**reference_rows[:, 1], 'E': 10.0**reference_rows[:, 2], 'S': 10.0**reference_rows[:, 3],
        'Cv': reference_rows[:, 7], 'Cp': reference_rows[:, 6], 'chiRho': reference_rows[:, 4], 'chiT': reference_rows[:, 5]}
    reference['P'] = reference['Pgas'] + CRAD * (10.0**withheld_t)**4 / 3
    for field in FIELDS:
        products['printed_reference_' + field] = reference[field]
    metrics = archive['metrics'].tolist()
    archived_coarse = {field: archive['thermo_0.1_withheld_values'][metrics[:11].index(field)] for field in FIELDS}
    products['archived_C4_thermo0.1_withheld_values'] = archive['thermo_0.1_withheld_values']
    products['archived_C4_thermo0.1_values'] = archive['thermo_0.1_values']
    results = {'claim': 'Explicit roughness-selected double-knot C3 potential on one archived active geometry patch; table-only fidelity, not fresh-source accuracy, active fraction or chosen runtime.',
        'table_sha256': TABLE_SHA, '015_archive_sha256': digest(archive_path), '015_results_sha256': digest(ARCHIVED / 'results.json'),
        '029_results_sha256': digest(WORKSPACE / 'freeeos_biquintic_progress/experiments/029_c3_spline/results.json'),
        'source_hashes': {str(path.relative_to(WORKSPACE)): digest(path) for path in [HERE / 'regularized.py', HERE / 'test_regularized.py',
            HERE / 'spline_control.py', BIQUINTIC / 'freeeos_potential.py']}, 'command': sys.argv,
        'timestamp': datetime.now().astimezone().isoformat(), 'versions': {'python': sys.version, 'numpy': np.__version__, 'scipy': version('scipy')},
        'patch': {'label': 'active_chiT', 'actual_bounds_from015': patch['actual_bounds'], 'selected_indices_from015': patch['selected_indices'],
                  'holdouts': len(withheld_q), 'offnodes': len(archive['logQ'])},
        'bound': {'threads': 16, 'per_fit_seconds': MAX_FIT_SECONDS, 'dense_QR_matrix_bytes': MAX_DENSE_BYTES,
                  'max_unknowns': 4500, 'max_fits': 6, 'no_new_solve_after_local_hour': 20},
        'objective': 'mean weighted relative six-combination source residual^2 + lambda*normalized integrated native normal-fourth roughness; data response importance .1 fixed; lambdas .001/.01/.1.',
        'previous_attempts': previous_attempts,
        'nullity_remedy': '029 nonuniqueness is removed by an explicit physical-coordinate-derivative roughness selector, not arbitrary coefficient minimum norm. Zero roughness is a tensor bicubic polynomial; source values eliminate that common homogeneous kernel.',
        'conventions': 'Total potential/current CRAD as015; no fitted entropy/energy gauge or source constant. Native coordinates u=ln10*logQ,tau=lnT; normal penalty is not a rotationally invariant or physical-prior optimum.',
        'C4_archived_no_refit': {'coarse_same_query_errors': errors(archived_coarse, reference),
                               'full': patch['candidates']['thermo_0.1']['full'], 'coarse': patch['candidates']['thermo_0.1']['coarse']},
        'fits': {}}
    for resolution in ['coarse', 'full']:
        grid_q, grid_t = (patch_q[::2], patch_t[::2]) if resolution == 'coarse' else (patch_q, patch_t)
        fit_low = patch_low[:, :, ::2, ::2] if resolution == 'coarse' else patch_low
        initial_coefficients = insert_double_knots(archive['thermo_0.1_full_coefficients'],
            patch_q * LN10, patch_t * LN10) if resolution == 'full' else None
        results['fits'][resolution] = {}
        for regularization in LAMBDAS:
            label = f'lambda_{regularization:g}'
            if datetime.now().hour >= 20 or (resolution == 'full' and args.coarse_only):
                results['fits'][resolution][label] = {'status': 'deferred_bound'}
                continue
            try:
                coefficients, solve = solve_regularized(fit_low, grid_q * LN10, grid_t * LN10, regularization,
                    verify_sparse=resolution == 'coarse' and regularization == .01, initial_coefficients=initial_coefficients)
                diagnostics = {'status': 'solved', 'solve': solve, 'continuity': continuity(coefficients, grid_q, grid_t)}
                node_q, node_t = np.meshgrid(grid_q, grid_t, indexing='ij')
                nodal, native = state_at(coefficients, grid_q, grid_t, node_q.ravel(), node_t.ravel())
                fitted_combinations, source_combinations = thermodynamic_combinations(native), thermodynamic_combinations(fit_low)
                diagnostics['six_combination_nodal_soft_reconciliation'] = {field: summary((fitted_combinations[field] - target.ravel()) /
                    np.maximum(abs(target.ravel()), .01 * np.median(abs(target)))) for field, target in source_combinations.items()}
                diagnostics['nodal_positivity'] = positive_status(nodal)
                if resolution == 'coarse':
                    state, _ = state_at(coefficients, grid_q, grid_t, withheld_q, withheld_t)
                    diagnostics['heldout_relative_errors_NOT_source_accuracy'] = errors(state, reference)
                    diagnostics['heldout_positivity'] = positive_status(state)
                    for field in FIELDS:
                        products[f'{resolution}_{label}_withheld_{field}'] = state[field]
                else:
                    state, _ = state_at(coefficients, grid_q, grid_t, archive['logQ'], archive['logT'])
                    bq = {field: archive['biquintic_values'][metrics.index(field)] for field in FIELDS}
                    diagnostics['offnode_relative_departures_fromBQ_NOT_accuracy'] = errors(state, bq)
                    diagnostics['offnode_response_partial_departures_fromBQ_NOT_accuracy'] = {}
                    for field in ['Cv', 'Cp', 'chiRho', 'chiT']:
                        for suffix in ['_x', '_t']:
                            partial = field + suffix
                            bq_partial = archive['biquintic_values'][metrics.index(partial)]
                            scale = np.maximum(1., np.maximum(abs(bq[field]), abs(bq_partial)))
                            diagnostics['offnode_response_partial_departures_fromBQ_NOT_accuracy'][partial] = summary((state[partial] - bq_partial) / scale)
                    diagnostics['offnode_positivity'] = positive_status(state)
                    for field in FIELDS + PARTIALS:
                        products[f'{resolution}_{label}_offnode_{field}'] = state[field]
                products[f'{resolution}_{label}_coefficients'] = coefficients
                results['fits'][resolution][label] = diagnostics
                print(f'{resolution} {label}: solved; {solve}', flush=True)
            except (ValueError, RuntimeError, TimeoutError, np.linalg.LinAlgError) as error:
                results['fits'][resolution][label] = {'status': 'failed', 'reason': str(error)}
                print(f'{resolution} {label}: FAILED {error}', flush=True)
            results['elapsed_seconds'] = time.perf_counter() - started
            results['peak_process_rss_kib_linux'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            (OUTPUT / 'results.json').write_text(json.dumps(results, indent=2, allow_nan=False) + '\n')
    np.savez_compressed(OUTPUT / 'paired_fits_and_queries.npz', **products)
    results['threadpools'] = threadpool_info()
    results['elapsed_seconds'] = time.perf_counter() - started
    results['peak_process_rss_kib_linux'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    (OUTPUT / 'results.json').write_text(json.dumps(results, indent=2, allow_nan=False) + '\n')
    lines = ['# Explicit-roughness regularized C3 control', '', results['claim'], '',
        'Active geometry only; source response importance .1; same 204 archived stride-2 holdouts.', '',
        '| Coarse candidate | Cv p99 / max | Pgas p99 / max | E p99 / max | S p99 / max |',
        '| --- | --- | --- | --- | --- |']
    for name, candidate in [('archived_C4_thermo0.1', {'heldout_relative_errors_NOT_source_accuracy': results['C4_archived_no_refit']['coarse_same_query_errors']})] + list(results['fits']['coarse'].items()):
        stats = candidate.get('heldout_relative_errors_NOT_source_accuracy')
        if stats:
            lines.append('| ' + name + ' | ' + ' | '.join(f"{stats[field]['p99']:.3e} / {stats[field]['max']:.3e}" for field in ['Cv', 'Pgas', 'E', 'S']) + ' |')
        else:
            lines.append(f'| {name} | {candidate["status"]} | - | - | - |')
    lines.extend(['', 'All selected coefficients use the stated derivative roughness penalty; no arbitrary deficient coefficient minimum-norm acceptance.',
        'Full-grid solver statuses, stationarity, all-field/response comparisons, nodal reconciliation and sampled positivity are in results.json.',
        'C3 does not imply C4, physical fidelity or global positivity. Lambda comparison is a declared selector sensitivity, not production parameter selection.'])
    (OUTPUT / 'README.md').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--coarse-only', action='store_true')
    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
