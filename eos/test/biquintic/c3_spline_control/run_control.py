"""Bounded double-knot feasibility audit on archived experiment-015 inputs and queries."""

import argparse
import hashlib
import json
import os
import resource
import sys
import time
from datetime import datetime
from importlib.metadata import version
from pathlib import Path

for variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[variable] = '16'

import numpy as np
from threadpoolctl import threadpool_info, threadpool_limits

from spline_control import (
    BIQUINTIC, HERE, LOW_ORDERS, RankDeficientFit, evaluate_fit,
    fit_thermodynamic_constraints, knots_for, null_witness, spline_basis,
)
from freeeos_potential import CRAD, LN10, read_table, source_constraints, thermodynamic_partials

WORKSPACE = BIQUINTIC.parents[3]
OUTPUT = WORKSPACE / 'freeeos_biquintic_progress/experiments/029_c3_spline'
ARCHIVED = WORKSPACE / 'freeeos_biquintic_progress/experiments/015_smooth_validation'
TABLE = WORKSPACE / 'mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data'
TABLE_SHA = 'e6ca7ce27bee8197bb257a8c4886a06366e6951e58958fbed3735c877ee68502'
WEIGHTS = (.01, .1)
DEADLINE_HOUR = 20
FIELDS = ('P', 'Pgas', 'E', 'S', 'Cv', 'Cp', 'chiRho', 'chiT')
PARTIALS = tuple(field + suffix for field in FIELDS for suffix in ('_x', '_t'))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def summary(values):
    values = np.asarray(values)
    finite = np.isfinite(values)
    result = {'count': int(values.size), 'nonfinite': int((~finite).sum())}
    if finite.any():
        result.update(zip(('median', 'p99', 'max'), map(float, np.quantile(abs(values[finite]), [.5, .99, 1]))))
    return result


def positivity(state):
    return {field: {'nonfinite': int((~np.isfinite(state[field])).sum()),
                    'nonpositive_finite': int((np.isfinite(state[field]) & (state[field] <= 0)).sum())}
            for field in FIELDS}


def geometry(grid_q, grid_t):
    dimensions = [len(knots_for(grid_q)) - 6, len(knots_for(grid_t)) - 6]
    unknowns = int(np.prod(dimensions))
    nullity = 2 * (len(grid_q) + len(grid_t)) + 4
    return {'native_shape': [len(grid_q), len(grid_t)], 'basis_shape': dimensions, 'unknowns': unknowns,
            'six_constraint_equations': 6 * len(grid_q) * len(grid_t),
            'predicted_structural_nullity': nullity, 'predicted_rank': unknowns - nullity}


def continuity(coefficients, grid_q, grid_t, multiplicity=2, physical_fields=True):
    grid_u, grid_tau = grid_q * LN10, grid_t * LN10
    result = {}
    for axis, grid in enumerate((grid_u, grid_tau)):
        face = grid[len(grid) // 2]
        query_u, query_tau = np.array([(grid_u[0] + grid_u[-1]) / 2]), np.array([(grid_tau[0] + grid_tau[-1]) / 2])
        traces, states = [], []
        basis_jumps = {}
        for order in range(5):
            left = spline_basis(grid, [np.nextafter(face, -np.inf)], order, multiplicity)
            right = spline_basis(grid, [np.nextafter(face, np.inf)], order, multiplicity)
            scale = np.median(np.diff(grid))**order
            basis_jumps[str(order)] = float(np.max(abs(left - right)) * scale)
        for direction in (-np.inf, np.inf):
            sample_u, sample_tau = query_u.copy(), query_tau.copy()
            (sample_u if axis == 0 else sample_tau)[0] = np.nextafter(face, direction)
            native = evaluate_fit(coefficients, grid_u, grid_tau, sample_u, sample_tau, multiplicity=multiplicity)
            traces.append(native)
            if physical_fields:
                states.append(thermodynamic_partials(native, sample_u / LN10, sample_tau / LN10))
        native_jumps = {str(orders): summary((traces[1][orders] - traces[0][orders]) *
            np.median(np.diff(grid_u))**orders[0] * np.median(np.diff(grid_tau))**orders[1])
            for orders in LOW_ORDERS + ((3, 0), (2, 1), (1, 2), (0, 3))}
        result[f'axis{axis}'] = {'face_native_coordinate': float(face),
            'method': 'nextafter(face,-inf)/nextafter(face,+inf), floating one-sided limits; not an exact rational identity',
            'spacing_scaled_basis_derivative_jumps': basis_jumps, 'spacing_scaled_native_jet_jumps': native_jumps}
        if physical_fields:
            result[f'axis{axis}']['physical_value_and_response_partial_jumps'] = {field: summary((states[1][field] - states[0][field]) /
                np.maximum(1., np.maximum(abs(states[0][field]), abs(states[1][field])))) for field in FIELDS + PARTIALS}
        else:
            result[f'axis{axis}']['physical_response_status'] = 'Not evaluated: homogeneous null witness is not a thermodynamic EOS and may have singular/zero pressure denominators.'
    return result


def null_mode_audit(grid_q, grid_t, withheld_q, withheld_t, sample_q, sample_t):
    grid_u, grid_tau = grid_q * LN10, grid_t * LN10
    witness, dimensions = null_witness(grid_u, grid_tau)
    node_u, node_tau = np.meshgrid(grid_u, grid_tau, indexing='ij')
    nodal = evaluate_fit(witness, grid_u, grid_tau, node_u.ravel(), node_tau.ravel())
    spacing_u, spacing_tau = np.diff(grid_u)[0], np.diff(grid_tau)[0]
    constraints = {str(orders): summary(nodal[orders] * spacing_u**orders[0] * spacing_tau**orders[1]) for orders in LOW_ORDERS}
    withheld = evaluate_fit(witness, grid_u, grid_tau, withheld_q * LN10, withheld_t * LN10)
    offnode = evaluate_fit(witness, grid_u, grid_tau, sample_q[:256] * LN10, sample_t[:256] * LN10)
    physical_temperature = withheld[0, 1] - 2 * withheld[1, 0]
    changes = {'E': withheld[0, 0] - physical_temperature, 'TS': -physical_temperature,
               'P_over_rho': withheld[1, 0],
               'TCv': physical_temperature - withheld[0, 2] + 4 * withheld[1, 1] - 4 * withheld[2, 0]}
    result = {'claim': 'Nonzero homogeneous tensor mode: coefficients can vary without changing any six retained-node constraints. NOT a fitted EOS or physical field prediction.',
        'construction': 'Outer product of one-dimensional modes that vanish in value and first derivative at every retained grid node.',
        'one_dimensional_value_first_kernel_dimensions': dimensions,
        'coefficient_frobenius_norm': float(np.linalg.norm(witness)), 'scaled_six_retained_node_residuals': constraints,
        'same_archived_withheld_homogeneous_direction': {field: summary(value) for field, value in changes.items()},
        'first256_archived_offnode_homogeneous_value': summary(offnode[0, 0]),
        'C3_continuity_of_homogeneous_witness': continuity(witness, grid_q, grid_t, physical_fields=False)}
    return witness, result, changes


def compare_solved(coefficients, grid_q, grid_t, query_q, query_t, reference):
    native = evaluate_fit(coefficients, grid_q * LN10, grid_t * LN10, query_q * LN10, query_t * LN10)
    state = thermodynamic_partials(native, query_q, query_t)
    return state, {'relative_departures': {field: summary((state[field] - reference[field]) / abs(reference[field])) for field in FIELDS},
                   'positivity': positivity(state)}


def run(args):
    started = time.perf_counter()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    if (OUTPUT / 'results.json').exists() and not args.overwrite:
        raise ValueError('Products exist; use --overwrite')
    if digest(TABLE) != TABLE_SHA:
        raise ValueError('Unexpected pinned source-table SHA-256')
    archived_result = json.loads((ARCHIVED / 'results.json').read_text())
    _, logq, logt, rows = read_table(TABLE)
    low, valid, _ = source_constraints(rows, logq, logt)
    snapshots = {}
    for path, expected in archived_result['provenance']['files'].items():
        current = digest(Path(path))
        snapshots[path] = {'archived_sha256': expected, 'current_sha256': current, 'matches_archived': current == expected}
    if not snapshots[str(TABLE)]['matches_archived']:
        raise ValueError('015 table fingerprint mismatch')
    results = {'claim': 'Patch-local double-interior-knot C3 tensor quintic feasibility; archived table fidelity, not fresh-source accuracy or selected runtime.',
        'table_sha256': TABLE_SHA, 'archived_015_results_sha256': digest(ARCHIVED / 'results.json'),
        'archived_source_snapshot_comparison': snapshots,
        'source_hashes': {str(path.relative_to(WORKSPACE)): digest(path) for path in [HERE / 'spline_control.py', HERE / 'run_control.py', HERE / 'test_spline_control.py']},
        'versions': {'python': sys.version, 'numpy': np.__version__, 'scipy': version('scipy')},
        'command': sys.argv, 'timestamp': datetime.now().astimezone().isoformat(),
        'bounds': {'max_unknowns': args.max_unknowns, 'threads': 16, 'stop_new_QR_at_local_hour': DEADLINE_HOUR,
                   'max_new_fits': 16, 'response_weights': list(WEIGHTS), 'execution_order': 'all eight coarse fits before any full fits'},
        'C4_baseline': 'Reuse saved single-knot coefficients and saved value/query arrays from experiment015; no C4 refits.',
        'radiation_gauge': 'Same total-potential source constraints and current-CRAD total pressure convention as 015; no gauge or radiation refit.',
        'source_floor_default': {'current_default': source_constraints.__defaults__,
            'snapshot_note': 'Archived core hash may differ after the added optional minimum_logT parameter. This runner retains the unchanged default 3.2 floor and records both hashes, rather than asserting identical source files.'},
        'rank_policy': 'Do not export minimum-norm deficient QR coefficients as a selected fit. Rank failures remain failures. No unapproved null-space penalty, extra constraints or spline selector.',
        'patches': {}}
    retained = {}
    new_fits = 0
    for name, patch in archived_result['patches'].items():
        archive_path = ARCHIVED / name / 'fits_and_queries.npz'
        archive = np.load(archive_path)
        start_q, start_t, count_q, count_t = patch['selected_indices']
        patch_q, patch_t = logq[start_q:start_q + count_q], logt[start_t:start_t + count_t]
        np.testing.assert_array_equal(patch_q, archive['patch_logQ'])
        np.testing.assert_array_equal(patch_t, archive['patch_logT'])
        if not valid[start_q:start_q + count_q, start_t:start_t + count_t].all():
            raise ValueError(f'{name} original patch is no longer all native valid')
        patch_low = low[:, :, start_q:start_q + count_q, start_t:start_t + count_t]
        coarse_q, coarse_t, coarse_low = patch_q[::2], patch_t[::2], patch_low[:, :, ::2, ::2]
        density, temperature = np.meshgrid(np.arange(1, count_q, 2), np.arange(1, count_t, 2), indexing='ij')
        withheld_q, withheld_t = patch_q[density.ravel()], patch_t[temperature.ravel()]
        np.testing.assert_array_equal(withheld_q, archive['withheld_logQ'])
        np.testing.assert_array_equal(withheld_t, archive['withheld_logT'])
        source = rows[start_q + density.ravel(), start_t + temperature.ravel()]
        pgas = 10.0**source[:, 1]
        reference = {'Pgas': pgas, 'P': pgas + CRAD * (10.0**withheld_t)**4 / 3,
                     'E': 10.0**source[:, 2], 'S': 10.0**source[:, 3], 'Cv': source[:, 7],
                     'Cp': source[:, 6], 'chiRho': source[:, 4], 'chiT': source[:, 5]}
        metrics = archive['metrics'].tolist()
        withheld_metrics = metrics[:11]
        patch_result = {'actual_bounds_from015': patch['actual_bounds'], 'selected_indices_from015': patch['selected_indices'],
            'archived_query_sha256': digest(archive_path), 'withheld_count': len(withheld_q), 'offnode_count': len(archive['logQ']),
            'source_nodes': len(patch_q) * len(patch_t), 'full_geometry': geometry(patch_q, patch_t),
            'coarse_geometry': geometry(coarse_q, coarse_t), 'C4_archived': {}, 'C3_double': {}}
        products = {'patch_logQ': patch_q, 'patch_logT': patch_t, 'withheld_logQ': withheld_q, 'withheld_logT': withheld_t,
            'withheld_original_iq': start_q + density.ravel(), 'withheld_original_it': start_t + temperature.ravel(),
            'logQ': archive['logQ'], 'logT': archive['logT'], 'metrics': archive['metrics'],
            'biquintic_values': archive['biquintic_values'], 'local_coarse_common': archive['local_coarse_common']}
        for field in FIELDS:
            products['printed_reference_' + field] = reference[field]
        for weight in WEIGHTS:
            label = f'thermo_{weight:g}'
            archived_candidate = patch['candidates'][label]
            archived_coarse = archive[label + '_withheld_values']
            archived_state = {field: archived_coarse[withheld_metrics.index(field)] for field in FIELDS}
            patch_result['C4_archived'][label] = {'full': archived_candidate['full'], 'coarse': archived_candidate['coarse'],
                'coarse_same_query_relative_errors': {field: summary((archived_state[field] - reference[field]) / abs(reference[field])) for field in FIELDS}}
            for suffix in ['_values', '_withheld_values', '_full_coefficients', '_coarse_coefficients']:
                products['archived_C4_' + label + suffix] = archive[label + suffix]
            grid_q, grid_t = coarse_q, coarse_t
            if datetime.now().hour >= DEADLINE_HOUR:
                diagnostics = {'status': 'deferred_deadline', 'reason': 'No new QR after local 20:00'}
            else:
                try:
                    with threadpool_limits(limits=16):
                        coefficients, solve = fit_thermodynamic_constraints(coarse_low, grid_q * LN10, grid_t * LN10,
                                                                         weight, args.max_unknowns)
                    new_fits += 1
                    state, comparison = compare_solved(coefficients, grid_q, grid_t, withheld_q, withheld_t, reference)
                    diagnostics = {'status': 'solved', 'solve': solve, 'heldout_NOT_source_accuracy': comparison,
                                   'continuity': continuity(coefficients, grid_q, grid_t)}
                    products[label + '_C3_coarse_coefficients'] = coefficients
                    products[label + '_C3_withheld_values'] = np.stack([state[field] for field in FIELDS])
                except RankDeficientFit as failure:
                    new_fits += 1
                    diagnostics = {'status': 'rank_deficient_no_selected_fit', 'solve': failure.diagnostics,
                                   'reason': 'Six nodal constraints do not determine the double-knot tensor coefficients; no deficient solution exported as EOS prediction.'}
                except (ValueError, np.linalg.LinAlgError) as error:
                    diagnostics = {'status': 'failed', 'reason': str(error)}
            patch_result['C3_double'][label] = {'coarse': diagnostics}
            if diagnostics['status'] == 'rank_deficient_no_selected_fit':
                if diagnostics['solve']['rank'] != patch_result['coarse_geometry']['predicted_rank']:
                    raise AssertionError('Observed QR rank differs from structural rank prediction')
            print(f'{name} {label} coarse: {diagnostics["status"]}; {diagnostics.get("solve", {})}', flush=True)
        witness, null_audit, changes = null_mode_audit(coarse_q, coarse_t, withheld_q, withheld_t, archive['logQ'], archive['logT'])
        patch_result['coarse_homogeneous_null_mode'] = null_audit
        products['homogeneous_null_coefficients'] = witness
        for field, value in changes.items():
            products['homogeneous_withheld_direction_' + field] = value
        retained[name] = (patch_q, patch_t, patch_low, products, reference, patch_result)
        results['patches'][name] = patch_result
    all_coarse_deficient = all(candidate['coarse']['status'] == 'rank_deficient_no_selected_fit'
        for patch in results['patches'].values() for candidate in patch['C3_double'].values())
    for name, (patch_q, patch_t, patch_low, products, reference, patch_result) in retained.items():
        for weight in WEIGHTS:
            label = f'thermo_{weight:g}'
            if all_coarse_deficient:
                diagnostics = {'status': 'deferred_structural_rank_failure',
                    'reason': 'All eight coarse fits deficient and an analytic nonzero six-constraint null mode exists on every full grid. No expensive full QR or arbitrary selector justified.'}
            elif args.coarse_only or datetime.now().hour >= DEADLINE_HOUR:
                diagnostics = {'status': 'deferred_bound', 'reason': 'Coarse-only flag or 20:00 new-QR deadline'}
            else:
                try:
                    with threadpool_limits(limits=16):
                        coefficients, solve = fit_thermodynamic_constraints(patch_low, patch_q * LN10, patch_t * LN10,
                                                                         weight, args.max_unknowns)
                    new_fits += 1
                    diagnostics = {'status': 'solved', 'solve': solve, 'continuity': continuity(coefficients, patch_q, patch_t)}
                    products[label + '_C3_full_coefficients'] = coefficients
                except RankDeficientFit as failure:
                    new_fits += 1
                    diagnostics = {'status': 'rank_deficient_no_selected_fit', 'solve': failure.diagnostics}
                except (ValueError, np.linalg.LinAlgError) as error:
                    diagnostics = {'status': 'failed', 'reason': str(error)}
            patch_result['C3_double'][label]['full'] = diagnostics
        patch_output = OUTPUT / name
        patch_output.mkdir(parents=True, exist_ok=True)
        (patch_output / 'results.json').write_text(json.dumps(patch_result, indent=2, allow_nan=False) + '\n')
        np.savez_compressed(patch_output / 'inputs_queries_and_null_audit.npz', **products)
    results['new_QR_fits_attempted'] = new_fits
    results['structural_interpretation'] = {
        'status': 'rank_failure' if all_coarse_deficient else 'consult_individual_status',
        'dimension': 'For Nq,Nt grid nodes, double-knot basis dimensions are (2Nq+2)*(2Nt+2). The six-node-jet map lacks the (1,2),(2,1),(2,2) tensor flag blocks.',
        'predicted_nullity': '2*Nq+2*Nt+4; predicted rank=4*Nq*Nt+2*Nq+2*Nt. Checked against every attempted patch QR.',
        'witness': 'Take nonzero one-dimensional splines with value and first derivative zero at all nodes; their tensor product has all six low nodal jets zero but changes offnode/withheld values.',
        'hard_E_TS_P': 'Enforcing E,TS,P/rho exactly and fitting second responses softly retains these same homogeneous null directions; constrained KKT/null-space solve alone cannot choose a unique potential.',
        'missing_selection_rule': 'Would require an explicit additional high-jet constraint, penalty or minimum-norm/gauge-dependent selector and separate validation; not silently added here.',
        'limits': 'More DOFs and C3 continuity do not imply better printed-field fidelity or truth. No nonunique QR vector is reported as a selected interpolator.'}
    results['threadpools'] = threadpool_info()
    results['elapsed_seconds'] = time.perf_counter() - started
    results['peak_process_rss_kib_linux'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    (OUTPUT / 'results.json').write_text(json.dumps(results, indent=2, allow_nan=False) + '\n')
    lines = ['# Double-knot C3 spline feasibility', '', results['claim'], '',
        'Same four experiment-015 patches, archived query IDs, native constraints, weights .01/.1; no C4 refits.', '',
        '| Patch | Coarse unknowns | Equations | QR rank | Nullity | C4 Cv p99 (.01 / .1) |',
        '| --- | ---: | ---: | ---: | ---: | --- |']
    for name, patch in results['patches'].items():
        solve = patch['C3_double']['thermo_0.1']['coarse'].get('solve', {})
        errors = [patch['C4_archived'][f'thermo_{weight:g}']['coarse_same_query_relative_errors']['Cv']['p99'] for weight in WEIGHTS]
        lines.append(f"| {name} | {patch['coarse_geometry']['unknowns']} | {patch['coarse_geometry']['six_constraint_equations']} | {solve.get('rank', 'n/a')} | {solve.get('nullity', 'n/a')} | {errors[0]:.3e} / {errors[1]:.3e} |")
    lines.extend(['', 'All eight attempted double-knot coarse QR fits are rank deficient; full fits are deferred structurally, not timed out.' if all_coarse_deficient else 'Consult per-fit status.',
        'No C3 heldout/positivity superiority claim or deficient minimum-norm coefficient output. Homogeneous null-mode diagnostics are explicitly not EOS predictions.',
        'Double knots give C3, not C4: partition/polynomial representation and one-sided third/physical-response continuity tests are retained.',
        'Hard E/TS/P enforcement plus soft second constraints alone does not eliminate the demonstrated null modes. Additional selector deferred.',
        f"New QR fits attempted: {new_fits}; elapsed {results['elapsed_seconds']:.2f}s; Linux peak RSS {results['peak_process_rss_kib_linux']} KiB."])
    (OUTPUT / 'README.md').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-unknowns', type=int, default=4500)
    parser.add_argument('--coarse-only', action='store_true')
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    if not 1 <= args.max_unknowns <= 4500:
        parser.error('Require 1..4500 maximum unknowns')
    run(args)


if __name__ == '__main__':
    main()
