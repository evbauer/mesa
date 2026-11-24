"""Bounded active-patch hard E/TS/P control versus frozen031 soft fits and archived015 C4."""

import argparse
import hashlib
import json
import resource
import signal
import sys
import time
from datetime import datetime
from importlib.metadata import version
from pathlib import Path

import constrained
import numpy as np
from constrained import BIQUINTIC, HERE, HardConstraintFailure, solve_constrained
from regularized import (
    ARCHIVED, FIELDS, LAMBDAS, TABLE, TABLE_SHA, WORKSPACE, continuity, errors,
    positive_status, state_at, thermodynamic_combinations,
)
from freeeos_potential import CRAD, LN10, read_table, source_constraints

OUTPUT = WORKSPACE / 'freeeos_biquintic_progress/experiments/033_first_constraint'
SOFT_OUTPUT = WORKSPACE / 'freeeos_biquintic_progress/experiments/031_regularized_c3'
FROZEN_HELPER_SHA = 'ba30d62f864e094bb52ad0b83f87554b8fff199c4146f1c1e1a6fc297b35edfd'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def summary(values):
    values = np.asarray(values)
    finite = np.isfinite(values)
    result = {'count': int(values.size), 'nonfinite': int((~finite).sum())}
    if finite.any():
        result.update(zip(('median', 'p99', 'max'), map(float, np.quantile(abs(values[finite]), [.5, .99, 1]))))
    return result


def timeout_handler(signum, frame):
    raise TimeoutError('033 per-fit 90-second limit')


def run(args):
    started = time.perf_counter()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    if (OUTPUT / 'results.json').exists() and not args.overwrite:
        raise ValueError('Products exist; use --overwrite')
    helper = BIQUINTIC / 'c3_spline_control/regularized.py'
    if digest(helper) != FROZEN_HELPER_SHA or digest(TABLE) != TABLE_SHA:
        raise ValueError('Frozen helper or pinned table mismatch')
    archived = json.loads((ARCHIVED / 'results.json').read_text())
    patch = archived['patches']['active_chiT']
    archive_path = ARCHIVED / 'active_chiT/fits_and_queries.npz'
    archive = np.load(archive_path)
    soft = json.loads((SOFT_OUTPUT / 'results.json').read_text())
    soft_arrays = np.load(SOFT_OUTPUT / 'paired_fits_and_queries.npz')
    _, logq, logt, rows = read_table(TABLE)
    low, valid, _ = source_constraints(rows, logq, logt)
    start_q, start_t, count_q, count_t = patch['selected_indices']
    patch_q, patch_t = logq[start_q:start_q + count_q], logt[start_t:start_t + count_t]
    np.testing.assert_array_equal(patch_q, archive['patch_logQ'])
    np.testing.assert_array_equal(patch_t, archive['patch_logT'])
    if not valid[start_q:start_q + count_q, start_t:start_t + count_t].all():
        raise ValueError('Active patch no longer native valid')
    grid_q, grid_t = patch_q[::2], patch_t[::2]
    coarse_low = low[:, :, start_q:start_q + count_q:2, start_t:start_t + count_t:2]
    coarse_rows = rows[start_q:start_q + count_q:2, start_t:start_t + count_t:2]
    node_q, node_t = np.meshgrid(grid_q, grid_t, indexing='ij')
    node_q, node_t = node_q.ravel(), node_t.ravel()
    source_nodes = coarse_rows.reshape(-1, 22)
    node_reference = {'Pgas': 10.0**source_nodes[:, 1], 'E': 10.0**source_nodes[:, 2], 'S': 10.0**source_nodes[:, 3],
        'Cv': source_nodes[:, 7], 'Cp': source_nodes[:, 6], 'chiRho': source_nodes[:, 4], 'chiT': source_nodes[:, 5]}
    node_reference['P'] = node_reference['Pgas'] + CRAD * (10.0**node_t)**4 / 3
    withheld_q, withheld_t = archive['withheld_logQ'], archive['withheld_logT']
    np.testing.assert_array_equal(withheld_q, soft_arrays['withheld_logQ'])
    np.testing.assert_array_equal(withheld_t, soft_arrays['withheld_logT'])
    reference = {field: soft_arrays['printed_reference_' + field] for field in FIELDS}
    metrics = archive['metrics'].tolist()[:11]
    c4 = {field: archive['thermo_0.1_withheld_values'][metrics.index(field)] for field in FIELDS}
    source_first = thermodynamic_combinations(coarse_low)
    products = {'coarse_logQ': grid_q, 'coarse_logT': grid_t, 'node_logQ': node_q, 'node_logT': node_t,
        'withheld_logQ': withheld_q, 'withheld_logT': withheld_t, 'withheld_iq': soft_arrays['withheld_iq'],
        'withheld_it': soft_arrays['withheld_it'], 'archived_C4_thermo0.1_withheld_values': archive['thermo_0.1_withheld_values']}
    for field in FIELDS:
        products['printed_reference_' + field] = reference[field]
        products['retained_reference_' + field] = node_reference[field]
    result = {'claim': 'One active geometry coarse control with exact E,TS,P/rho constraints, soft second responses and explicit roughness selector. Table-only fidelity, not fresh-source accuracy or selected runtime.',
        'table_sha256': TABLE_SHA, '015_archive_sha256': digest(archive_path), '031_results_sha256': digest(SOFT_OUTPUT / 'results.json'),
        '031_paired_array_sha256': digest(SOFT_OUTPUT / 'paired_fits_and_queries.npz'),
        'source_hashes': {str(path.relative_to(WORKSPACE)): digest(path) for path in [HERE / 'constrained.py', HERE / 'run_control.py',
            HERE / 'test_constrained.py', helper, BIQUINTIC / 'c3_spline_control/spline_control.py', BIQUINTIC / 'freeeos_potential.py']},
        'frozen031_helper_sha256': FROZEN_HELPER_SHA, 'command': sys.argv, 'timestamp': datetime.now().astimezone().isoformat(),
        'versions': {'python': sys.version, 'numpy': np.__version__, 'scipy': version('scipy')},
        'patch': {'label': 'active_chiT', 'bounds_from015': patch['actual_bounds'], 'source_coarse_shape': list(coarse_rows.shape[:2]),
                  'source_nodes': len(node_q), 'same_heldout_count': len(withheld_q), 'active_fraction_note': 'Geometry only, not measured active MESA fraction.'},
        'objective': 'Hard exact E/TS/P/rho at retained nodes; remaining weighted-relative response residual sum/(six*nodes) + lambda*031 normalized integrated normal-fourth roughness, response importance .1.',
        'rank_policy': 'Independent hard-row full rank required, reduced-selector full rank required; hard E/S/Pgas discrepancy <=1e-9 relative and reduced gradient <=1e-10. Stop all remaining fits on any hard constraint failure.',
        'bounds': {'unknowns': 1064, 'hard_equations': 702, 'lambdas': list(LAMBDAS), 'max_fits': 3, 'per_fit_seconds': 90, 'no_new_solve_after_local_hour': 20, 'threads': 16},
        'conventions': 'Same total potential/current CRAD and printed E/S gauges as015/031; no fitted constant, gauge adjustment, hardsecond-constraint claim or source accuracy.',
        'C4_archived_heldout_errors': errors(c4, reference), '031_soft_heldout_errors': {label:fit['heldout_relative_errors_NOT_source_accuracy'] for label,fit in soft['fits']['coarse'].items()},
        'fits': {}}
    for regularization in LAMBDAS:
        label = f'lambda_{regularization:g}'
        if datetime.now().hour >= 20:
            result['fits'][label] = {'status': 'deferred_deadline'}
            continue
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.setitimer(signal.ITIMER_REAL, 90)
        try:
            coefficients, solve = solve_constrained(coarse_low, grid_q * LN10, grid_t * LN10, regularization)
        except (HardConstraintFailure, ValueError, RuntimeError, TimeoutError, np.linalg.LinAlgError) as error:
            result['fits'][label] = {'status': 'failed_stop', 'reason': str(error),
                'diagnostics': error.diagnostics if isinstance(error, HardConstraintFailure) else None}
            break
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        nodal, native = state_at(coefficients, grid_q, grid_t, node_q, node_t)
        node_error = errors(nodal, node_reference)
        first_error = {field:node_error[field] for field in ['E', 'S', 'P', 'Pgas']}
        if max(check['max'] for check in first_error.values()) > 1e-9:
            result['fits'][label] = {'status': 'failed_stop', 'reason': 'Physical hard E/S/Pgas guard >1e-9', 'diagnostics': first_error}
            break
        fitted = thermodynamic_combinations(native)
        reconciliation = {field:summary((fitted[field] - target.ravel()) /
            np.maximum(abs(target.ravel()), .01 * np.median(abs(target)))) for field,target in source_first.items()}
        state, _ = state_at(coefficients, grid_q, grid_t, withheld_q, withheld_t)
        heldout_error = errors(state, reference)
        pressure = state['P']
        rho = 10.0**(withheld_q + 2 * withheld_t - 12)
        temperature = 10.0**withheld_t
        source_pressure = reference['P']
        source_pressure_rho = source_pressure / rho * reference['chiRho']
        source_pressure_t = source_pressure / temperature * reference['chiT']
        partial_errors = {'P_rho': summary((state['P_rho'] - source_pressure_rho) / abs(source_pressure_rho)),
                          'P_T': summary((state['P_T'] - source_pressure_t) / abs(source_pressure_t))}
        paired = {}
        for field in FIELDS:
            soft_state = soft_arrays[f'coarse_{label}_withheld_{field}']
            hard_error = abs((state[field] - reference[field]) / reference[field])
            soft_error = abs((soft_state - reference[field]) / reference[field])
            paired[field] = {'fraction_improved': float(np.mean(hard_error < soft_error)),
                             'fraction_worsened': float(np.mean(hard_error > soft_error))}
            products[label + '_withheld_' + field] = state[field]
            products[label + '_node_' + field] = nodal[field]
            products[label + '_soft031_withheld_' + field] = soft_state
        diagnostics = {'status': 'solved', 'solve': solve, 'physical_first_nodal_reconciliation': first_error,
            'six_combination_nodal_reconciliation': reconciliation, 'heldout_relative_errors_NOT_source_accuracy': heldout_error,
            'heldout_pressure_partial_relative_errors': partial_errors,
            'heldout_positivity': positive_status(state), 'nodal_positivity': positive_status(nodal),
            'continuity': continuity(coefficients, grid_q, grid_t), 'paired_hard_vs031soft': paired}
        result['fits'][label] = diagnostics
        products[label + '_coefficients'] = coefficients
        products[label + '_withheld_P_rho'] = state['P_rho']
        products[label + '_withheld_P_T'] = state['P_T']
        print(f'{label}: rank {solve["hard_rank"]}/{solve["hard_rows"]}, reduced {solve["reduced_rank"]}/{solve["nullspace_columns"]}, first_max={max(check["max"] for check in first_error.values()):.3e}, Cv p99={heldout_error["Cv"]["p99"]:.3e}', flush=True)
    result['elapsed_seconds'] = time.perf_counter() - started
    result['peak_process_rss_kib_linux'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if digest(helper) != FROZEN_HELPER_SHA:
        raise ValueError('031 frozen helper changed during run')
    (OUTPUT / 'results.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    np.savez_compressed(OUTPUT / 'paired_controls.npz', **products)
    lines = ['# Exact-first C3 roughness-selected control', '', result['claim'], '',
        'Same204 archived active-geometry heldouts; no new C4/031 fits and no gauge changes.', '',
        '| Candidate | Cv p99 / max | Pgas p99 / max | E p99 / max | S p99 / max |',
        '| --- | --- | --- | --- | --- |']
    entries = [('archived_C4_thermo0.1', result['C4_archived_heldout_errors'])]
    for label, fit in result['fits'].items():
        if fit['status']=='solved':
            entries.extend([(label + '_031soft', result['031_soft_heldout_errors'][label]),
                            (label + '_033hard', fit['heldout_relative_errors_NOT_source_accuracy'])])
    for label, metrics in entries:
        lines.append('| ' + label + ' | ' + ' | '.join(f"{metrics[field]['p99']:.3e} / {metrics[field]['max']:.3e}" for field in ['Cv','Pgas','E','S']) + ' |')
    lines.extend(['', 'Hard retained-node E/S/Pgas checks, full hard/reduced ranks, exact-constraint roundoff, pressure-partial holdouts and paired improvement fractions are in JSON.',
                  'No universal winner, runtime choice, offnode source truth, fresh-source accuracy or global positivity bound.'])
    (OUTPUT / 'README.md').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
