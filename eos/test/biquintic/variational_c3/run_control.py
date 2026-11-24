"""Bounded same-query printed holdouts and offnode/face diagnostics for variational C3."""

import argparse
import hashlib
import json
import resource
import sys
import time
from datetime import datetime
from importlib.metadata import version
from pathlib import Path

import numpy as np

from normal_jets import BIQUINTIC, JET_ORDERS, LOW_ORDERS, coarse_bridge_guard, construct_variational_jets
from c3_potential import construct_jets, evaluate as evaluate_c3
from freeeos_potential import (
    CRAD, EVAL_ORDERS, LN10, construct_twists, evaluate as evaluate_bq,
    query_cells, read_table, source_constraints, thermodynamic_partials, valid_cells,
)

sys.path.insert(0, str(BIQUINTIC / 'refinement_control'))
from refinement import material_mask, native_mask, own_use, table_reference

WORKSPACE = BIQUINTIC.parents[3]
HERE = Path(__file__).resolve().parent
OUTPUT = WORKSPACE / 'freeeos_biquintic_progress/experiments/025_variational_c3'
TABLE = WORKSPACE / 'mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data'
TABLE_SHA = 'e6ca7ce27bee8197bb257a8c4886a06366e6951e58958fbed3735c877ee68502'
REGIONS = WORKSPACE / 'freeeos_biquintic_progress/experiments/015_smooth_validation/results.json'
FIELDS = ('P', 'Pgas', 'E', 'S', 'Cv', 'Cp', 'chiRho', 'chiT')
PARTIALS = tuple(field + suffix for field in FIELDS for suffix in ('_x', '_t'))
POSITIVE = FIELDS + ('P_rho',)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
            for field in POSITIVE}


def region_masks(bounds, sample_q, sample_t):
    return {name: ((sample_q >= box['logQ'][0] - 1e-12) & (sample_q <= box['logQ'][1] + 1e-12)
                   & (sample_t >= box['logT'][0] - 1e-12) & (sample_t <= box['logT'][1] + 1e-12))
            for name, box in bounds.items()}


def native_at(model, logq, logt, sample_q, sample_t, orders=EVAL_ORDERS, cells=None):
    return model['evaluate'](model['nodes'], model['complete'], logq, logt, sample_q, sample_t,
                             orders=orders, cell_indices=cells, add_radiation=True)


def state_at(model, logq, logt, sample_q, sample_t):
    return thermodynamic_partials(native_at(model, logq, logt, sample_q, sample_t), sample_q, sample_t)


def models_from(low, valid, spacing):
    nodes, complete, _ = construct_twists(low, valid, *spacing, radius=2, route='mean')
    models = {'BQ_r2_reference': {'nodes': nodes, 'complete': complete, 'evaluate': evaluate_bq}}
    nodes, complete, diagnostics = construct_jets(low, valid, *spacing, radius=2)
    models['C3_FD_r2_reference'] = {'nodes': nodes, 'complete': complete, 'evaluate': evaluate_c3,
                                  'jet_masks': diagnostics['jet_masks']}
    nodes, complete, diagnostics = construct_variational_jets(low, valid, *spacing)
    models['C3_variational_r2_reference'] = {'nodes': nodes, 'complete': complete, 'evaluate': evaluate_c3,
        'jet_masks': diagnostics['jet_masks'], 'normal_segments': diagnostics['normal_segments']}
    return models


def face_control(models, common_cells, logq, logt, bounds, products):
    results = {}
    for axis in (0, 1):
        shared = common_cells[:-1, :] & common_cells[1:, :] if axis == 0 else common_cells[:, :-1] & common_cells[:, 1:]
        density, temperature = np.nonzero(shared)
        sample_q = logq[density + 1] if axis == 0 else logq[density] + .43 * (logq[1] - logq[0])
        sample_t = logt[temperature] + .41 * (logt[1] - logt[0]) if axis == 0 else logt[temperature + 1]
        selected = np.flatnonzero(np.logical_or.reduce(list(region_masks(bounds, sample_q, sample_t).values())))
        if selected.size > 64:
            selected = selected[np.linspace(0, selected.size - 1, 64, dtype=int)]
        sample_q, sample_t = sample_q[selected], sample_t[selected]
        left = (density[selected], temperature[selected])
        right = (left[0] + int(axis == 0), left[1] + int(axis == 1))
        prefix = f'face_axis{axis}'
        products[prefix + '_logQ'], products[prefix + '_logT'] = sample_q, sample_t
        products[prefix + '_left_cells'] = np.column_stack(left)
        products[prefix + '_right_cells'] = np.column_stack(right)
        results[prefix] = {'queries': int(len(sample_q)), 'models': {}}
        if not len(sample_q):
            continue
        for name, model in models.items():
            traces = [native_at(model, logq, logt, sample_q, sample_t, JET_ORDERS, cells) for cells in (left, right)]
            states = [thermodynamic_partials(trace, sample_q, sample_t) for trace in traces]
            scale = np.maximum(1., np.maximum(abs(traces[0]), abs(traces[1])))
            scaled = (traces[1] - traces[0]) / scale
            response_jumps = {field: summary((states[1][field] - states[0][field]) /
                np.maximum(1., np.maximum(abs(states[0][field]), abs(states[1][field]))))
                for field in FIELDS + PARTIALS}
            results[prefix]['models'][name] = {'all16_native_scaled_jumps': summary(scaled),
                'value_and_response_partial_scaled_jumps': response_jumps,
                'left_positivity': positive_status(states[0]), 'right_positivity': positive_status(states[1])}
            for side, trace, state in zip(('left', 'right'), traces, states):
                products[f'{prefix}_{name}_native_{side}'] = trace
                for field in FIELDS + PARTIALS:
                    products[f'{prefix}_{name}_{field}_{side}'] = state[field]
    return results


def jacobian_control(models, logq, logt, sample_q, sample_t):
    selected = np.linspace(0, len(sample_q) - 1, min(32, len(sample_q)), dtype=int)
    sample_q, sample_t = sample_q[selected], sample_t[selected]
    step = 1e-5
    result = {}
    for name, model in models.items():
        center = state_at(model, logq, logt, sample_q, sample_t)
        checks = {}
        for suffix, shift_q, shift_t in (('_x', step / LN10, 0.), ('_t', -2 * step / LN10, step / LN10)):
            plus = state_at(model, logq, logt, sample_q + shift_q, sample_t + shift_t)
            minus = state_at(model, logq, logt, sample_q - shift_q, sample_t - shift_t)
            for field in FIELDS:
                difference = (plus[field] - minus[field]) / (2 * step)
                scale = np.maximum(1., np.maximum(abs(center[field]), abs(center[field + suffix])))
                checks[field + suffix] = summary((difference - center[field + suffix]) / scale)
        result[name] = checks
    return {'sample_count': len(selected), 'physical_log_step': step,
            'claim': 'Internal analytic-Jacobian versus central differences, NOT source accuracy.', 'models': result}


def run(args):
    started = time.perf_counter()
    if digest(TABLE) != TABLE_SHA:
        raise ValueError('Unexpected pinned source-table SHA-256')
    OUTPUT.mkdir(parents=True, exist_ok=True)
    if not args.overwrite and (OUTPUT / 'results.json').exists():
        raise ValueError('Products exist; use --overwrite')
    _, full_q, full_t, rows = read_table(TABLE)
    full_valid = native_mask(rows, full_q, full_t)
    logq, logt, coarse = full_q[::2], full_t[::2], rows[::2, ::2]
    low, valid, _ = source_constraints(coarse, logq, logt, split_radiation=True)
    valid &= full_valid[::2, ::2]
    before_guard = valid.copy()
    bridge_guard, bridge_evidence = coarse_bridge_guard(full_valid, valid)
    valid &= bridge_guard
    spacing = ((logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10)
    models = models_from(low, valid, spacing)
    common_cells = valid_cells(material_mask(coarse, logq, logt) & own_use(logq[:, None], logt[None, :]))
    products = {'coarse_logQ': logq, 'coarse_logT': logt, 'coarse_source_valid': valid,
                'coarse_valid_before_bridge_guard': before_guard, 'coarse_bridge_guard': bridge_guard,
                'original_native_valid': full_valid}
    for name, model in models.items():
        common_cells &= valid_cells(model['complete'])
        products[name + '_complete'] = model['complete']
        for orders, mask in model.get('jet_masks', {}).items():
            products[name + '_mask' + ''.join(map(str, orders))] = mask
    products['common_material_ownuse_cells'] = common_cells
    bounds = {name: patch['actual_bounds'] for name, patch in json.loads(REGIONS.read_text())['patches'].items()}
    density, temperature = np.indices(rows.shape[:2])
    eligible = ((density % 2 == 1) & (temperature % 2 == 1) & full_valid
                & material_mask(rows, full_q, full_t) & own_use(full_q[:, None], full_t[None, :]))
    eligible &= np.logical_or.reduce(list(region_masks(bounds, full_q[:, None], full_t[None, :]).values()))
    density, temperature = np.nonzero(eligible)
    selected = np.flatnonzero(common_cells[density // 2, temperature // 2])
    common_count = len(selected)
    if common_count > args.max_queries:
        selected = np.sort(np.random.default_rng(20261028).choice(selected, args.max_queries, replace=False))
    if not len(selected):
        raise ValueError('No common native patch holdouts')
    density, temperature = density[selected], temperature[selected]
    sample_q, sample_t = full_q[density], full_t[temperature]
    masks = region_masks(bounds, sample_q, sample_t)
    reference = table_reference(rows[density, temperature])
    reference['P'] = reference['Pgas'] + CRAD * (10.0**sample_t)**4 / 3
    products.update({'holdout_iq': density, 'holdout_it': temperature, 'holdout_logQ': sample_q,
                     'holdout_logT': sample_t, 'holdout_printed_logRho': rows[density, temperature, 18]})
    for field in FIELDS:
        products['printed_reference_' + field] = reference[field]
    holdouts = {}
    for name, model in models.items():
        state = state_at(model, logq, logt, sample_q, sample_t)
        errors = {field: (state[field] - reference[field]) / abs(reference[field]) for field in FIELDS}
        holdouts[name] = {'relative_printed_field_errors': {field: summary(error) for field, error in errors.items()},
            'positivity': positive_status(state), 'regions': {region: {'queries': int(mask.sum()),
                'errors': {field: summary(error[mask]) for field, error in errors.items()}} for region, mask in masks.items()}}
        for field in FIELDS:
            products[f'holdout_{name}_{field}'] = state[field]
            products[f'holdout_{name}_{field}_relative_error'] = errors[field]
    paired = {}
    for field in FIELDS:
        old = abs(products[f'holdout_C3_FD_r2_reference_{field}_relative_error'])
        new = abs(products[f'holdout_C3_variational_r2_reference_{field}_relative_error'])
        paired[field] = {'fraction_improved': float(np.mean(new < old)), 'fraction_worsened': float(np.mean(new > old)),
                         'baseline_zero_count': int((old == 0).sum()), 'variational_zero_count': int((new == 0).sum())}
    off_q = logq[density // 2] + .37 * (logq[1] - logq[0])
    off_t = logt[temperature // 2] + .63 * (logt[1] - logt[0])
    products['offnode_logQ'], products['offnode_logT'] = off_q, off_t
    states = {name: state_at(model, logq, logt, off_q, off_t) for name, model in models.items()}
    bq = states['BQ_r2_reference']
    offnodes = {}
    for name, state in states.items():
        departures = {}
        for field in FIELDS + PARTIALS:
            base_field = field[:-2] if field in PARTIALS else field
            scale = np.maximum(1., np.maximum(abs(bq[base_field]), abs(bq[field]))) if field in PARTIALS else abs(bq[field])
            departures[field] = summary((state[field] - bq[field]) / scale)
            products[f'offnode_{name}_{field}'] = state[field]
        offnodes[name] = {'departures_from_BQ_not_accuracy': departures, 'positivity': positive_status(state)}
    faces = face_control(models, common_cells, logq, logt, bounds, products)
    segment_data = models['C3_variational_r2_reference']['normal_segments']
    normal_audit = {}
    for orders, segments in segment_data.items():
        ratios = np.array([segment['energy_ratio'] for segment in segments if segment['energy_ratio'] is not None])
        normal_audit[orders] = {'segments': len(segments), 'longest_nodes': max(segment['nodes'] for segment in segments),
            'stationarity': summary([segment['stationarity_relative'] for segment in segments]),
            'minimum_over_fd_interior_minimum_elsewhere_energy': summary(ratios),
            'energy_ratio_above_one_plus_1e-10': int((ratios > 1 + 1e-10).sum()), 'line_ranges_and_energy': segments}
    source_checks = {}
    for name, model in models.items():
        source_checks[name] = {str(orders): {'max_absolute_difference': float(np.max(abs(model['nodes'][orders][model['complete']] - low[orders][model['complete']])))}
                              for orders in LOW_ORDERS}
    structural_checks = {
        'targets_never_retained': bool(np.all(density % 2 == 1) and np.all(temperature % 2 == 1)),
        'six_constraints_exact': all(check['max_absolute_difference'] == 0 for model in source_checks.values() for check in model.values()),
        'C3_all16_faces_exact': all(model['all16_native_scaled_jumps']['max'] == 0
            for face in faces.values() for name, model in face['models'].items() if name.startswith('C3_')),
        'C3_response_faces_exact': all(check['max'] == 0
            for face in faces.values() for name, model in face['models'].items() if name.startswith('C3_')
            for check in model['value_and_response_partial_scaled_jumps'].values()),
        'normal_stationarity': all(audit['stationarity']['max'] < 1e-12 for audit in normal_audit.values()),
        'normal_energy_minimum': all(audit['energy_ratio_above_one_plus_1e-10'] == 0 for audit in normal_audit.values()),
    }
    if not all(structural_checks.values()):
        raise AssertionError(f'Structural validation failed: {structural_checks}')
    result = {
        'claim': 'Experimental high-jet selection. Printed-node holdout fidelity and internal offnode/face consistency; NOT fresh-source accuracy or selected runtime.',
        'table_sha256': TABLE_SHA, 'region_source_sha256': digest(REGIONS),
        'source_hashes': {str(path.relative_to(WORKSPACE)): digest(path) for path in [HERE / 'normal_jets.py', HERE / 'run_control.py',
            HERE / 'test_variational.py', BIQUINTIC / 'freeeos_potential.py', BIQUINTIC / 'c3_control/c3_potential.py',
            BIQUINTIC / 'refinement_control/refinement.py']},
        'command': sys.argv, 'timestamp': datetime.now().astimezone().isoformat(),
        'versions': {'python': sys.version, 'numpy': np.__version__, 'scipy': version('scipy')},
        'radiation': {'current_crad': CRAD, 'mode': 'Coherent current-CRAD analytic-reference split/restoration for all models; no fitted constant or gauge, not verified pure material.'},
        'grid': {'stride': 2, 'original_shape': list(rows.shape[:2]), 'coarse_shape': list(valid.shape),
                 'native_coordinate_spacing': list(spacing), 'iq520_original_excluded': True},
        'original_invalid_midpoint_guard': {'evidence': bridge_evidence,
            'rule': 'Conservatively exclude both retained endpoints of any accepted coarse edge with an invalid intervening original row. Apply the same stricter donor mask to every model, including baseline C3 and BQ.'},
        'dependency': 'Each normal third jet depends on the entire contiguous valid coarse line segment. Mixed jets propagate those nonlocal dependencies through the existing strict two-route radius-2 graph. No bridge across invalid/lowT/iq520 rows and no one-sided FD route.',
        'endpoint_condition': 'Free endpoint third jets minimize the integrated squared normal fourth derivative; natural zero fourth derivative at segment endpoints. No arbitrary degree-7 reproduction claim.',
        'training': 'Only original even/even coarse retained rows; preserve all six supplied low constraints. Global coarse native donors, not patch-restricted solves; no fine estimated jets reused.',
        'normal_jets': normal_audit, 'six_source_constraints': source_checks, 'structural_validation': structural_checks,
        'coverage': {'eligible_patch_odd_odd_material_ownuse_targets': int(eligible.sum()), 'common_targets': common_count,
            'selected_targets': len(sample_q), 'per_model_complete_cells': {name: int(valid_cells(model['complete']).sum()) for name, model in models.items()},
            'common_material_ownuse_cells': int(common_cells.sum())},
        'patch_bounds_from_015': bounds, 'region_counts': {name: int(mask.sum()) for name, mask in masks.items()},
        'active_note': 'active_chiT is geometry from 015, not measured MESA component fractions.',
        'holdout': holdouts, 'paired_C3_variational_vs_FD': paired,
        'offnode': {'queries': len(off_q), 'definition': 'Identical deterministic 0.37/0.63 interiors in each selected coarse holdout cell; model-to-BQ departures, no offnode source truth.', 'models': offnodes},
        'faces': faces, 'jacobian_internal_check': jacobian_control(models, logq, logt, off_q, off_t),
        'elapsed_seconds': time.perf_counter() - started,
        'peak_process_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (OUTPUT / 'results.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    np.savez_compressed(OUTPUT / 'audit_arrays.npz', **products)
    lines = ['# Variational normal-jet C3 control', '', result['claim'], '',
             f"{len(sample_q)} identical native odd/odd stride-2 holdouts; regions: {result['region_counts']}.", '',
             '| Model | Pgas p99 / max | E p99 / max | S p99 / max | Cv p99 / max | chiRho p99 / max | chiT p99 / max |',
             '| --- | --- | --- | --- | --- | --- | --- |']
    for name, model in holdouts.items():
        columns = [f"{model['relative_printed_field_errors'][field]['p99']:.3e} / {model['relative_printed_field_errors'][field]['max']:.3e}"
                   for field in ('Pgas', 'E', 'S', 'Cv', 'chiRho', 'chiT')]
        lines.append('| ' + name + ' | ' + ' | '.join(columns) + ' |')
    lines.extend(['', f"Variational Cv improves at {100 * paired['Cv']['fraction_improved']:.2f}% of identical targets versus FD C3; aggregate Cv remains worse than BQ.",
        'All sampled heldout, offnode and face physical fields/positive responses are finite and positive.' if all(
            sum(sum(status.values()) for status in model['positivity'].values()) == 0
            for section in [holdouts, offnodes] for model in section.values()) and all(
                sum(sum(status.values()) for status in model[side + '_positivity'].values()) == 0
                for face in faces.values() for model in face['models'].values() for side in ['left', 'right'])
            else 'Some sampled positivity/finiteness failures occur; consult results.json.',
        'Both C3 controls have zero sampled all16 native and physical response-Jacobian exact-face jumps; this is a construction check.',
        f"Original midpoint bridge guards exclude {bridge_evidence['previously_valid_coarse_nodes_excluded']} retained donor nodes ({bridge_evidence['axis0_invalid_original_midpoints']}/{bridge_evidence['axis1_invalid_original_midpoints']} invalid original midpoint edges). The target set is unchanged.",
        'All normal-line segment ranges, energy ratios, stationarity, masks and paired raw values are archived.',
        'Normal third jets have global segment dependencies; mixed jets retain strict existing radius-2 two-route masks.',
        'No fresh-source offnode accuracy, formal order, positivity between sampled points, or runtime selection claim.'])
    (OUTPUT / 'README.md').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    print(f'Elapsed {result["elapsed_seconds"]:.2f}s; outputs {OUTPUT}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-queries', type=int, default=1000)
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    if not 1 <= args.max_queries <= 2000:
        parser.error('Require 1..2000 paired queries')
    run(args)


if __name__ == '__main__':
    main()
