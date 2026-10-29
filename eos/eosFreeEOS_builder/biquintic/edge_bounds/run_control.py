"""New native shifted-edge cells: bounded split-potential whole-cell sign audit."""

import argparse
from datetime import datetime
import json
import os
from pathlib import Path
import resource
import sys
import time
from zoneinfo import ZoneInfo

sys.dont_write_bytecode = True
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import numpy as np

HERE = Path(__file__).resolve().parent
BIQUINTIC = HERE.parent
WORKSPACE = BIQUINTIC.parents[3]
sys.path.insert(0, str(BIQUINTIC / 'positivity_bounds'))
sys.path.insert(0, str(BIQUINTIC / 'edge_control'))
sys.path.insert(0, str(BIQUINTIC))

from bounds import (DTYPE, FIELDS, STATUS, bernstein_bounds, cell_power, domain_masks,
                    normalize_closed_cell, polynomial_value, sign_powers)
from run_bounds import PHASES, TABLE, TABLE_SHA, collect_jets, location
from bound_split import adaptive_split, radiation_cv_numerator, radiation_guard, split_witness
from clean_floor import build_control
from freeeos_potential import CRAD, LN10, checksum, evaluate, read_table, thermodynamics, valid_cells

OUTPUT = WORKSPACE / 'freeeos_biquintic_progress/experiments/035_edge_bounds'


def finish_fence(now, requested=None):
    local = now.astimezone(ZoneInfo('US/Pacific'))
    if requested is None:
        fence = local.replace(hour=20, minute=0, second=0, microsecond=0)
    else:
        fence = datetime.fromisoformat(requested)
        if fence.tzinfo is None:
            raise ValueError('--finish-by requires a timezone-aware ISO timestamp')
        fence = fence.astimezone(ZoneInfo('US/Pacific'))
    if (local.year, local.month, local.day) == (2026, 10, 28):
        fence = min(fence, local.replace(hour=20, minute=0, second=0, microsecond=0))
    return fence


def check_limits(max_cells, max_depth, max_nodes, seconds):
    if not 1 <= max_cells <= 20000 or not 0 <= max_depth <= 3 or not 1 <= max_nodes <= 85 or not 1 <= seconds <= 180:
        raise ValueError('Limits: cells 1..20000, depth 0..3, nodes 1..85, seconds 1..180')


def new_cell_mask(complete, source_valid, baseline_complete, max_cells):
    if complete.shape != source_valid.shape or complete.shape != baseline_complete.shape:
        raise ValueError('Candidate/source/baseline node-mask shape mismatch')
    candidate = valid_cells(complete)
    if np.any(candidate & ~valid_cells(source_valid)):
        raise ValueError('Complete candidate cell lacks all four valid native source corners')
    selected = candidate & ~valid_cells(baseline_complete)
    if int(selected.sum()) > max_cells:
        raise ValueError(f'{int(selected.sum())} new cells exceed {max_cells}; no silent truncation')
    return selected


def support_guard(support, logt, floor):
    if floor not in (3.1, 3.2) or not np.isfinite(floor):
        raise ValueError('Only declared clean source floors 3.1/3.2 allowed')
    if support.size == 0:
        return None
    if np.any(support < 0) or np.any(support[..., 2] > support[..., 3]) or np.any(support[..., 3] >= len(logt)):
        raise ValueError('Invalid donor-support bounds')
    minimum = float(logt[support[..., 2]].min())
    if minimum < floor - 1e-12 or minimum < 3.1 - 1e-12:
        raise ValueError('Recursive support crosses declared floor or known below-3.1 exclusion')
    return minimum


def guarded_powers(jets, du, dtau, density_stretch, temperature_stretch):
    if not np.isfinite(du) or not np.isfinite(dtau) or min(du, dtau) <= 0 or not np.isfinite(jets).all():
        raise ValueError('Nonfinite jets or unsafe native spacing; no certificate')
    power, magnitude = cell_power(jets, du, dtau)
    powers, scales = sign_powers(power, magnitude, du, dtau)
    for stretch in (density_stretch, temperature_stretch):
        if not np.isfinite(stretch).all() or np.any(stretch <= 0):
            raise ValueError('Invalid closed-cell stretch')
    powers = {field: normalize_closed_cell(values, density_stretch, temperature_stretch) for field, values in powers.items()}
    scales = {field: normalize_closed_cell(values, density_stretch, temperature_stretch) for field, values in scales.items()}
    if any(not np.isfinite(values).all() for values in list(powers.values()) + list(scales.values())):
        raise ValueError('Nonfinite assembled polynomial or absolute magnitude guard')
    return powers, scales


def field_summary(record, mask):
    status = record['status'][mask]
    return {'cells': int(mask.sum()), 'positive': int((status == 1).sum()),
            'negative_witness': int((status == 2).sum()), 'unresolved': int((status == 0).sum()),
            'visited_nodes': int(record['visited_nodes'][mask].sum()),
            'unresolved_leaves': int(record['unresolved_leaves'][mask].sum()),
            'maximum_depth': int(record['maximum_depth'][mask].max()) if mask.any() else 0,
            'root_representation_positive_but_restored_unresolved': int(((status == 0) & record['root_representation_positive'][mask]).sum())}


def audit_candidate(control, baseline_complete, logq, logt, rows, floor, args, deadline):
    selected = new_cell_mask(control['complete'], control['valid'], baseline_complete, args.max_cells)
    density_indices, temperature_indices = np.nonzero(selected)
    count = len(density_indices)
    ratio = CRAD * (10.0**logt[None, :])**4 / (3 * 10.0**rows[:, :, 1])
    own_grid, material_grid, branches, ratio_max = domain_masks(logq, logt, ratio)
    own = own_grid[selected]
    support = np.stack([control['metadata']['support'][density_indices + density_end, temperature_indices + temperature_end]
                        for density_end in range(2) for temperature_end in range(2)], axis=1)
    minimum_donor = support_guard(support, logt, floor)
    if count and float(logt[temperature_indices].min()) < floor - 1e-12:
        raise ValueError('Selected cell crosses source floor')
    du, dtau = (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10
    archive = {'density_cell_index': density_indices, 'temperature_cell_index': temperature_indices,
               'source_floor': np.array(floor), 'baseline_floor': np.array(3.2),
               'full_one_branch_own_use': own, 'material_corner_screen': material_grid[selected],
               'source_corner_valid': np.stack([control['valid'][density_indices + density_end, temperature_indices + temperature_end]
                    for density_end in range(2) for temperature_end in range(2)], axis=1),
               'donor_corner_support': support, 'logQ_axis': logq, 'logT_axis': logt,
               'selected_new_cell_mask': selected, 'candidate_complete_node_mask': control['complete'],
               'candidate_source_valid_node_mask': control['valid'], 'baseline_centered_complete_node_mask': baseline_complete,
               'PHASES': np.array(PHASES), 'sample_logQ': np.full((count, 4), np.nan), 'sample_logT': np.full((count, 4), np.nan)}
    records = {}
    witnesses = {field: {} for field in FIELDS}
    agreement = {field: {'queries': 0, 'outside_effective_margin': 0, 'nonfinite_core_queries': 0,
                         'negative_core_queries': 0, 'maximum_difference_over_effective_margin': 0., 'worst': None} for field in FIELDS}
    for field in FIELDS:
        records[field] = {name: np.zeros(count, dtype=dtype) for name, dtype in (
            ('status', np.int8), ('visited_nodes', np.int16), ('maximum_depth', np.int8),
            ('unresolved_leaves', np.int16), ('positive_leaves', np.int16), ('witness_core_queries', np.int16),
            ('root_representation_positive', bool), ('processed', bool), ('root_lower_bound', DTYPE),
            ('final_lower_bound', DTYPE), ('polynomial_margin', DTYPE), ('analytic_margin', DTYPE),
            ('restoration_margin', DTYPE), ('effective_margin', DTYPE), ('root_upper_bound', DTYPE))}
        for name in ('core_value', 'polynomial_plus_radiation_value', 'difference_over_effective_margin'):
            archive[field + '_' + name] = np.full((count, 4), np.nan, dtype=DTYPE)
    started = time.monotonic()
    completed = 0
    for start in range(0, count, args.batch_size):
        if time.monotonic() >= deadline:
            break
        stop = min(count, start + args.batch_size)
        density_batch, temperature_batch = density_indices[start:stop], temperature_indices[start:stop]
        density_width = logq[density_batch + 1].astype(DTYPE) - logq[density_batch].astype(DTYPE)
        temperature_width = logt[temperature_batch + 1].astype(DTYPE) - logt[temperature_batch].astype(DTYPE)
        powers, scales = guarded_powers(collect_jets(control['nodes'], density_batch, temperature_batch), du, dtau,
                                       density_width / DTYPE(logq[1] - logq[0]), temperature_width / DTYPE(logt[1] - logt[0]))
        radiation_min = radiation_cv_numerator(logq[density_batch + 1], logt[temperature_batch])
        radiation_max = radiation_cv_numerator(logq[density_batch], logt[temperature_batch + 1])
        effective_margins = {}
        for field in FIELDS:
            coefficients, margin, assembly, core_margin, absolute_scale = bernstein_bounds(powers[field], scales[field])
            analytic = radiation_guard(field, radiation_max, False)
            restoration = radiation_guard(field, radiation_max, True) - analytic
            effective = margin + analytic + restoration
            if not np.isfinite(coefficients).all() or not np.isfinite(effective).all() or np.any(effective <= 0):
                raise ValueError('Nonfinite/zero bounds margin; refuse numerical certificate')
            addition_min = radiation_min if field == 'Cv_numerator' else 0
            addition_max = radiation_max if field == 'Cv_numerator' else 0
            representation_lower = coefficients.min(axis=(1, 2)) + addition_min - margin - analytic
            lower = representation_lower - restoration
            record = records[field]
            for name, values in [('polynomial_margin', margin), ('analytic_margin', analytic), ('restoration_margin', restoration),
                                 ('effective_margin', effective), ('root_lower_bound', lower), ('final_lower_bound', lower),
                                 ('root_upper_bound', coefficients.max(axis=(1, 2)) + addition_max + effective)]:
                record[name][start:stop] = values
            effective_margins[field] = effective
            record['status'][start:stop] = (lower > 0).astype(np.int8)
            record['root_representation_positive'][start:stop] = representation_lower > 0
            record['visited_nodes'][start:stop] = 1
            record['positive_leaves'][start:stop] = (lower > 0).astype(np.int16)
            record['processed'][start:stop] = True
            for local_index in np.flatnonzero(lower <= 0):
                index = start + local_index
                if time.monotonic() >= deadline:
                    record['unresolved_leaves'][index] = 1
                    continue
                density_index, temperature_index = int(density_batch[local_index]), int(temperature_batch[local_index])
                box = tuple(map(DTYPE, (logq[density_index], logq[density_index + 1], logt[temperature_index], logt[temperature_index + 1])))
                callback = split_witness(field, density_index, temperature_index, control['nodes'], control['complete'], logq, logt)
                bounded = adaptive_split(coefficients[local_index], margin[local_index], box, field, callback,
                                         restored_core=True, max_depth=args.max_depth, max_nodes=args.max_nodes)
                record['status'][index] = STATUS[bounded['status']]
                record['final_lower_bound'][index] = bounded['lower_bound']
                for name in ('visited_nodes', 'maximum_depth', 'positive_leaves', 'unresolved_leaves', 'witness_core_queries'):
                    record[name][index] = bounded[name]
                if bounded['status'] == 'negative_witness':
                    witnesses[field][index] = {**location(index, density_indices, temperature_indices, logq, logt),
                                              'source_floor': floor, 'full_one_branch_own_use': bool(own[index]), **bounded}
        for phase_index, (density_phase, temperature_phase) in enumerate(PHASES):
            sample_q = np.asarray(logq[density_batch] + density_width * density_phase, dtype=float)
            sample_t = np.asarray(logt[temperature_batch] + temperature_width * temperature_phase, dtype=float)
            archive['sample_logQ'][start:stop, phase_index] = sample_q
            archive['sample_logT'][start:stop, phase_index] = sample_t
            actual_density = (sample_q.astype(DTYPE) - logq[density_batch].astype(DTYPE)) / density_width
            actual_temperature = (sample_t.astype(DTYPE) - logt[temperature_batch].astype(DTYPE)) / temperature_width
            native = evaluate(control['nodes'], control['complete'], logq, logt, sample_q, sample_t, add_radiation=True,
                              cell_indices=(density_batch, temperature_batch))
            state = thermodynamics(native, sample_q, sample_t)
            for field in FIELDS:
                prediction = polynomial_value(powers[field], actual_density, actual_temperature)
                if field == 'Cv_numerator':
                    prediction += radiation_cv_numerator(sample_q, sample_t)
                core_value = state['Cv'] * state['T'] if field == 'Cv_numerator' else state['P_rho']
                difference = abs(prediction - core_value) / effective_margins[field]
                archive[field + '_core_value'][start:stop, phase_index] = core_value
                archive[field + '_polynomial_plus_radiation_value'][start:stop, phase_index] = prediction
                archive[field + '_difference_over_effective_margin'][start:stop, phase_index] = difference
                stats = agreement[field]
                stats['queries'] += len(sample_q)
                stats['outside_effective_margin'] += int((difference > 1).sum())
                stats['nonfinite_core_queries'] += int((~np.isfinite(core_value)).sum())
                stats['negative_core_queries'] += int((core_value < 0).sum())
                if np.any((records[field]['status'][start:stop] == 1) & (core_value <= 0)):
                    raise AssertionError('Restored guarded certificate contradicts a core sign query')
                worst = int(np.argmax(difference))
                if difference[worst] > stats['maximum_difference_over_effective_margin']:
                    stats['maximum_difference_over_effective_margin'] = float(difference[worst])
                    stats['worst'] = {**location(start + worst, density_indices, temperature_indices, logq, logt),
                                      'logQ': float(sample_q[worst]), 'logT': float(sample_t[worst]),
                                      'polynomial_plus_radiation': float(prediction[worst]), 'restored_core': float(core_value[worst]),
                                      'effective_margin': float(effective_margins[field][worst])}
                for local_index in np.flatnonzero(np.isfinite(core_value) & (core_value < 0)):
                    index = start + int(local_index)
                    records[field]['status'][index] = 2
                    witnesses[field].setdefault(index, {**location(index, density_indices, temperature_indices, logq, logt),
                        'source_floor': floor, 'full_one_branch_own_use': bool(own[index]), 'status': 'negative_witness',
                        'witness': {'logQ': float(sample_q[local_index]), 'logT': float(sample_t[local_index]),
                                    'core_field': 'Cv' if field == 'Cv_numerator' else 'P_rho',
                                    'core_value': float(state['Cv'][local_index] if field == 'Cv_numerator' else core_value[local_index]),
                                    'origin': 'one of four restored-core calibration queries'}})
        completed = stop
        print(f'floor {floor:g}: {completed}/{count} new cells', flush=True)
    calibrated = completed == count and all(stats['outside_effective_margin'] == 0 and stats['nonfinite_core_queries'] == 0 for stats in agreement.values())
    if not calibrated:
        for record in records.values():
            record['status'][record['status'] == 1] = 0
    subsets = {'all_new_cells': np.ones(count, dtype=bool), 'full_one_branch_own_use': own,
               'not_full_one_branch_own_use': ~own}
    summaries = {name: {field: field_summary(records[field], mask) for field in FIELDS} for name, mask in subsets.items()}
    for name, mask in subsets.items():
        summaries[name]['both_positive'] = int((mask & (records['Cv_numerator']['status'] == 1) & (records['P_rho']['status'] == 1)).sum())
    for field, record in records.items():
        for name, values in record.items():
            archive[field + '_' + name] = values
    unresolved = {field: [{**location(int(index), density_indices, temperature_indices, logq, logt),
                           'full_one_branch_own_use': bool(own[index]), 'source_floor': floor,
                           'visited_nodes': int(records[field]['visited_nodes'][index]),
                           'unresolved_leaves': int(records[field]['unresolved_leaves'][index]),
                           'root_representation_positive': bool(records[field]['root_representation_positive'][index]),
                           'root_lower_bound': float(records[field]['root_lower_bound'][index]),
                           'final_lower_bound': float(records[field]['final_lower_bound'][index])}
                          for index in np.flatnonzero(records[field]['status'] == 0)] for field in FIELDS}
    output = OUTPUT / f'floor_{floor:g}'
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output / 'cell_bounds_and_core_queries.npz', **archive)
    result = {'source_floor': floor, 'radiation_split': True, 'new_cells': count, 'fully_sampled_cells': completed,
              'status': 'complete' if completed == count else 'incomplete_runtime_cap', 'calibration_passed': calibrated,
              'support': control['support'], 'minimum_selected_donor_logT': minimum_donor,
              'baseline': 'fixed centered radius2 mean source floor3.2, all four complete/native source corners',
              'classification': 'exact 024 domain_masks full_one_branch_own_use; geometry only, no actual MESA fraction inference',
              'subsets': summaries, 'restored_core_agreement': agreement,
              'negative_witnesses': {field: list(values.values()) for field, values in witnesses.items()},
              'unresolved_cells': unresolved, 'elapsed_bound_seconds': time.monotonic() - started,
              'archive_sha256': checksum(output / 'cell_bounds_and_core_queries.npz')}
    (output / 'results.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    return result


def write_summary(result):
    lines = ['# Experiment 035: new shifted native edge-cell sign bounds', '',
             'Whole-cell numerical sign bounds on fixed float64 jets, not source accuracy, actual EOS activity, or runtime authorization.', '',
             '| Source floor | New-cell subset | Cells | Cv: positive / witness / unresolved | Prho: positive / witness / unresolved | Both positive |',
             '|---|---|---:|---:|---:|---:|']
    for floor, candidate in result['candidates'].items():
        for subset in ('all_new_cells', 'full_one_branch_own_use'):
            row = candidate['subsets'][subset]
            values = [' / '.join(str(row[field][name]) for name in ('positive', 'negative_witness', 'unresolved')) for field in FIELDS]
            lines.append(f'| {floor} | {subset} | {row["Cv_numerator"]["cells"]} | {values[0]} | {values[1]} | {row["both_positive"]} |')
    lines += ['', '- Floor3.2 is the default-source candidate; floor3.1 is explicitly experimental. Both use radiation-split shifted native jets and the SAME centered floor3.2 baseline. No source below3.1 or invalid four-corner cell is admitted.',
              '- Cv numerator adds monotone analytic radiation; Prho has zero analytic radiation sign term but retains the restored-core cancellation guard. Existing longdouble Hermite/Bernstein margins and adaptive dyadic bounds are reused unchanged.',
              '- A negative coefficient/lower bound is not a witness. Witnesses require finite negative restored-core Cv or Prho at an actual forced-cell query. Depth/node/time exhaustion remains explicitly unresolved.',
              '- Four restored-core queries per new cell compare polynomial+analytic values against the effective assembly/core/radiation/restoration margin. Any failed calibration suppresses positive certificates.',
              '- full_one_branch_own_use is the exact geometric 024 helper predicate, NOT an actual MESA contribution fraction. Experimental cells below3.2 cannot enter that predicate.',
              f'- Run status {result["status"]}; elapsed {result["elapsed_seconds"]:.2f}s; depth {result["limits"]["max_depth"]}, nodes/field/cell {result["limits"]["max_nodes"]}, new-cell cap/candidate {result["limits"]["max_cells"]}. Detailed unresolved cells, witnesses, sample margins, source/API hashes and indices are in JSON/NPZ.']
    (OUTPUT / 'README.md').write_text('\n'.join(lines) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-cells', type=int, default=20000)
    parser.add_argument('--max-depth', type=int, default=3)
    parser.add_argument('--max-nodes', type=int, default=85)
    parser.add_argument('--batch-size', type=int, default=512)
    parser.add_argument('--seconds', type=float, default=180)
    parser.add_argument('--finish-by', help='optional timezone-aware ISO finish fence; default local 20:00 each day')
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    check_limits(args.max_cells, args.max_depth, args.max_nodes, args.seconds)
    if not 1 <= args.batch_size <= 1024:
        raise ValueError('Batch size must be 1..1024')
    now, started = datetime.now(ZoneInfo('US/Pacific')), time.monotonic()
    fence = finish_fence(now, args.finish_by)
    if now >= fence:
        raise RuntimeError('Requested/effective finish fence reached')
    deadline = started + min(args.seconds, (fence - now).total_seconds())
    OUTPUT.mkdir(parents=True, exist_ok=True)
    if (OUTPUT / 'results.json').exists() and not args.overwrite:
        raise ValueError('Own output exists; --overwrite required')
    inputs = [TABLE, BIQUINTIC / 'freeeos_potential.py'] + [BIQUINTIC / 'positivity_bounds' / name for name in ('bounds.py', 'run_bounds.py', 'bound_split.py')]
    inputs += [BIQUINTIC / 'edge_control' / name for name in ('edge_control.py', 'clean_floor.py')]
    before = {str(path.relative_to(WORKSPACE)): checksum(path) for path in inputs}
    if checksum(TABLE) != TABLE_SHA:
        raise ValueError('Pinned source changed')
    _, logq, logt, rows = read_table(TABLE)
    default = build_control(rows, logq, logt, 3.2, True)
    baseline_complete = default['centered_complete'].copy()
    result = {'status': 'running', 'started_local': now.isoformat(), 'command': sys.argv,
              'effective_finish_by_local': fence.isoformat(), 'requested_finish_by': args.finish_by,
              'limits': vars(args), 'candidates': {}, 'input_hashes_before': before,
              'sidecar_source_hashes': {path.name: checksum(path) for path in HERE.iterdir() if path.is_file()},
              'precision': {'bounds': str(np.dtype(DTYPE)), 'core_jets': 'float64', 'longdouble_eps': float(np.finfo(DTYPE).eps),
                            'float64_eps': float(np.finfo(float).eps), 'guard_safety_factor': 4096},
              'fixed_baseline_complete_cells': int(valid_cells(baseline_complete).sum()),
              'guard_policy': 'Reuse existing absolute-Hermite and float64 core margins plus analytic/restoration guards; no guard removed or reduced. Numerical fixed-polynomial scope, not directed-rounding formal proof or source accuracy.'}
    for floor in (3.2, 3.1):
        if time.monotonic() >= deadline:
            result['candidates'][f'{floor:g}'] = {'status': 'not_started_runtime_cap'}
            break
        control = default if floor == 3.2 else build_control(rows, logq, logt, floor, True)
        result['candidates'][f'{floor:g}'] = audit_candidate(control, baseline_complete, logq, logt, rows, floor, args, deadline)
        (OUTPUT / 'results.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    result['status'] = 'complete' if len(result['candidates']) == 2 and all(value['status'] == 'complete' and value['calibration_passed'] for value in result['candidates'].values()) else 'incomplete_or_guard_failure'
    result['elapsed_seconds'] = time.monotonic() - started
    result['input_hashes_after'] = {str(path.relative_to(WORKSPACE)): checksum(path) for path in inputs}
    result['read_only_inputs_unchanged'] = before == result['input_hashes_after']
    result['peak_process_rss_kib_linux'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    (OUTPUT / 'results.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    if all('subsets' in candidate for candidate in result['candidates'].values()):
        write_summary(result)
    print(f'{result["status"]}; {result["elapsed_seconds"]:.2f}s; {OUTPUT / "results.json"}', flush=True)
    if result['status'] != 'complete' or not result['read_only_inputs_unchanged']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
