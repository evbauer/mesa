"""Three regional coarse regularized-C3 gates using the unchanged regularized API."""

import argparse
from datetime import datetime
from importlib.metadata import version
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time
from zoneinfo import ZoneInfo

for variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[variable] = '16'

import numpy as np
from threadpoolctl import threadpool_info, threadpool_limits

from helpers import (BIQUINTIC, FIELDS, HERE, IMPORTANCES, LAMBDAS, PATCHES, WORKSPACE, errors,
                     paired_improvement, positivity, printed_reference, retained_constraints, summary,
                     validate_query_identity, value_array)

API_DIR = BIQUINTIC / 'c3_spline_control'
sys.path.insert(0, str(API_DIR))

from regularized import (LN10, LOW_ORDERS, digest, evaluate_fit, knots_for, solve_regularized, state_at,
                         thermodynamic_combinations, thermodynamic_partials)
from freeeos_potential import read_table, source_constraints

OUTPUT = WORKSPACE / 'freeeos_biquintic_progress/experiments/032_regional_c3'
ARCHIVED = OUTPUT.parent / '015_smooth_validation'
TABLE = WORKSPACE / 'mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data'
TABLE_SHA = 'e6ca7ce27bee8197bb257a8c4886a06366e6951e58958fbed3735c877ee68502'
DEADLINE = datetime(2026, 10, 28, 20, 0, tzinfo=ZoneInfo('US/Pacific'))


def label_for(importance, regularization):
    return f'importance_{importance:g}_lambda_{regularization:g}'


def objective_normalization(solve, importance, regularization):
    return {'source_rows': solve['source_rows'], 'source_mean_divisor': solve['source_rows'],
            'source_operator_and_target_multiplier': 1 / np.sqrt(solve['source_rows']),
            'potential_coefficient_scale': solve['potential_scale'], 'source_response_importance': importance,
            'source_response_importance_squared': importance**2,
            'importance_one_fields': ['E', 'TS', 'P_over_rho'],
            'response_importance_fields': ['P_rho', 'P_over_rho_chiT', 'TCv'],
            'relative_source_floors': solve['relative_source_floors'],
            'lambda': regularization, 'roughness_operator_multiplier': float(np.sqrt(regularization)),
            'data_objective': solve['data_objective'], 'unmultiplied_roughness_objective': solve['roughness_objective'],
            'lambda_times_roughness_objective': regularization * solve['roughness_objective'],
            'total_augmented_objective': solve['data_objective'] + regularization * solve['roughness_objective'],
            'roughness_normalization': solve['roughness']['formula']}


def annotate_objective(candidate, importance, regularization):
    solve = candidate['solve']
    solve.setdefault('API_source_objective_text', solve['source_objective'])
    solve['source_objective'] = f'Mean over {solve["source_rows"]} six-combination/source-node rows of weighted relative residual squared; response importance {importance:g} enters squared.'
    candidate['objective_normalization'] = objective_normalization(solve, importance, regularization)


def worker(patch, importance, regularization):
    if patch not in PATCHES or importance not in IMPORTANCES or regularization not in LAMBDAS:
        raise ValueError('Undeclared patch/importance/lambda')
    output = OUTPUT / patch
    with np.load(output / 'fit_inputs.npz', allow_pickle=False) as archive:
        low, coarse_q, coarse_t = archive['low'], archive['coarse_logQ'], archive['coarse_logT']
    unknowns = (len(knots_for(coarse_q)) - 6) * (len(knots_for(coarse_t)) - 6)
    if unknowns >= 2000:
        raise ValueError('Regional coarse unknown bound is strictly below 2000')
    started = time.perf_counter()
    with threadpool_limits(limits=16):
        coefficients, solve = solve_regularized(low, coarse_q * LN10, coarse_t * LN10, regularization,
                                                response_weight=importance, verify_sparse=False)
    label = label_for(importance, regularization)
    np.savez_compressed(output / f'{label}_coefficients.npz', coefficients=coefficients)
    result = {'status': 'solved', 'solve': solve, 'elapsed_build_solve_write_seconds': time.perf_counter() - started,
              'worker_peak_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              'actual_source_objective': f'Mean weighted relative six-combination residual squared; declared response importance {importance:g} enters squared. sqrt(lambda) multiplies the normalized integrated normal-fourth roughness rows.',
              'API_diagnostic_text_context': 'The API source_objective wording may name its .1 default. The actual response_weight argument and recorded solve.response_weight control this fit.',
              'training_inputs_sha256': digest(output / 'fit_inputs.npz'), 'regularized_API_sha256_at_worker': digest(API_DIR / 'regularized.py')}
    annotate_objective(result, importance, regularization)
    (output / f'{label}_solve.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(f'{patch} {label}: rank {solve["rank"]}/{unknowns}; solve {solve["elapsed_seconds"]:.3f} s', flush=True)


def checkpoint(results):
    (OUTPUT / 'results.json').write_text(json.dumps(results, indent=2, allow_nan=False) + '\n')


def make_table(results):
    lines = ['# Experiment 032: regional coarse regularized C3', '',
             'Printed-source reconciliation and paired holdout fidelity only; not physical truth, measured EOS activity, or a selected runtime.', '',
             '| Region | Candidate | Cv p99 / max | Pgas p99 / max | E p99 / max | S p99 / max | Cv improved vs matching C4 |',
             '|---|---|---:|---:|---:|---:|---:|']
    for patch, result in results['patches'].items():
        baseline = result['BQ_global_donor_archived_context']['archived_withheld_departure']['all']['relative_departure']
        lines.append('| ' + patch + ' | archived BQ global donors | ' + ' | '.join(f'{baseline[field]["p99"]:.3e} / {baseline[field]["max"]:.3e}' for field in ['Cv', 'Pgas', 'E', 'S']) + ' | not individually archived |')
        for importance in IMPORTANCES:
            baseline = result['C4_archived_matching_importance'][f'{importance:g}']['heldout_errors']
            lines.append('| ' + patch + f' | archived C4 importance {importance:g} | ' + ' | '.join(f'{baseline[field]["p99"]:.3e} / {baseline[field]["max"]:.3e}' for field in ['Cv', 'Pgas', 'E', 'S']) + ' | baseline |')
            for regularization in LAMBDAS:
                label = label_for(importance, regularization)
                candidate = result['coarse_fits'].get(label, {'status': 'not_started'})
                if candidate['status'] != 'solved':
                    lines.append(f'| {patch} | {label} | {candidate["status"]} | — | — | — | — |')
                    continue
                errors_by_field = candidate['heldout_errors']
                improvement = candidate['paired_vs_matching_C4']['Cv']
                lines.append('| ' + patch + f' | C3 w={importance:g}, lambda={regularization:g} | ' + ' | '.join(f'{errors_by_field[field]["p99"]:.3e} / {errors_by_field[field]["max"]:.3e}' for field in ['Cv', 'Pgas', 'E', 'S']) + f' | {improvement["strictly_improved"]}/{improvement["count"]} |')
    lines += ['', '- C4 comparisons use identical importance, archived holdout order, and one printed reference. No C4 refit or target-based parameter selection.',
              '- BQ is the 015 archived global-donor control; it uses clean source donors outside the local patch. Its summaries are context, not a same-training-data win or individual paired ranking.',
              '- Explicit integrated normal-fourth roughness resolves 029 source-fit null modes through a declared augmented objective; full-rank QR and stationarity are recorded. It is not physical truth or a universal winning weight.',
              '- No regional geometry is labeled active: no actual installed-EOS fraction measurements are supplied for these three holdout arrays.',
              '- All eleven-field nodal/holdout errors, individual paired reductions, positivity, runtime caps, fingerprints, and ranks are in JSON/NPZ. Positivity is sampled, not a whole-patch certificate.',
              f'- Whole run status: {results["status"]}; elapsed {results.get("elapsed_seconds", 0):.2f} s; optional full fits deferred (coarse-only gate).']
    (OUTPUT / 'RESULTS.md').write_text('\n'.join(lines) + '\n')


def run(args):
    started = time.perf_counter()
    now = datetime.now(ZoneInfo('US/Pacific'))
    if now >= DEADLINE:
        raise RuntimeError('20:00 PDT finish fence reached; no new fit starts')
    if not 1 <= args.max_seconds <= 280 or not 1 <= args.fit_cap_seconds <= 30:
        raise ValueError('Whole-run cap must be 1..280 s and per-fit process cap 1..30 s')
    OUTPUT.mkdir(parents=True, exist_ok=True)
    if (OUTPUT / 'results.json').exists() and not args.overwrite:
        raise ValueError('Own results exist; --overwrite required')
    deadline = started + min(args.max_seconds, (DEADLINE - now).total_seconds())
    if digest(TABLE) != TABLE_SHA:
        raise ValueError('Pinned printed source table changed')
    reference_result = json.loads((ARCHIVED / 'results.json').read_text())
    _, logq, logt, rows = read_table(TABLE)
    tracked = [TABLE, ARCHIVED / 'results.json', API_DIR / 'regularized.py', API_DIR / 'spline_control.py', BIQUINTIC / 'freeeos_potential.py',
               OUTPUT.parent / '031_regularized_c3/results.json', OUTPUT.parent / '029_c3_spline/results.json']
    tracked += [ARCHIVED / patch / 'fits_and_queries.npz' for patch in PATCHES]
    before = {str(path.relative_to(WORKSPACE)): digest(path) for path in tracked}
    results = {'status': 'running', 'started_local': now.isoformat(), 'command': sys.argv, 'patches': {},
               'claim': 'Regional table-only roughness-selected C3 feasibility; no physical truth, universal winning weight, measured activity, or integrated runtime claim.',
               'versions': {'python': sys.version, 'numpy': np.__version__, 'scipy': version('scipy')},
               'read_only_input_hashes_before': before,
               'sidecar_source_hashes': {path.name: digest(path) for path in HERE.iterdir() if path.is_file()},
               'bounds': {'max_new_fits': 18, 'resolution': 'coarse_only', 'max_unknowns_exclusive': 2000,
                          'fit_process_seconds': args.fit_cap_seconds, 'whole_run_seconds': args.max_seconds,
                          'deadline_local': DEADLINE.isoformat(), 'threads': 16, 'importance': list(IMPORTANCES), 'lambda': list(LAMBDAS),
                          'optional_full_fits': 'deferred: this task executes only the first coarse gate'},
               'objective': 'Same six weighted relative E,TS,P/rho,P_rho,P/rho*chiT,T*Cv source combinations as 015 plus lambda normalized integrated normal-fourth potential roughness. Source importance enters squared, sqrt(lambda) multiplies roughness rows.',
               'zero_leak_policy': 'Fit worker receives only independently converted even/even retained source rows, never printed holdout targets, C4 predictions, BQ donors, or outside-patch values.',
               'gauge_radiation': 'Same unsplit total specific potential and current CRAD total-pressure convention as 015; default source floor 3.2; no radiation/gauge refit.',
               'core_provenance_context': {'015_core_declared': reference_result['provenance']['files'].get(str(BIQUINTIC / 'freeeos_potential.py')),
                                           'current_core': before[str((BIQUINTIC / 'freeeos_potential.py').relative_to(WORKSPACE))],
                                           'meaning': 'Both source snapshots recorded; retained pointwise constraints use the current strict header reader and unchanged default 3.2 floor.'}}
    attempted, solved = 0, 0
    for patch in PATCHES:
        archived_patch = reference_result['patches'][patch]
        patch_output = OUTPUT / patch
        patch_output.mkdir(parents=True, exist_ok=True)
        with np.load(ARCHIVED / patch / 'fits_and_queries.npz', allow_pickle=False) as source_archive:
            archive = {key: source_archive[key] for key in source_archive.files}
        selection = archived_patch['selected_indices']
        retained_q, retained_t, withheld_q, withheld_t = validate_query_identity(archive, logq, logt, selection)
        start_q, start_t, count_q, count_t = selection
        _, valid, _ = source_constraints(rows[start_q:start_q + count_q, start_t:start_t + count_t],
                                         archive['patch_logQ'], archive['patch_logT'], minimum_logT=3.2)
        if not valid.all():
            raise ValueError('015 patch no longer entirely native valid')
        low, coarse_q, coarse_t, _, _ = retained_constraints(rows, logq, logt, selection)
        np.savez_compressed(patch_output / 'fit_inputs.npz', low=low, coarse_logQ=coarse_q, coarse_logT=coarse_t,
                            retained_original_iq=retained_q, retained_original_it=retained_t)
        node_q, node_t = np.meshgrid(coarse_q, coarse_t, indexing='ij')
        node_reference = printed_reference(rows[retained_q[:, None], retained_t[None, :]])
        node_reference = {field: values.ravel() for field, values in node_reference.items()}
        reference = printed_reference(rows[withheld_q, withheld_t])
        sample_q, sample_t = archive['withheld_logQ'], archive['withheld_logT']
        products = {'fields': np.array(FIELDS), 'withheld_logQ': sample_q, 'withheld_logT': sample_t,
                    'withheld_original_iq': withheld_q, 'withheld_original_it': withheld_t,
                    'local_coarse_common': archive['local_coarse_common'], 'printed_reference': value_array(reference),
                    'retained_nodal_printed_reference': value_array(node_reference), 'retained_nodal_logQ': node_q.ravel(),
                    'retained_nodal_logT': node_t.ravel()}
        result = {'selection_from015': selection, 'bounds_from015': archived_patch['actual_bounds'],
                  'geometry_activity': 'not measured; regional label only', 'holdout_count': len(sample_q),
                  'retained_node_count': node_q.size, 'coarse_unknowns': (2 * len(coarse_q) + 2) * (2 * len(coarse_t) + 2),
                  'query_identity': 'Exact archived order and original odd/odd source indices verified; no query regeneration.',
                  'BQ_global_donor_archived_context': archived_patch['coarse_biquintic_global_donor_control'],
                  'BQ_caveat': archived_patch['controls_caveat'] + '; individual global-donor holdout values were not stored in 015, so no new pointwise BQ ranking or donor regeneration.',
                  'local_donor_subset_count_archived': int(archive['local_coarse_common'].sum()),
                  'C4_archived_matching_importance': {}, 'coarse_fits': {}}
        results['patches'][patch] = result
        baselines = {}
        for importance in IMPORTANCES:
            archive_label = f'thermo_{importance:g}'
            baseline = dict(zip(FIELDS, archive[archive_label + '_withheld_values']))
            baselines[importance] = baseline
            c4_native = evaluate_fit(archive[archive_label + '_coarse_coefficients'], coarse_q * LN10, coarse_t * LN10,
                                     node_q.ravel() * LN10, node_t.ravel() * LN10, multiplicity=1)
            c4_nodal = thermodynamic_partials(c4_native, node_q.ravel(), node_t.ravel())
            result['C4_archived_matching_importance'][f'{importance:g}'] = {
                'importance': importance, 'source': '015 saved coarse coefficients and withheld state values; no refit',
                'archived_solve': archived_patch['candidates'][archive_label]['coarse']['solve'],
                'heldout_errors': errors(baseline, reference, sample_q, sample_t, withheld_q, withheld_t),
                'retained_nodal_errors': errors(c4_nodal, node_reference, node_q.ravel(), node_t.ravel()),
                'heldout_positivity': positivity(baseline), 'retained_nodal_positivity': positivity(c4_nodal)}
            products[f'C4_importance_{importance:g}_heldout_values'] = value_array(baseline)
            products[f'C4_importance_{importance:g}_retained_nodal_values'] = value_array(c4_nodal)
        for importance in IMPORTANCES:
            for regularization in LAMBDAS:
                label = label_for(importance, regularization)
                remaining = deadline - time.perf_counter()
                if remaining <= 0:
                    result['coarse_fits'][label] = {'status': 'not_started_runtime_cap'}
                    continue
                attempted += 1
                command = [sys.executable, '-B', str(HERE / 'run_control.py'), '--worker', patch, str(importance), str(regularization)]
                cap = min(args.fit_cap_seconds, remaining)
                fit_started = time.perf_counter()
                try:
                    with (patch_output / f'{label}.stdout').open('w') as stdout, (patch_output / f'{label}.stderr').open('w') as stderr:
                        subprocess.run(command, cwd=OUTPUT, stdout=stdout, stderr=stderr, timeout=cap, check=True)
                    candidate = json.loads((patch_output / f'{label}_solve.json').read_text())
                    with np.load(patch_output / f'{label}_coefficients.npz', allow_pickle=False) as coefficient_archive:
                        coefficients = coefficient_archive['coefficients']
                    nodal, native = state_at(coefficients, coarse_q, coarse_t, node_q.ravel(), node_t.ravel())
                    heldout, _ = state_at(coefficients, coarse_q, coarse_t, sample_q, sample_t)
                    candidate.update({'process_wall_seconds': time.perf_counter() - fit_started, 'process_cap_seconds': cap,
                        'retained_nodal_errors': errors(nodal, node_reference, node_q.ravel(), node_t.ravel()),
                        'heldout_errors': errors(heldout, reference, sample_q, sample_t, withheld_q, withheld_t),
                        'retained_nodal_positivity': positivity(nodal), 'heldout_positivity': positivity(heldout)})
                    combinations, target_combinations = thermodynamic_combinations(native), thermodynamic_combinations(low)
                    candidate['six_combination_retained_soft_reconciliation'] = {field: summary((combinations[field] - target.ravel()) /
                        np.maximum(abs(target.ravel()), .01 * np.median(abs(target)))) for field, target in target_combinations.items()}
                    paired, individual = paired_improvement(heldout, baselines[importance], reference, sample_q, sample_t)
                    candidate['paired_vs_matching_C4'] = paired
                    products[label + '_heldout_values'] = value_array(heldout)
                    products[label + '_retained_nodal_values'] = value_array(nodal)
                    products[label + '_individual_error_reductions_vs_matching_C4'] = individual
                    solved += 1
                    print(f'{patch} {label}: solved; Cv p99 {candidate["heldout_errors"]["Cv"]["p99"]:.3e}', flush=True)
                except (subprocess.TimeoutExpired, subprocess.CalledProcessError, ValueError, RuntimeError) as error:
                    candidate = {'status': 'failed_or_runtime_cap', 'reason': str(error), 'process_cap_seconds': cap,
                                 'process_wall_seconds': time.perf_counter() - fit_started}
                    print(f'{patch} {label}: {candidate["status"]}: {error}', flush=True)
                result['coarse_fits'][label] = candidate
                results['elapsed_seconds'] = time.perf_counter() - started
                checkpoint(results)
        np.savez_compressed(patch_output / 'paired_queries_and_predictions.npz', **products)
        (patch_output / 'results.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    results.update(status='complete' if solved == 18 else 'incomplete_bounded', attempted_fits=attempted, solved_fits=solved,
                   elapsed_seconds=time.perf_counter() - started, peak_parent_rss_kib_linux=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                   threadpools=threadpool_info())
    after = {str(path.relative_to(WORKSPACE)): digest(path) for path in tracked}
    results['read_only_input_hashes_after'] = after
    results['read_only_inputs_unchanged'] = {path: before[path] == after[path] for path in before}
    checkpoint(results)
    make_table(results)
    print(f'{solved}/18 solved; whole run {results["elapsed_seconds"]:.2f} s; {OUTPUT / "RESULTS.md"}')
    if solved != 18:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--max-seconds', type=float, default=240)
    parser.add_argument('--fit-cap-seconds', type=float, default=30)
    parser.add_argument('--worker', nargs=3, metavar=('PATCH', 'IMPORTANCE', 'LAMBDA'))
    args = parser.parse_args()
    if args.worker:
        worker(args.worker[0], float(args.worker[1]), float(args.worker[2]))
    else:
        run(args)


if __name__ == '__main__':
    main()
