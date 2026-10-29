"""Stable recorded-output regression; no new fits, source processing, or cache access."""

from datetime import datetime
import json
from pathlib import Path
import sys
import time

import numpy as np

from run_control import (API_DIR, ARCHIVED, FIELDS, HERE, IMPORTANCES, LAMBDAS, LN10, OUTPUT, PATCHES,
                         WORKSPACE, annotate_objective, digest, errors, label_for, paired_improvement,
                         positivity, state_at, thermodynamic_combinations, value_array)
from regularized import roughness_operator


def check_metrics(stored, expected):
    for field in FIELDS:
        for metric in ('median', 'p99', 'max'):
            np.testing.assert_allclose(stored[field][metric], expected[field][metric], rtol=1e-13, atol=1e-15)
        if stored[field]['nonfinite'] != 0 or stored[field]['count'] != expected[field]['count']:
            raise AssertionError('Nonfinite/count regression')


def write_small_summary(results):
    lines = ['# Regional coarse regularized C3 gate — completed 032', '',
             '18 declared coarse cases; printed-table reconciliation/holdouts, not physical truth, measured EOS activity, or chosen runtime.', '',
             '**Observed decision:** the explicit roughness selector resolves the 029 source-fit null modes, but none of these 18 cases improves Cv holdout p99 over its matching archived C4 importance. Individual queries can improve while the tail worsens. No universal weight or runtime representation is selected.', '',
             'The small table shows lambda=.001, the first declared value, **not a holdout-selected winner**. All lambdas and all eleven fields are in RESULTS.md/results.json.', '',
             '| Region | Importance | Retained nodal Cv p99 | C3 heldout Cv p99 / max | Matching C4 Cv p99 / max | Individual Cv improvements |',
             '|---|---:|---:|---:|---:|---:|']
    solve_times, process_times, ranks, stationarity = [], [], [], []
    failed_positivity = 0
    for patch in PATCHES:
        result = results['patches'][patch]
        for importance in IMPORTANCES:
            candidate = result['coarse_fits'][label_for(importance, .001)]
            baseline = result['C4_archived_matching_importance'][f'{importance:g}']['heldout_errors']['Cv']
            cv = candidate['heldout_errors']['Cv']
            improved = candidate['paired_vs_matching_C4']['Cv']
            lines.append(f'| {patch} | {importance:g} | {candidate["retained_nodal_errors"]["Cv"]["p99"]:.3e} | {cv["p99"]:.3e} / {cv["max"]:.3e} | {baseline["p99"]:.3e} / {baseline["max"]:.3e} | {improved["strictly_improved"]}/{improved["count"]} |')
        for candidate in result['coarse_fits'].values():
            solve_times.append(candidate['solve']['elapsed_seconds'])
            process_times.append(candidate['process_wall_seconds'])
            ranks.append(candidate['solve']['rank'])
            stationarity.append(candidate['solve']['augmented_stationarity_relative'])
            failed_positivity += sum(field['nonfinite'] + field['nonpositive_finite']
                                     for subset in ('retained_nodal_positivity', 'heldout_positivity') for field in candidate[subset].values())
    lines += ['', f'- Run: {results["elapsed_seconds"]:.2f} s; 18/18 coarse solves, unknown counts 952/1008/912; solve times {min(solve_times):.3f}–{max(solve_times):.3f} s, process/evaluation times {min(process_times):.3f}–{max(process_times):.3f} s; per-fit subprocess cap 30 s and outer run cap 280 s.',
              f'- All augmented ranks equal unknown counts. Maximum relative stationarity residual {max(stationarity):.3e}. Sampled eleven-field nonfinite/nonpositive counts: {failed_positivity}; this is not a whole-patch certificate.',
              '- Holdouts: 180 ionization, 192 helium_ionization, 170 dense_tail; exact 015 odd/odd order and one printed reference. Independently converted even/even training rows exclude withheld, cross, and outside rows.',
              '- Matching C4 values/coefficients are archived, not refit. BQ numbers in RESULTS.md are archived **global-donor context** from 015, with outside-patch donors and no stored individual holdout BQ array; no fabricated pointwise BQ ranking.',
              '- Degree-five double-knot C3 and normalized integrated normal-fourth roughness are the unchanged 031 API. Source importance=.01/.1 and lambda=.001/.01/.1 are declared sensitivities, not physical truth or an accuracy/positivity guarantee.',
              '- No region is labeled active: actual installed-EOS fractions on these regional queries are not measured. Optional full fits remain deferred; no new long experiment is authorized by these results.',
              '- Reproduce tests, bounded coarse fits, and recorded-only regression using the code README under `mesa/eos/eosFreeEOS_builder/biquintic/regional_c3/`. All generated artifacts are owned by 032; no source/runtime/cache/other-worker write.']
    (OUTPUT / 'README.md').write_text('\n'.join(lines) + '\n')


def main():
    started = time.perf_counter()
    source_paths = list(HERE.glob('*.py')) + [HERE / 'README.md', API_DIR / 'regularized.py', API_DIR / 'spline_control.py']
    frozen = {str(path.relative_to(WORKSPACE)): digest(path) for path in source_paths}
    results = json.loads((OUTPUT / 'results.json').read_text())
    if results['status'] != 'complete' or results['solved_fits'] != 18:
        raise AssertionError('Expected 18 complete coarse solves')
    patch_report, full_rank_count, cv_p99_improvement_count = {}, 0, 0
    for patch in PATCHES:
        result = results['patches'][patch]
        with np.load(OUTPUT / patch / 'paired_queries_and_predictions.npz', allow_pickle=False) as predictions:
            products = {name: predictions[name] for name in predictions.files}
        with np.load(ARCHIVED / patch / 'fits_and_queries.npz', allow_pickle=False) as original:
            np.testing.assert_array_equal(products['withheld_logQ'], original['withheld_logQ'])
            np.testing.assert_array_equal(products['withheld_logT'], original['withheld_logT'])
            np.testing.assert_array_equal(products['local_coarse_common'], original['local_coarse_common'])
            for importance in IMPORTANCES:
                np.testing.assert_array_equal(products[f'C4_importance_{importance:g}_heldout_values'], original[f'thermo_{importance:g}_withheld_values'])
        with np.load(OUTPUT / patch / 'fit_inputs.npz', allow_pickle=False) as training:
            if set(training.files) != {'low', 'coarse_logQ', 'coarse_logT', 'retained_original_iq', 'retained_original_it'}:
                raise AssertionError('Unexpected training payload: potential holdout leak')
            coarse_q, coarse_t = training['coarse_logQ'], training['coarse_logT']
            fit_low = training['low']
            if np.any(training['retained_original_iq'] % 2) or np.any(training['retained_original_it'] % 2):
                raise AssertionError('Training is not original even/even rows')
            if np.intersect1d(training['retained_original_iq'], products['withheld_original_iq']).size or np.intersect1d(training['retained_original_it'], products['withheld_original_it']).size:
                raise AssertionError('Retained/holdout overlap')
        reference = dict(zip(FIELDS, products['printed_reference']))
        node_reference = dict(zip(FIELDS, products['retained_nodal_printed_reference']))
        checked = 0
        for importance in IMPORTANCES:
            baseline = dict(zip(FIELDS, products[f'C4_importance_{importance:g}_heldout_values']))
            for regularization in LAMBDAS:
                label = label_for(importance, regularization)
                candidate = result['coarse_fits'][label]
                solve = candidate['solve']
                annotate_objective(candidate, importance, regularization)
                if solve['rank'] != solve['unknowns'] or solve['unknowns'] >= 2000 or solve['response_weight'] != importance or solve['lambda'] != regularization:
                    raise AssertionError('Rank/size/declared-weight regression')
                full_rank_count += 1
                if solve['augmented_stationarity_relative'] > 1e-8:
                    raise AssertionError('Augmented stationarity regression')
                with np.load(OUTPUT / patch / f'{label}_coefficients.npz', allow_pickle=False) as stored:
                    coefficients = stored['coefficients']
                heldout, _ = state_at(coefficients, coarse_q, coarse_t, products['withheld_logQ'], products['withheld_logT'])
                nodal, native = state_at(coefficients, coarse_q, coarse_t, products['retained_nodal_logQ'], products['retained_nodal_logT'])
                np.testing.assert_allclose(value_array(heldout), products[label + '_heldout_values'], rtol=1e-12, atol=1e-12)
                np.testing.assert_allclose(value_array(nodal), products[label + '_retained_nodal_values'], rtol=1e-12, atol=1e-12)
                check_metrics(candidate['heldout_errors'], errors(heldout, reference, products['withheld_logQ'], products['withheld_logT']))
                check_metrics(candidate['retained_nodal_errors'], errors(nodal, node_reference, products['retained_nodal_logQ'], products['retained_nodal_logT']))
                paired, individual = paired_improvement(heldout, baseline, reference, products['withheld_logQ'], products['withheld_logT'])
                np.testing.assert_array_equal(individual, products[label + '_individual_error_reductions_vs_matching_C4'])
                if paired != candidate['paired_vs_matching_C4']:
                    raise AssertionError('Individual pairing regression')
                if positivity(heldout) != candidate['heldout_positivity'] or positivity(nodal) != candidate['retained_nodal_positivity']:
                    raise AssertionError('Positivity regression')
                fitted, targets = thermodynamic_combinations(native), thermodynamic_combinations(fit_low)
                replayed_data = 0.
                for field, target in targets.items():
                    target = target.ravel()
                    scale = np.maximum(abs(target), solve['relative_source_floors'][field])
                    multiplier = importance if field in ('P_rho', 'P_over_rho_chiT', 'TCv') else 1
                    replayed_data += np.sum((multiplier * (fitted[field] - target) / scale)**2)
                if solve['source_rows'] != 6 * len(coarse_q) * len(coarse_t):
                    raise AssertionError('Source objective mean normalization differs from six node combinations')
                replayed_data /= solve['source_rows']
                rough, _ = roughness_operator(coarse_q * LN10, coarse_t * LN10)
                replayed_roughness = float(np.sum((rough @ (coefficients.ravel(order='F') / solve['potential_scale']))**2))
                np.testing.assert_allclose(replayed_data, solve['data_objective'], rtol=1e-5, atol=1e-16)
                np.testing.assert_allclose(replayed_roughness, solve['roughness_objective'], rtol=1e-5, atol=1e-18)
                candidate['objective_normalization']['recorded_coefficient_replay'] = {
                    'data_objective': float(replayed_data), 'unmultiplied_roughness_objective': replayed_roughness,
                    'total_augmented_objective': float(replayed_data + regularization * replayed_roughness),
                    'normalization_verified': True}
                cv_p99_improvement_count += int(candidate['heldout_errors']['Cv']['p99'] < result['C4_archived_matching_importance'][f'{importance:g}']['heldout_errors']['Cv']['p99'])
                checked += 1
        patch_report[patch] = {'checked_fits': checked, 'checked_holdout_count': result['holdout_count'],
                               'query_identity': True, 'paired_improvement_identity': True, 'coefficient_replay': True}
        (OUTPUT / patch / 'results.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    after = {str(path.relative_to(WORKSPACE)): digest(path) for path in source_paths}
    if frozen != after:
        raise AssertionError('Code/API changed during final stable regression')
    report = {'status': 'passed', 'timestamp_local': datetime.now().astimezone().isoformat(), 'command': sys.argv,
              'elapsed_seconds': time.perf_counter() - started, 'new_fits': 0, 'source_table_reprocessed': False,
              'cache_access': False, 'frozen_source_hashes': frozen, 'code_unchanged_during_regression': True,
              'objective_mean_importance_squared_and_lambda_normalization_verified': True,
              'checked_full_rank_fits': full_rank_count, 'Cv_holdout_p99_improvements_vs_matching_C4': cv_p99_improvement_count,
              'patches': patch_report}
    (OUTPUT / 'validation_results.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    results['final_stable_regression'] = report
    (OUTPUT / 'results.json').write_text(json.dumps(results, indent=2, allow_nan=False) + '\n')
    write_small_summary(results)
    print(f'Stable regression passed: {full_rank_count} fits, no new solves/source/cache access; {report["elapsed_seconds"]:.2f} s.')


if __name__ == '__main__':
    main()
