"""Freeze benchmark provenance and diagnose recorded outputs without rereading any cache."""

from datetime import datetime
import json
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np

from run_benchmark import FIELD_NAMES, OUTPUT, PARTIAL_NAMES, STATE_NAMES, WORKSPACE, parse_report, sha256


def diagnostic_distribution(values):
    values = np.asarray(values, dtype=float)
    return {'count': len(values), 'p50_p90_p99_max': np.quantile(values, [.5, .9, .99, 1]).tolist(),
            'counts_above_diagnostic_thresholds': {str(threshold): int(np.sum(values > threshold))
                                                  for threshold in [1e-8, 1e-6, 1e-4, 1e-2]}}


def response_parent_comparison(real_values, ad_values):
    real_values = np.asarray(real_values, dtype=float)
    ad_values = np.asarray(ad_values, dtype=float)
    if real_values.shape != ad_values.shape or real_values.ndim != 2 or real_values.shape[1] != 43:
        raise ValueError('Expected paired arrays with 43 fields')
    if not np.all(np.isfinite(real_values)) or not np.all(np.isfinite(ad_values)):
        raise ValueError('Nonfinite captured field')
    scales = np.maximum(1, np.maximum(np.abs(real_values), np.abs(ad_values)))
    response_values = np.maximum(np.abs(real_values[:, 10:21]), np.abs(ad_values[:, 10:21]))
    scales[:, 21:32] = np.maximum(scales[:, 21:32], response_values)
    scales[:, 32:43] = np.maximum(scales[:, 32:43], response_values)
    return scales, np.abs(real_values - ad_values) / scales


def field_record(real_values, ad_values, scales, differences, coordinates, query_index, field_index):
    physical_x, physical_t = coordinates[query_index]
    real_value, ad_value = real_values[query_index, field_index], ad_values[query_index, field_index]
    own_scale = max(1, abs(real_value), abs(ad_value))
    record = {'query_index_1based': query_index + 1, 'field': FIELD_NAMES[field_index], 'column_offset_0based': field_index,
              'lnrho': float(physical_x), 'lnT': float(physical_t),
              'logQ': float((physical_x - 2 * physical_t) / np.log(10) + 12), 'logT': float(physical_t / np.log(10)),
              'real': float(real_value), 'ad': float(ad_value), 'absolute_difference': float(abs(real_value - ad_value)),
              'raw_own_column_scale': float(own_scale), 'raw_own_column_relative_difference': float(abs(real_value - ad_value) / own_scale),
              'response_parent_scale': float(scales[query_index, field_index]),
              'response_parent_scaled_difference': float(differences[query_index, field_index])}
    if field_index >= 21:
        parent_index = 10 + (field_index - 21) % 11
        record.update(parent_response=FIELD_NAMES[parent_index],
                      real_parent_response=float(real_values[query_index, parent_index]),
                      ad_parent_response=float(ad_values[query_index, parent_index]))
    if field_index < 10:
        record['potential_context'] = {'real_f': float(real_values[query_index, 0]), 'ad_f': float(ad_values[query_index, 0]),
                                       'real_f_xx': float(real_values[query_index, 3]), 'ad_f_xx': float(ad_values[query_index, 3]),
                                       'note': 'Potential own-column relative differences are cancellation diagnostics, not physical error estimates.'}
    return record


def captured_parent_diagnostics(results, coordinates, maximum):
    manifest_path = OUTPUT / 'preflight/capture_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    trace_path = WORKSPACE / manifest['trace']
    if manifest['cache_sha256'] != results['inputs']['cache_sha256'] or manifest['queries_sha256'] != results['inputs']['queries_sha256']:
        raise ValueError('Capture uses different benchmark inputs')
    if sha256(trace_path) != manifest['trace_sha256']:
        raise ValueError('Recorded capture changed')
    with trace_path.open() as stream:
        trace = np.array([np.fromstring(line.split('trace', 1)[1], sep=' ') for line in stream
                          if line.lstrip().startswith('trace ')])
    if trace.shape != (len(coordinates), 87) or not np.array_equal(trace[:, 0], np.arange(1, len(coordinates) + 1)):
        raise ValueError('Expected every recorded query and 43 paired fields')
    real_values, ad_values = trace[:, 1:44], trace[:, 44:87]
    own_scale = np.maximum(1, np.maximum(abs(real_values), abs(ad_values)))
    raw_difference = abs(real_values - ad_values) / own_scale
    np.testing.assert_allclose(raw_difference.max(axis=0), [maximum[name] for name in FIELD_NAMES], rtol=1e-12, atol=1e-15)
    scales, differences = response_parent_comparison(real_values, ad_values)
    worst_query, worst_field = np.unravel_index(differences[:, 21:43].argmax(), differences[:, 21:43].shape)
    worst_field += 21
    raw_cp_query = int(raw_difference[:, 26].argmax())
    raw_partial_query, raw_partial_field = np.unravel_index(raw_difference[:, :10].argmax(), raw_difference[:, :10].shape)
    diagnostic = {'scaling': 'For each state derivative: abs(real_partial-AD_partial)/max(1,abs(real_state),abs(AD_state),abs(real_partial),abs(AD_partial)). State values use their own value scale. Potential columns retain the raw own-column diagnostic.',
                  'column_offsets_0based': {'potential_partials': [0, 10], 'state_values': [10, 11], 'dlnrho': [21, 11], 'dlnT': [32, 11]},
                  'capture_manifest': manifest,
                  'per_field_maxima_43_columns': dict(zip(FIELD_NAMES, differences.max(axis=0).tolist())),
                  'state_values_per_query_maxima': diagnostic_distribution(differences[:, 10:21].max(axis=1)),
                  'state_jacobians_per_query_maxima': diagnostic_distribution(differences[:, 21:43].max(axis=1)),
                  'worst_response_parent_scaled_jacobian': field_record(real_values, ad_values, scales, differences, coordinates, int(worst_query), int(worst_field)),
                  'raw_Cp_x_relative_worst_with_parent_scale': field_record(real_values, ad_values, scales, differences, coordinates, raw_cp_query, 26),
                  'raw_potential_relative_worst_with_context': field_record(real_values, ad_values, scales, differences, coordinates, int(raw_partial_query), int(raw_partial_field)),
                  'per_field_worst_43_columns': {name: field_record(real_values, ad_values, scales, differences, coordinates, int(differences[:, field].argmax()), field)
                                               for field, name in enumerate(FIELD_NAMES)},
                  'interpretation': 'Raw relative-to-partial tails near cancellation are not project response-parent-scaled errors or physical/source errors. These are numerical diagnostics, not a blanket field-scaled failure gate.'}
    np.savez_compressed(OUTPUT / 'preflight/response_parent_traces.npz', coordinates=coordinates, real=real_values, ad=ad_values,
                        response_parent_scales=scales, response_parent_scaled_differences=differences, raw_own_column_relative_differences=raw_difference,
                        field_names=np.array(FIELD_NAMES))
    (OUTPUT / 'response_parent_scaling.json').write_text(json.dumps(diagnostic, indent=2, allow_nan=False) + '\n')
    return diagnostic


def cache_manifest(inputs, core_path):
    header = inputs['cache_header']
    declared_preprocessing = header[4].split()[0]
    current_core = sha256(core_path)
    return {'benchmarked_cache_path': inputs['cache'], 'benchmarked_cache_sha256': inputs['cache_sha256'],
            'benchmarked_cache_bytes': inputs['cache_bytes'], 'recorded_cache_header': header,
            'declared_source_sha256': header[3], 'declared_preprocessing_sha256': declared_preprocessing,
            'current_core_path': str(core_path.relative_to(WORKSPACE)), 'current_core_sha256': current_core,
            'declared_preprocessing_matches_current_core': declared_preprocessing == current_core,
            'cache_was_reread_by_finalizer': False,
            'interpretation': 'Manifest identifies the exact preexisting cache benchmarked, not a later reexport at the same path. '
                              'Declared preprocessing SHA and current core source SHA are distinct provenance facts. '
                              'The user reports a strict-header-only core change with unchanged solar jets; '
                              'this control does not regenerate or verify numerical equivalence of a later export.'}


def main():
    result_path = OUTPUT / 'results.json'
    results = json.loads(result_path.read_text())
    if results['status'] != 'passed':
        raise ValueError('Only complete benchmark results can be finalized')
    query_path = OUTPUT / 'queries.txt'
    if sha256(query_path) != results['inputs']['queries_sha256']:
        raise ValueError('Archived benchmark query copy changed')
    coordinates = np.loadtxt(query_path, skiprows=1, ndmin=2)
    report = WORKSPACE / results['table']['processes'][0]['report']
    records = parse_report(report.read_text())
    validation = np.array([list(map(float, row)) for row in records['validation']])
    if len(validation) != len(coordinates) or np.any(validation[:, 1:3] != 0):
        raise ValueError('Recorded validation status/count mismatch')
    maximum = results['table']['processes'][0]['max_scaled_real_ad_differences']
    numerical = {'scaling': 'RAW own-column diagnostic: abs(real-AD)/max(1,abs(real),abs(AD)); not project response-parent scaling',
                 'partials_per_query_maxima': diagnostic_distribution(validation[:, 3]),
                 'state_and_jacobian_per_query_maxima': diagnostic_distribution(validation[:, 4]),
                 'max_value_through_second_order_partial_difference': max(maximum[name] for name in PARTIAL_NAMES[:6]),
                 'max_third_order_partial_difference': max(maximum[name] for name in PARTIAL_NAMES[6:]),
                 'max_physical_state_value_difference': max(maximum[name] for name in STATE_NAMES),
                 'raw_max_own_column_state_jacobian_difference': max(maximum[name] for name in FIELD_NAMES[21:]),
                 'interpretation': 'Own-column relative derivative tails are not physical error estimates or a response-scaled blanket failure. '
                                   'Project response-parent diagnostics are reported separately. No source-accuracy or integrated EOS speedup claim follows.'}
    parent_diagnostic = captured_parent_diagnostics(results, coordinates, maximum)
    numerical['project_response_parent'] = parent_diagnostic
    for name, column in [('partials', 3), ('state_jacobians', 4)]:
        row = validation[np.argmax(validation[:, column])]
        query = int(row[0])
        physical_x, physical_t = coordinates[query - 1]
        numerical['worst_' + name] = {'query_index_1based': query, 'lnrho': float(physical_x), 'lnT': float(physical_t),
                                     'logQ': float((physical_x - 2 * physical_t) / np.log(10) + 12),
                                     'logT': float(physical_t / np.log(10)), 'scaled_difference': float(row[column]),
                                     'source': 'recorded untimed validation only; no solar-cache reevaluation'}
    checksum_checks = {}
    for backend in ['table', 'mesa']:
        if backend not in results:
            continue
        for mode, summary in results[backend]['summary'].items():
            lower, upper = summary['checksum_per_pass_min_max']
            relative_spread = abs(upper - lower) / max(1, abs(lower), abs(upper))
            checksum_checks[mode] = {'checksum_per_pass_range': [lower, upper], 'scaled_repeat_spread': relative_spread,
                                     'repeat_consistency_tolerance': 1e-8, 'repeat_consistency_passed': relative_spread < 1e-8}
    if not all(check['repeat_consistency_passed'] for check in checksum_checks.values()):
        raise ValueError('Recorded timing checksum repeat-consistency check failed')
    manifest = cache_manifest(results['inputs'], OUTPUT.parents[2] / 'mesa/eos/eosFreeEOS_builder/biquintic/freeeos_potential.py')
    manifest['recorded_at_local'] = datetime.now(ZoneInfo('US/Pacific')).isoformat()
    manifest['benchmark_started_local'] = results['started_local']
    results['cache_manifest'] = manifest
    results['numerical_diagnostics'] = numerical
    results['timing_checksum_repeat_checks'] = checksum_checks
    results['status_interpretation'] = 'timing/status/cache-safety tests passed; raw and response-parent numerical diagnostics are separate, not a physical-error or blanket numerical-failure gate'
    results['finalizer_provenance'] = {'script': str(Path(__file__).relative_to(WORKSPACE)), 'script_sha256': sha256(__file__),
                                     'recorded_at_local': manifest['recorded_at_local'],
                                     'solar_cache_reread_or_benchmark_rerun': False,
                                     'separate_untimed_capture': 'preflight/capture_manifest.json; cache access finished before recorded-only finalization'}
    (OUTPUT / 'cache_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    result_path.write_text(json.dumps(results, indent=2, allow_nan=False) + '\n')
    results_md = OUTPUT / 'RESULTS.md'
    original = results_md.read_text().split('\n## Final provenance and numerical caution')[0]
    original = '\n'.join('- Raw own-column real/AD maxima: potential partials 1.277e-02; state derivatives 7.670e-02. These divide by each partial itself (floor 1), NOT the parent-response scale and NOT a physical-error or blanket failure estimate; corrected project scaling follows.'
                         if line.startswith('- Real/AD maximum per-field difference') else line for line in original.splitlines())
    lines = ['', '## Final provenance and numerical caution',
             f"- Benchmarked cache SHA-256: `{manifest['benchmarked_cache_sha256']}`.",
             f"- Cache-declared preprocessing SHA-256: `{manifest['declared_preprocessing_sha256']}`.",
             f"- Current core source SHA-256 at finalization: `{manifest['current_core_sha256']}`; equals cache declaration: **{manifest['declared_preprocessing_matches_current_core']}**.",
             '- `cache_manifest.json` freezes the original measured cache identity. A later 008 reexport is not retrospectively the benchmarked cache. Finalization reads recorded 017 artifacts and the current core, never the solar cache.',
             f"- State values reach maximum scaled difference {numerical['max_physical_state_value_difference']:.3e}. Raw own-column third-partial and Jacobian diagnostics retain their tails, but are **not** parent-response-scaled errors or a blanket physical-error/failure claim.",
             '- Timing ratios compare the same requested output contract; no exact production speedup or independent physical-accuracy claim follows.',
             '- All recorded checksum/pass repeat spreads meet the diagnostic 1e-8 relative tolerance; field-level/raw sums remain in JSON.']
    for name, key in [('Potential partials', 'partials_per_query_maxima'), ('States/Jacobians', 'state_and_jacobian_per_query_maxima')]:
        distribution = numerical[key]
        quantiles = ', '.join(f'{value:.3e}' for value in distribution['p50_p90_p99_max'])
        lines.append(f"- RAW own-column {name}: p50/p90/p99/max of per-query maximum relative diagnostic = {quantiles}; counts above diagnostic thresholds {distribution['counts_above_diagnostic_thresholds']}.")
    worst = numerical['worst_state_jacobians']
    lines.append(f"- Worst recorded state/Jacobian query: index {worst['query_index_1based']}, logQ={worst['logQ']:.12g}, logT={worst['logT']:.12g}. This is a numerical-path diagnostic only, not a new active-face physics analysis.")
    lines += ['', '### Project response-parent scaling',
              '- **Lead-review resolution:** the earlier interpretation of raw `Cp_x` .0767 as a project response-scaled failure is withdrawn. Its correct parent-scaled difference is **2.474902236e-9**; the full-corpus parent-scaled Jacobian maximum is **4.898390989e-7**. Raw own-column tails remain cancellation-sensitive diagnostics, not physical-error estimates. Benchmark/capture work is complete; no further benchmark or solar-cache read is pending, and 008 may be reexported.',
              '- State offsets are value=10, dlnrho=21, dlnT=32, eleven fields/block (zero-based). Each state derivative uses max(1, |real response|, |AD response|, |real derivative|, |AD derivative|). Potential columns retain the explicitly labeled own-column diagnostic.',
              f"- The {len(coordinates):,}-query untimed-only capture takes {parent_diagnostic['capture_manifest']['elapsed_seconds']:.2f} seconds including harness build/load/output; original timed binaries/results and cache contents remain unchanged. Finalization uses only recorded traces.",
              '- All 43 columns, scales, differences, and paired raw outputs are archived in `preflight/response_parent_traces.npz`; per-field maxima/worst records are in `response_parent_scaling.json`.']
    for name, key in [('State values', 'state_values_per_query_maxima'), ('State Jacobians', 'state_jacobians_per_query_maxima')]:
        quantiles = ', '.join(f'{value:.3e}' for value in parent_diagnostic[key]['p50_p90_p99_max'])
        lines.append(f'- RESPONSE-PARENT {name} p50/p90/p99/max = {quantiles}.')
    for label, key in [('Former raw-relative Cp_x worst', 'raw_Cp_x_relative_worst_with_parent_scale'),
                       ('Response-parent Jacobian worst', 'worst_response_parent_scaled_jacobian')]:
        record = parent_diagnostic[key]
        lines.append(f"- {label}: query {record['query_index_1based']}, {record['field']}, logQ={record['logQ']:.12g}, logT={record['logT']:.12g}; real/AD derivative={record['real']:.16e}/{record['ad']:.16e}; real/AD {record['parent_response']}={record['real_parent_response']:.16e}/{record['ad_parent_response']:.16e}; denominator={record['response_parent_scale']:.16e}; raw-relative difference={record['raw_own_column_relative_difference']:.9e}; response-parent difference={record['response_parent_scaled_difference']:.9e}.")
    record = parent_diagnostic['raw_potential_relative_worst_with_context']
    context = record['potential_context']
    lines.append(f"- Potential raw-relative worst {record['field']}: real/AD={record['real']:.16e}/{record['ad']:.16e}; parent context f_xx={context['real_f_xx']:.16e}/{context['ad_f_xx']:.16e}, f={context['real_f']:.16e}/{context['ad_f']:.16e}. The {record['raw_own_column_relative_difference']:.9e} own-column ratio is not a physical-error estimate.")
    results_md.write_text(original.rstrip() + '\n' + '\n'.join(lines) + '\n')
    print('Finalized recorded benchmark; solar cache not reread.')
    print('Cache/preprocessing/current-core SHA:', manifest['benchmarked_cache_sha256'],
          manifest['declared_preprocessing_sha256'], manifest['current_core_sha256'])


if __name__ == '__main__':
    main()
