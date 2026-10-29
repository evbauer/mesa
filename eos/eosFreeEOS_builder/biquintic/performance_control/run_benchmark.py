"""Bounded, read-only-input throughput control; all products belong to experiment 017."""

import argparse
from datetime import datetime
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
from zoneinfo import ZoneInfo

import numpy as np


CONTROL = Path(__file__).resolve().parent
WORKSPACE = CONTROL.parents[4]
OUTPUT = WORKSPACE / 'freeeos_biquintic_progress/experiments/017_performance_control'
TABLE = CONTROL.parent / 'table_fortran'
MESA = WORKSPACE / 'mesa'
SDK = WORKSPACE / 'mesasdk-x86_64-linux-26.6.1'
MASSFILE = MESA / 'eos/eosFreeEOS_builder/mass_frac_solar.txt'
MASS_SHA = 'bd62922804e073a2a410480994b8a28e85b0bc53e469729d6ed5136ea9da10a6'
STOP_LONG_TASKS = datetime(2026, 10, 28, 20, 18, tzinfo=ZoneInfo('US/Pacific'))
PARTIAL_NAMES = ['f', 'f_x', 'f_t', 'f_xx', 'f_xt', 'f_tt', 'f_xxx', 'f_xxt', 'f_xtt', 'f_ttt']
STATE_NAMES = ['P', 'Pgas', 'E', 'S', 'Cv', 'Cp', 'chiRho', 'chiT', 'gamma1', 'gamma3', 'grad_ad']
FIELD_NAMES = PARTIAL_NAMES + STATE_NAMES + [name + '_x' for name in STATE_NAMES] + [name + '_t' for name in STATE_NAMES]


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def saved_faces(payload, limit, seed):
    with np.load(io.BytesIO(payload), allow_pickle=False) as archive:
        logq = np.concatenate([archive['u_q'], archive['tau_q']])
        logt = np.concatenate([archive['u_t'], archive['tau_t']])
        probes = np.concatenate([archive['u_probe'], archive['tau_probe']])
        directions = np.concatenate([np.zeros(len(archive['u_q']), dtype=int), np.ones(len(archive['tau_q']), dtype=int)])
    if not np.all(np.isfinite(logq)) or not np.all(np.isfinite(logt)):
        raise ValueError('Nonfinite saved face coordinates')
    indices = np.random.default_rng(seed).permutation(len(logq))[:limit or len(logq)]
    coordinates = np.column_stack(((logq[indices] + 2 * logt[indices] - 12) * np.log(10), logt[indices] * np.log(10)))
    metadata = {'available': len(logq), 'selected': len(indices), 'seed': seed,
                'order': 'one deterministic shuffled corpus, repeated identically by every mode/process',
                'u_faces': int(np.sum(directions[indices] == 0)), 'tau_faces': int(np.sum(directions[indices] == 1)),
                'logQ_range': [float(logq[indices].min()), float(logq[indices].max())],
                'logT_range': [float(logt[indices].min()), float(logt[indices].max())],
                'saved_MESA_success_count': int(np.sum(probes[indices, 2] == 0)),
                'saved_FreeEOS_positive_count': int(np.sum(probes[indices, 5] > 0)),
                'saved_Skye_positive_count': int(np.sum(probes[indices, 7] > 0))}
    return coordinates, indices, metadata


def cache_inventory(directory):
    return {path.name: {'size': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns}
            for path in sorted(directory.glob('*.bin'))}


def prepare_mesa_cache(output):
    if sha256(MASSFILE) != MASS_SHA:
        raise ValueError('Scaled-solar pattern differs from experiment 014')
    installed = MESA / 'data/eosDT_data/cache'
    before = cache_inventory(installed)
    if not before:
        raise ValueError('No installed MESA caches; generation is prohibited')
    cache = output / 'runtime_cache/eosDT_cache'
    temporary = output / 'runtime_temp_cache'
    cache.mkdir(parents=True, exist_ok=True)
    temporary.mkdir(parents=True, exist_ok=True)
    for name in before:
        destination = cache / name
        source = installed / name
        if destination.is_symlink():
            if destination.resolve() != source.resolve():
                raise ValueError(f'Unexpected cache link: {destination}')
        elif destination.exists():
            raise ValueError(f'Regenerated private cache is not permitted: {destination}')
        else:
            destination.symlink_to(source)
    return installed, cache, temporary, before


def parse_report(text):
    records = {}
    for line in text.splitlines():
        tokens = line.split()
        if tokens:
            records.setdefault(tokens[0], []).append(tokens[1:])
    return records


def summarize(samples):
    if not samples:
        raise ValueError('No timing samples')
    seconds = np.array([sample['wall_seconds'] / sample['queries'] for sample in samples])
    process_medians = np.array([np.median([sample['wall_seconds'] / sample['queries'] for sample in samples
                                         if sample['process'] == process])
                                for process in sorted({sample['process'] for sample in samples})])
    resamples = np.random.default_rng(20261028).choice(process_medians, (10000, len(process_medians)), replace=True)
    interval = np.quantile(np.median(resamples, axis=1), [.025, .975])
    return {'samples': len(samples), 'processes': len(process_medians),
            'median_seconds_per_query': float(np.median(seconds)),
            'median_queries_per_second': float(1 / np.median(seconds)),
            'p10_p90_seconds_per_query': np.quantile(seconds, [.1, .9]).tolist(),
            'min_max_seconds_per_query': [float(seconds.min()), float(seconds.max())],
            'sample_relative_stddev': float(np.std(seconds, ddof=1) / np.mean(seconds)),
            'approximate_process_bootstrap_95pct_median_seconds_per_query': interval.tolist(),
            'process_median_seconds_per_query': process_medians.tolist(),
            'median_cpu_to_wall_ratio': float(np.median([sample['cpu_seconds'] / sample['wall_seconds'] for sample in samples])),
            'checksum_per_pass_min_max': [min(sample['checksum_per_pass'] for sample in samples),
                                          max(sample['checksum_per_pass'] for sample in samples)],
            'failures': sum(sample['failures'] for sample in samples)}


def analyze_reports(paths, backend):
    samples = {}
    processes = []
    for process, path in enumerate(paths, 1):
        records = parse_report(path.read_text())
        detail = {'process': process, 'report': str(path.relative_to(WORKSPACE)),
                  'phases_seconds': {row[0]: float(row[1]) for row in records['phase_seconds']},
                  'memory_kib': {row[0]: {'rss': int(row[1]), 'peak_rss': int(row[2])} for row in records['memory_kib']},
                  'timer_resolution_seconds': float(records['timer_resolution_seconds'][0][0]),
                  'query_count': int(records['query_count'][0][0]),
                  'warmups': records.get('warmup', []), 'calibration': records.get('calibration', [])}
        for row in records['timing']:
            mode, repeat, passes, count, wall, cpu, checksum, failures = row
            sample = {'process': process, 'repeat': int(repeat), 'passes': int(passes),
                      'queries': int(passes) * int(count), 'wall_seconds': float(wall), 'cpu_seconds': float(cpu),
                      'checksum_per_pass': float(checksum) / int(passes), 'failures': int(failures)}
            if sample['wall_seconds'] <= 0 or sample['failures']:
                raise ValueError('Invalid timing sample')
            samples.setdefault(mode, []).append(sample)
        if backend == 'table':
            detail['metadata'] = {key: records[key][0] for key in ['metadata_grid', 'metadata_composition', 'metadata_preprocessing', 'metadata_source']}
            detail['payload_bytes_excluding_descriptors'] = int(records['payload_bytes'][0][0])
            detail['max_scaled_real_ad_differences'] = dict(zip(FIELD_NAMES, map(float, records['max_scaled_difference'][0])))
            detail['raw_checksums_real'] = dict(zip(FIELD_NAMES, map(float, records['raw_sum_real'][0])))
            detail['raw_checksums_ad'] = dict(zip(FIELD_NAMES, map(float, records['raw_sum_ad'][0])))
            detail['validation_success_count'] = len(records['validation'])
        else:
            detail['output_counts_basic_abundance_species_total'] = list(map(int, records['mesa_output_counts'][0]))
            detail['fraction_names'] = ['HELM', 'OPAL_SCVH', 'FreeEOS', 'PC', 'Skye', 'CMS', 'ideal', 'WDM_AV']
            detail['fraction_means'] = list(map(float, records['fraction_means'][0]))
            detail['fraction_positive_counts'] = list(map(int, records['fraction_positive_counts'][0]))
            detail['raw_checksums'] = list(map(float, records['raw_sum_mesa'][0]))
        processes.append(detail)
    return {'processes': processes, 'samples': samples, 'summary': {mode: summarize(rows) for mode, rows in samples.items()}}


def analytic_state(physical_x, physical_t):
    gas_constant, radiation_constant = 2e7, 1e-6
    density, temperature = np.exp(physical_x), np.exp(physical_t)
    radiation_pressure = radiation_constant * temperature ** 4 / 3
    gas_pressure = density * gas_constant * temperature
    pressure = gas_pressure + radiation_pressure
    cv = 1.5 * gas_constant + 12 * radiation_pressure / (density * temperature)
    chirho = gas_pressure / pressure
    chit = (gas_pressure + 4 * radiation_pressure) / pressure
    gamma3_minus1 = (gas_pressure + 4 * radiation_pressure) / (density * temperature * cv)
    gamma1 = chirho + chit * gamma3_minus1
    return np.array([pressure, gas_pressure, 1.5 * gas_constant * temperature + 3 * radiation_pressure / density,
                     gas_constant * (6.5 + 1.5 * physical_t - physical_x) + 4 * radiation_pressure / (density * temperature),
                     cv, cv + pressure * chit ** 2 / (density * temperature * chirho),
                     chirho, chit, gamma1, 1 + gamma3_minus1, gamma3_minus1 / gamma1])


def check_analytic_smoke(report, coordinate):
    records = parse_report(report.read_text())
    physical_x, physical_t = coordinate
    expected = analytic_state(physical_x, physical_t)
    step = 1e-5
    derivative_x = (analytic_state(physical_x + step, physical_t) - analytic_state(physical_x - step, physical_t)) / (2 * step)
    derivative_t = (analytic_state(physical_x, physical_t + step) - analytic_state(physical_x, physical_t - step)) / (2 * step)
    expected_all = np.concatenate([expected, derivative_x, derivative_t])
    differences = {}
    for mode in ['first_real', 'first_ad']:
        actual = np.array(list(map(float, records[mode][0])))[10:]
        scaled = np.abs(actual - expected_all) / np.maximum(1, np.abs(expected_all))
        differences[mode] = float(scaled.max())
        if differences[mode] > 3e-5:
            raise AssertionError(f'Analytic gas+radiation smoke failed: {mode}: {differences[mode]}')
    return {'passed': True, 'maximum_scaled_state_and_jacobian_errors': differences,
            'tolerance': 3e-5, 'analytic_jacobian_reference': 'central differences, natural-log step 1e-5',
            'fixture': str((TABLE / 'build/ideal_rad.cache').relative_to(WORKSPACE)),
            'fixture_generation': 'none: existing fixture loaded read-only'}


def write_results(output, results):
    (output / 'results.json').write_text(json.dumps(results, indent=2, allow_nan=False) + '\n')
    lines = ['# Experiment 017: bounded in-process performance control', '',
             f"Run: {results['started_local']}; status: {results['status']}; elapsed: {results.get('elapsed_seconds', 0):.2f} s.", '',
             'These are harness throughputs, **not an exact production speedup**. Different outputs and EOS dispatch prevent that claim.', '',
             '| Mode | Median queries/s | Median us/query | p10--p90 us/query | Samples/processes |',
             '|---|---:|---:|---:|---:|']
    for backend in ['table', 'mesa']:
        for mode, summary in results.get(backend, {}).get('summary', {}).items():
            lower, upper = np.array(summary['p10_p90_seconds_per_query']) * 1e6
            lines.append(f"| {mode} | {summary['median_queries_per_second']:.0f} | {summary['median_seconds_per_query'] * 1e6:.3f} | {lower:.3f}--{upper:.3f} | {summary['samples']}/{summary['processes']} |")
    if 'table' in results:
        lines += ['', '## Table checks and memory']
        details = results['table']['processes']
        first = details[0]
        differences = first['max_scaled_real_ad_differences']
        lines += [f"- Identical shuffled corpus: {first['query_count']} queries; all four table modes succeed on every query.",
                  f"- Raw own-column real/AD diagnostic / max(1, |real column|, |AD column|): partials {max(differences[name] for name in PARTIAL_NAMES):.3e}; states/Jacobians {max(differences[name] for name in FIELD_NAMES[10:]):.3e}. This is NOT project response-parent scaling or a physical-error estimate.",
                  f"- Loaded numeric payload: {first['payload_bytes_excluding_descriptors']:,} bytes (nine binary64 jets/node, logical mask, two axes; excludes descriptors/metadata/allocator).",
                  f"- Text cache: {results['inputs']['cache_bytes']:,} bytes. Valid and invalid nodes remain allocated alike."]
        for backend in ['table', 'mesa']:
            if backend not in results:
                continue
            details = results[backend]['processes']
            rss = [detail['memory_kib']['warmed']['rss'] for detail in details]
            peaks = [detail['memory_kib']['finished']['peak_rss'] for detail in details]
            startup = [detail['phases_seconds']['startup_to_validated'] for detail in details]
            lines += [f"- {backend}: warmed RSS {min(rss):,}--{max(rss):,} KiB; peak RSS {min(peaks):,}--{max(peaks):,} KiB; startup through validation {min(startup):.3f}--{max(startup):.3f} s."]
    if 'table_real_ad_throughput_ratios' in results:
        lines += ['', '## Like-output real/AD comparison']
        for name, ratios in results['table_real_ad_throughput_ratios'].items():
            lines.append(f"- {name}: median real/AD throughput ratio {np.median(ratios):.3f}; per-process ratios {', '.join(f'{ratio:.3f}' for ratio in ratios)}.")
    if 'mesa' in results:
        first = results['mesa']['processes'][0]
        lines += ['', '## Installed MESA context',
                  '- Same coordinates and X=.70/Y=.28/Z=.02, normalized scaled-solar metals, 20 representative isotopes, WDM_AV disabled as in 014; otherwise installed defaults.',
                  f"- eosDT_get returns {first['output_counts_basic_abundance_species_total'][0]} basic fields and rho/T derivatives plus 2 abundance derivative fields for 20 isotopes ({first['output_counts_basic_abundance_species_total'][3]} returned scalars). Table evaluation returns 10 potential partials; table state paths additionally return 11 states and 22 rho/T derivatives (43 scalars).",
                  '- MESA computes composition interpolation, component routing/blending, and additional fields; table paths do neither composition interpolation nor dispatch.',
                  f"- Mean blend fractions: {dict(zip(first['fraction_names'], first['fraction_means']))}.",
                  f"- Installed EOS cache inventory unchanged: {results['mesa_cache_safety']['installed_inventory_unchanged']}; private regenerated links: {results['mesa_cache_safety']['replaced_links']}."]
    lines += ['', '## Timing and limitations',
              '- Every mode has a complete untimed warm pass, independent pass-count calibration, rotating mode order, and repeated warmed batches; no per-query text I/O in timed loops.',
              '- Timed work includes status checking, output packing, AD input seeding where applicable, and a normalized sum checksum over every returned scalar. State derivation uses the same metadata radiation constant and derive_state routine on both paths.',
              '- Inputs and MESA linear/log coordinates are prepared outside loops. Checksums are numerical regression evidence, not physical accuracy or an independent EOS validation.',
              '- Startup is fresh-process wall time, not cold disk-cache time. OS caches are not flushed; after-load RSS and first full-domain preflight expose lazy loading separately.',
              '- p10--p90 and between-process medians describe observed scatter. The process-bootstrap interval in JSON is approximate with few processes, not a confidence guarantee against shared-node load, turbo, thermals, or systematic compiler differences.',
              '- Existing table module/kernel objects are linked unchanged. Their original optimization/check flags are not reconstructed; only harness compilation flags are controlled. Installed MESA libraries can have different compiler options.',
              '- The corpus consists of saved valid source-domain face queries, not a stellar-track workload. Exact-face cell ownership is inherited unchanged from the table adapter.',
              '- See README for reproduction, fairness requirements, hard time limits, artifact definitions, and full scope.']
    if 'error' in results:
        lines += ['', 'Error: ' + results['error']]
    (output / 'RESULTS.md').write_text('\n'.join(lines) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--no-mesa', action='store_true')
    parser.add_argument('--smoke', action='store_true', help='load existing analytic fixture; no fixture generation')
    parser.add_argument('--limit', type=int, default=0, help='0 means all saved 014 faces')
    parser.add_argument('--repeats', type=int, default=7)
    parser.add_argument('--process-repeats', type=int, default=3)
    parser.add_argument('--target-seconds', type=float, default=.35)
    parser.add_argument('--time-limit', type=float, default=240, help='whole wrapper wall bound, at most 280 s')
    parser.add_argument('--seed', type=int, default=20261028)
    parser.add_argument('--cpu', type=int, help='default: first allowed CPU; pin both executables identically')
    args = parser.parse_args()
    if not 0 <= args.limit <= 100000 or not 3 <= args.repeats <= 15 or not 1 <= args.process_repeats <= 5:
        parser.error('limit/repeat controls outside bounded ranges')
    if not .001 <= args.target_seconds <= 1 or not 1 <= args.time_limit <= 280:
        parser.error('target/time limit outside bounded ranges')
    started = time.monotonic()
    now = datetime.now(ZoneInfo('US/Pacific'))
    if now >= STOP_LONG_TASKS:
        parser.error('requested 20:18 PDT stop fence reached; do not launch long tasks')
    available = sorted(os.sched_getaffinity(0))
    cpu = args.cpu if args.cpu is not None else available[0]
    if cpu not in available or shutil.which('taskset') is None:
        parser.error('requested CPU is unavailable, or taskset is absent')
    output = OUTPUT / 'smoke' if args.smoke else OUTPUT
    output.mkdir(parents=True, exist_ok=True)
    environment = dict(os.environ)
    environment.update(MESA_DIR=str(MESA), MESASDK_ROOT=str(SDK), OMP_NUM_THREADS='1',
                       OPENBLAS_NUM_THREADS='1', OMP_PROC_BIND='true', OMP_PLACES='cores',
                       LD_LIBRARY_PATH=f'{SDK}/lib:{SDK}/lib64:' + environment.get('LD_LIBRARY_PATH', ''))
    results = {'status': 'running', 'started_local': now.isoformat(), 'arguments': vars(args),
               'command': [sys.executable, *sys.argv],
               'context': {'host': platform.node(), 'platform': platform.platform(), 'allowed_cpus': available,
                           'pinned_cpu': cpu, 'load_average_before': list(os.getloadavg()), 'OMP_NUM_THREADS': 1,
                           'compiler': subprocess.check_output([str(SDK / 'bin/gfortran'), '--version'], env=environment, text=True).splitlines()[0],
                           'cpuinfo': [line.strip() for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith(('model name', 'cpu MHz'))][:2],
                           'memory_info': Path('/proc/meminfo').read_text().splitlines()[:3]},
               'scope': 'only sidecar harness compiled; existing modules, objects, libraries and table/cache inputs are read-only; no exact production speedup claim',
               'harness_hashes': {path.name: sha256(path) for path in sorted(CONTROL.iterdir()) if path.is_file()},
               'linked_table_input_hashes': {str(path.relative_to(WORKSPACE)): sha256(path) for path in
                   [TABLE / 'Makefile', TABLE / 'freeeos_table.f90', TABLE.parent / 'test_fortran/freeeos_biquintic.f90',
                    TABLE / 'build/freeeos_table.o', TABLE / 'build/freeeos_biquintic.o', TABLE / 'build/helm_polynomials.o',
                    TABLE / 'build/freeeos_table.mod']}}
    deadline = min(started + args.time_limit, started + (STOP_LONG_TASKS - now).total_seconds())

    def run(command, label, pinned=False):
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise TimeoutError('bounded wrapper deadline reached')
        actual = ['taskset', '-c', str(cpu), *command] if pinned else command
        with (output / f'{label}.stdout').open('w') as stdout, (output / f'{label}.stderr').open('w') as stderr:
            completed = subprocess.run(actual, cwd=output, env=environment, stdout=stdout, stderr=stderr, timeout=remaining)
        results.setdefault('executed_commands', []).append(actual)
        if completed.returncode:
            raise RuntimeError(f'{label} failed with status {completed.returncode}; see {label}.stderr')
        return output / f'{label}.stdout'

    safety = None
    try:
        cache = TABLE / 'build/ideal_rad.cache' if args.smoke else OUTPUT.parent / '008_standalone_table/solar_split.cache'
        cache_hash = sha256(cache)
        with cache.open() as stream:
            header = [next(stream).strip() for unused in range(6)]
        if header[0] != 'FREEEOS_BIQUINTIC_V1':
            raise ValueError('unexpected existing cache contract')
        if args.smoke:
            probes = np.arange(1, 6)
            natural_t = 5 + .015 * probes
            coordinates = np.column_stack((.14 * probes + 2 * natural_t - 12 * np.log(10), natural_t))
            indices = np.arange(5)
            query_metadata = {'selected': 5, 'purpose': 'existing analytic gas+radiation fixture smoke'}
            args.process_repeats, args.repeats, args.target_seconds = 1, 3, .001
        else:
            archive = OUTPUT.parent / '014_active_domain/active_faces.npz'
            payload = archive.read_bytes()
            coordinates, indices, query_metadata = saved_faces(payload, args.limit, args.seed)
            query_metadata.update(archive=str(archive.relative_to(WORKSPACE)), archive_sha256=hashlib.sha256(payload).hexdigest())
        query_path = output / 'queries.txt'
        with query_path.open('w') as stream:
            stream.write(f'{len(coordinates)}\n')
            np.savetxt(stream, coordinates, fmt='%.17e')
        np.save(output / 'selected_saved_indices.npy', indices, allow_pickle=False)
        results['inputs'] = {'cache': str(cache.relative_to(WORKSPACE)), 'cache_sha256': cache_hash,
                             'cache_bytes': cache.stat().st_size, 'cache_header': header,
                             'queries_sha256': sha256(query_path), 'query_metadata': query_metadata}
        target = 'table' if args.no_mesa or args.smoke else 'all'
        run(['make', '-C', str(CONTROL), '-j2', target], 'build')
        results['binary_hashes'] = {'table': sha256(OUTPUT / 'build/table_benchmark')}
        table_reports = []
        for process in range(1, args.process_repeats + 1):
            table_reports.append(run([str(OUTPUT / 'build/table_benchmark'), str(cache), str(query_path), header[3],
                                      str(args.repeats), str(args.target_seconds)], f'table_{process:02}', pinned=True))
        results['table'] = analyze_reports(table_reports, 'table')
        ratios = {}
        for real_mode, ad_mode in [('real_eval', 'ad_eval'), ('real_state', 'ad_state')]:
            ratios[real_mode + '_vs_' + ad_mode] = (np.array(results['table']['summary'][ad_mode]['process_median_seconds_per_query']) /
                                                  np.array(results['table']['summary'][real_mode]['process_median_seconds_per_query'])).tolist()
        results['table_real_ad_throughput_ratios'] = ratios
        if args.smoke:
            results['analytic_smoke'] = check_analytic_smoke(table_reports[0], coordinates[0])
        if not args.no_mesa and not args.smoke:
            safety = prepare_mesa_cache(output)
            installed, cache_links, temporary, before = safety
            environment.update(MESA_CACHES_DIR=str(cache_links.parent), MESA_TEMP_CACHES_DIR=str(temporary))
            environment.pop('MESA_TEMP_CACHES_DISABLE', None)
            results['binary_hashes']['mesa'] = sha256(OUTPUT / 'build/mesa_benchmark')
            results['mesa_inputs'] = {'inlist': str((CONTROL / 'inlist_benchmark').relative_to(WORKSPACE)),
                                     'inlist_sha256': sha256(CONTROL / 'inlist_benchmark'), 'massfile_sha256': sha256(MASSFILE),
                                     'installed_cache_count': len(before), 'installed_disk_cache_bytes': sum(item['size'] for item in before.values())}
            mesa_reports = []
            for process in range(1, args.process_repeats + 1):
                mesa_reports.append(run([str(OUTPUT / 'build/mesa_benchmark'), str(MESA), str(CONTROL / 'inlist_benchmark'),
                                         str(MASSFILE), str(query_path), str(args.repeats), str(args.target_seconds)], f'mesa_{process:02}', pinned=True))
            results['mesa'] = analyze_reports(mesa_reports, 'mesa')
        if sha256(cache) != cache_hash:
            raise RuntimeError('Read-only standalone cache content changed during run')
        results['readonly_table_inputs_unchanged'] = all(sha256(WORKSPACE / path) == fingerprint
                                                       for path, fingerprint in results['linked_table_input_hashes'].items())
        if not results['readonly_table_inputs_unchanged']:
            raise RuntimeError('linked read-only table inputs changed during run')
        if not args.smoke:
            results['saved_archive_changed_during_run'] = sha256(archive) != query_metadata['archive_sha256']
        results['status'] = 'passed'
    except (OSError, ValueError, RuntimeError, AssertionError, KeyError, IndexError, TimeoutError, subprocess.TimeoutExpired) as error:
        results['status'] = 'failed_or_bounded_stop'
        results['error'] = str(error)
    finally:
        if safety is not None:
            installed, cache_links, temporary, before = safety
            after = cache_inventory(installed)
            replaced = [path.name for path in cache_links.glob('*.bin') if not path.is_symlink()]
            temporary_files = [str(path.relative_to(output)) for path in temporary.rglob('*') if path.is_file()]
            results['mesa_cache_safety'] = {'installed_inventory_before': before, 'installed_inventory_after': after,
                                           'installed_inventory_unchanged': before == after, 'replaced_links': replaced,
                                           'temporary_files': temporary_files}
            if before != after or replaced or temporary_files:
                results['status'] = 'failed_or_bounded_stop'
                results['error'] = 'Cache safety audit failed; no regeneration is permitted'
        results['elapsed_seconds'] = time.monotonic() - started
        results['context']['load_average_after'] = list(os.getloadavg())
        write_results(output, results)
    print(f"{results['status']}: {output / 'RESULTS.md'} ({results['elapsed_seconds']:.2f} s)")
    if results['status'] != 'passed':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
