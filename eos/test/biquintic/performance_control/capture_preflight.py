"""Capture existing real/AD outputs once, untimed; never regenerate a table or cache."""

from datetime import datetime
import json
import os
from pathlib import Path
import subprocess
import time
from zoneinfo import ZoneInfo

from run_benchmark import CONTROL, OUTPUT, SDK, STOP_LONG_TASKS, WORKSPACE, sha256


def main():
    started = time.monotonic()
    now = datetime.now(ZoneInfo('US/Pacific'))
    if now >= STOP_LONG_TASKS:
        raise RuntimeError('20:18 PDT stop fence reached')
    deadline = started + min(120, (STOP_LONG_TASKS - now).total_seconds())
    results = json.loads((OUTPUT / 'results.json').read_text())
    cache = WORKSPACE / results['inputs']['cache']
    query = OUTPUT / 'queries.txt'
    if sha256(cache) != results['inputs']['cache_sha256']:
        raise ValueError('Cache changed since benchmark; do not mix exports')
    if sha256(query) != results['inputs']['queries_sha256']:
        raise ValueError('Recorded query copy changed')
    output = OUTPUT / 'preflight'
    output.mkdir(parents=True, exist_ok=True)
    environment = dict(os.environ)
    environment.update(MESASDK_ROOT=str(SDK), OMP_NUM_THREADS='1',
                       LD_LIBRARY_PATH=f'{SDK}/lib:{SDK}/lib64:' + environment.get('LD_LIBRARY_PATH', ''))
    commands = [['make', '-C', str(CONTROL), '-j2', f'BUILD_DIR={output / "build"}', 'table'],
                ['taskset', '-c', str(results['context']['pinned_cpu']), str(output / 'build/table_benchmark'),
                 str(cache), str(query), results['inputs']['cache_header'][3], '3', '.001', '--preflight-only']]
    for command, label in zip(commands, ['build', 'state_trace']):
        with (output / f'{label}.stdout').open('w') as stdout, (output / f'{label}.stderr').open('w') as stderr:
            subprocess.run(command, cwd=output, env=environment, stdout=stdout, stderr=stderr,
                           timeout=max(.001, deadline - time.monotonic()), check=True)
    unchanged = sha256(cache) == results['inputs']['cache_sha256']
    if not unchanged:
        raise ValueError('Cache changed during preflight')
    manifest = {'purpose': 'one untimed real_state/ad_state evaluation of each recorded query; no warm/timed loops',
                'started_local': now.isoformat(), 'elapsed_seconds': time.monotonic() - started,
                'commands': commands, 'cache_sha256': results['inputs']['cache_sha256'],
                'declared_preprocessing_sha256': results['inputs']['cache_header'][4].split()[0],
                'queries_sha256': sha256(query), 'cache_unchanged_before_after': unchanged,
                'trace': str((output / 'state_trace.stdout').relative_to(WORKSPACE)),
                'trace_sha256': sha256(output / 'state_trace.stdout'),
                'binary_sha256': sha256(output / 'build/table_benchmark'),
                'table_driver_sha256': sha256(CONTROL / 'table_benchmark.f90'),
                'support_sha256': sha256(CONTROL / 'benchmark_support.f90'),
                'capture_script_sha256': sha256(__file__),
                'original_timed_benchmark_binary_unchanged': sha256(OUTPUT / 'build/table_benchmark') == results['binary_hashes']['table']}
    (output / 'capture_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Untimed preflight captured in {manifest["elapsed_seconds"]:.2f} s; solar-cache access finished.')


if __name__ == '__main__':
    main()
