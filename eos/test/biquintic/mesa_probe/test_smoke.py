"""Bounded streaming, composition, fraction, and invalid-input probe checks."""

import json
import math
import select
import subprocess
import sys
from pathlib import Path


PROBE = Path(__file__).resolve().parent
REPORT = PROBE.parents[4] / 'freeeos_biquintic_progress/experiments/014_active_domain'
HEADER = ('logRho logT ierr frac_HELM frac_OPAL_SCVH frac_FreeEOS frac_PC '
          'frac_Skye frac_CMS frac_ideal frac_WDM_AV Cv dCv_dlnRho dCv_dlnT')


def read_line(process):
    ready, _, _ = select.select([process.stdout], [], [], 60)
    if not ready:
        raise TimeoutError('Probe did not flush a response within 60 seconds')
    line = process.stdout.readline()
    if not line:
        raise RuntimeError('Probe closed stdout unexpectedly')
    return line


def main():
    REPORT.mkdir(parents=True, exist_ok=True)
    cases = [('3.53 6.51', 0), ('2.33 6.765', 0), ('NaN 6', -1002),
             ('1 Infinity', -1002), ('400 6', -1003), ('bad record', -1001),
             ('1 2 3', -1001), ('2.33 6.765', 0)]
    output = []
    with (REPORT / 'smoke.stderr').open('w') as diagnostics:
        process = subprocess.Popen([sys.executable, '-B', str(PROBE / 'run_probe.py')],
                                   stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                   stderr=diagnostics, text=True, bufsize=1)
        try:
            header = read_line(process)
            assert header.strip() == HEADER
            output.append(header)
            for query, expected_status in cases:
                process.stdin.write(query + '\n')
                process.stdin.flush()
                record = read_line(process)
                fields = record.split()
                assert len(fields) == 14
                assert int(fields[2]) == expected_status
                values = [float(value) for value in fields[3:]]
                if expected_status == 0:
                    assert all(math.isfinite(value) for value in values)
                    assert abs(sum(values[:8]) - 1) < 1e-12
                    assert values[7] == 0
                    assert values[8] > 0
                else:
                    assert all(math.isnan(value) for value in values)
                output.append(record)
            assert output[2] == output[-1]
            process.stdin.close()
            assert process.wait(timeout=60) == 0
        finally:
            if process.poll() is None:
                process.kill()
                process.wait()
    (REPORT / 'smoke.tsv').write_text(''.join(output))
    metadata = json.loads((REPORT / 'probe_metadata.json').read_text())
    assert metadata['installed_cache_inventory_unchanged']
    assert not metadata['private_replaced_cache_links']
    assert len(metadata['installed_cache_inventory_before']) == 30
    (REPORT / 'smoke_results.json').write_text(json.dumps({
        'streaming_cases': len(cases), 'passed': True, 'header': HEADER,
        'installed_caches_unchanged': True, 'private_cache_regeneration': False,
        'initialization': 'one process, queries supplied individually before EOF',
        'representative_outputs': [dict(zip(HEADER.split(), record.split())) for record in output[1:3]],
    }, indent=2) + '\n')
    print('PASS: 8 individually flushed queries, default active fractions/Cv, invalid-input guards, recovery, and cache safety.')


if __name__ == '__main__':
    main()
