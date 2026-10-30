"""High-precision endpoint and within-cell finite-difference validation."""

import argparse
import json
import os
import resource
import socket
from decimal import Decimal, localcontext
from pathlib import Path

import numpy as np

from native import NativeHelm, REAL, decimal_cv
from run_control import MESA, REPORT, distribution, dump, environment, run_kernel


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=REPORT/'experiments/048_helm_control')
    args = parser.parse_args()
    if socket.gethostname().split('.')[0] != 'rzwhippet14' or os.environ.get('SLURM_JOB_ID') != '5528075':
        raise RuntimeError('TASK048 placement guard failed')
    resource.setrlimit(resource.RLIMIT_AS, (8*1024**3, 8*1024**3))
    output = args.output.resolve()
    result = json.loads((output/'results.json').read_text())
    table = NativeHelm(MESA/'data/eosDT_data/cache/helm_table.bin',
                       output/'density_axis.txt', output/'temperature_axis.txt')
    ye = result['composition']['Ye']
    precision_rows = []
    for face in result['faces']:
        values = []
        for side in ('left', 'right'):
            reference = decimal_cv(table, face['din'], face['T'], ye, face[side+'_cell'])
            computed = table.cv(face['din'], face['T'], ye, face[side+'_cell'])
            reference_array = np.asarray([str(value) for value in reference], dtype=REAL)
            scale = max(1, *np.abs(face['total_'+side]))
            precision_rows.append({'face_id': face['face_id'], 'side': side,
                'error_over_total_response_scale': float(np.max(np.abs(computed-reference_array))/scale),
                'reference_Cv_Cv_x_Cv_t': [str(value) for value in reference]})
            values.append(reference)
        with localcontext() as context:
            context.prec = 80
            for partial_index, field in ((1, 'Cv_x'), (2, 'Cv_t')):
                reference_jump = values[0][partial_index]-values[1][partial_index]
                declared = face['electron_jump'][field]
                precision_rows[-1][field+'_jump_error_over_scale'] = float(abs(
                    Decimal.from_float(declared['left_minus_right']) - reference_jump) / Decimal.from_float(declared['scale']))
    queries, records = [], []
    steps = (1e-4, 3e-5, 1e-5)
    for face in result['faces']:
        cell = face['right_cell']
        din = float((table.density[cell[0]] + table.density[cell[0]+1])/2)
        temperature = float((table.temperature[cell[1]] + table.temperature[cell[1]+1])/2)
        density = din / ye
        center_index = len(queries)
        queries.append((density, temperature))
        for axis_number, field in ((0, 'Cv_x'), (1, 'Cv_t')):
            for step in steps:
                minus, plus = [density, temperature], [density, temperature]
                minus[axis_number] *= np.exp(-step)
                plus[axis_number] *= np.exp(step)
                records.append({'face_id': face['face_id'], 'group': face['group'], 'field': field,
                                'step': step, 'center': center_index, 'minus': len(queries), 'plus': len(queries)+1})
                queries.extend((minus, plus))
    kernel = run_kernel(output, queries, 'finite_difference', environment(output))
    errors = []
    for record in records:
        center, minus, plus = [kernel[record[key]] for key in ('center', 'minus', 'plus')]
        partial_index = 16 if record['field'] == 'Cv_x' else 17
        numerical = (plus[15]-minus[15])/(2*record['step'])
        scale = max(1, abs(center[15]), abs(center[partial_index]))
        record.update(analytic=float(center[partial_index]), numerical=float(numerical),
                      error_over_response_scale=float(abs(numerical-center[partial_index])/scale))
        errors.append(record['error_over_response_scale'])
    checks = {'endpoint_count': len(precision_rows), 'decimal_precision': 80,
              'decimal_vs_extended_precision': distribution([row['error_over_total_response_scale'] for row in precision_rows]),
              'jump_error_over_response_scale': {field: distribution([row[field+'_jump_error_over_scale']
                  for row in precision_rows if field+'_jump_error_over_scale' in row]) for field in ('Cv_x', 'Cv_t')},
              'finite_difference_count': len(records), 'finite_difference_queries': len(queries),
              'finite_difference_by_step': {str(step): distribution([record['error_over_response_scale']
                  for record in records if record['step'] == step]) for step in steps},
              'finite_difference_by_group': {group: distribution([record['error_over_response_scale']
                  for record in records if record['group'] == group]) for group in ('frozen_FreeEOS_24', 'HELM_context_12')},
              'endpoint_records': precision_rows, 'finite_difference_records': records}
    dump(output/'independent_validation.json', checks)
    print(json.dumps({key: value for key, value in checks.items() if not key.endswith('_records')}, indent=2))
    if checks['decimal_vs_extended_precision']['max'] > 1e-7:
        raise RuntimeError('extended-precision endpoint error exceeds declared validation budget')


if __name__ == '__main__':
    main()
