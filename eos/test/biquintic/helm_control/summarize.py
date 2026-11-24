"""Produce compact review artifacts from the bounded HELM results."""

import argparse
import csv
import json
import os
import resource
import socket
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from run_control import REPORT, dump


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=REPORT/'experiments/048_helm_control')
    args = parser.parse_args()
    if socket.gethostname().split('.')[0] != 'rzwhippet14' or os.environ.get('SLURM_JOB_ID') != '5528075':
        raise RuntimeError('TASK048 placement guard failed')
    resource.setrlimit(resource.RLIMIT_AS, (8*1024**3, 8*1024**3))
    output = args.output.resolve()
    result = json.loads((output/'results.json').read_text())
    seeds = [seed for seed in result['seeds'] if seed['group']=='frozen_FreeEOS_24']
    faces = result['faces']
    case_rows = []
    for seed in seeds:
        neighbors = [face for face in faces if face['seed_id']==seed['case_id'] and face['direction']=='T']
        values = [abs(face['total_jump']['Cv_t']['signed_scaled']) for face in neighbors]
        case_rows.append({'case_id': seed['case_id'],
                          'FreeEOS_J_Cv_t': abs(seed['FreeEOS_jump']['Cv_t']['signed_scaled']),
                          'nearby_HELM_T_min_J_Cv_t': min(values), 'nearby_HELM_T_max_J_Cv_t': max(values),
                          'HELM_to_FreeEOS_Cv_at_exact_seed': seed['HELM_at_exact_seed']['Cv']/seed['split_left']['Cv']})
    with (output/'case_comparison.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(case_rows[0]))
        writer.writeheader()
        writer.writerows(case_rows)
    face_rows = []
    for face in faces:
        row = {key: face[key] for key in ('face_id', 'group', 'direction', 'logRho', 'logT')}
        row.update(frac_HELM=face['default_MESA']['frac_HELM'], frac_FreeEOS=face['default_MESA']['frac_FreeEOS'])
        for component in ('electron', 'total'):
            for side in ('left', 'right'):
                row.update({component+'_'+side+'_'+field: value for field, value in
                            zip(('Cv', 'Cv_x', 'Cv_t'), face[component+'_'+side])})
            for field in ('Cv_x', 'Cv_t'):
                row[component+'_J_'+field] = face[component+'_jump'][field]['signed_scaled']
        face_rows.append(row)
    with (output/'face_traces.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(face_rows[0]))
        writer.writeheader()
        writer.writerows(face_rows)
    figure, axes = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [1.25, 1]})
    positions = np.arange(len(case_rows))
    axes[0].scatter(positions, [row['FreeEOS_J_Cv_t'] for row in case_rows], color='black', s=22,
                    label='Frozen FreeEOS face (one per case)', zorder=3)
    for position, row in zip(positions, case_rows):
        axes[0].plot([position, position], [row['nearby_HELM_T_min_J_Cv_t'], row['nearby_HELM_T_max_J_Cv_t']],
                     color='tab:blue', marker='_', markersize=8, linewidth=2,
                     label='Two neighboring native HELM T faces' if position == 0 else None)
    axes[0].set_yscale('log')
    axes[0].set_ylabel(r'$J(C_{V,\ln T})$, total EOS')
    axes[0].set_xticks(positions)
    axes[0].set_xticklabels([row['case_id'].replace('ordinary', 'ord') for row in case_rows], rotation=75, ha='right', fontsize=7)
    axes[0].legend(fontsize=8, loc='upper left', bbox_to_anchor=(.01, .78))
    axes[0].set_title('Selected FreeEOS faces and nearby HELM faces have different coordinates')
    for group, label, color in (('frozen_FreeEOS_24', 'Near 24 FreeEOS cases', 'tab:blue'),
                                ('HELM_context_12', '12 hotter context states', 'tab:orange')):
        epsilons = sorted({check['epsilon_fraction_of_linear_cell'] for face in faces for check in face['near_checks']})
        medians, maxima = [], []
        for epsilon in epsilons:
            values = [check['kernel_total_limit_error_over_scale'] for face in faces if face['group']==group
                      for check in face['near_checks'] if check['epsilon_fraction_of_linear_cell']==epsilon]
            medians.append(np.median(values))
            maxima.append(max(values))
        axes[1].loglog(epsilons, medians, 'o-', color=color, label=label+' median')
        axes[1].loglog(epsilons, maxima, '--', color=color, label=label+' maximum')
    axes[1].invert_xaxis()
    axes[1].set_xlabel('Distance from face / adjacent linear cell width')
    axes[1].set_ylabel('Kernel-to-exact-limit error / response scale')
    axes[1].legend(fontsize=8, loc='lower left')
    axes[1].set_title('Finite-distance convergence reaches a floating-point floor')
    for axis in axes:
        axis.grid(alpha=.2)
    figure.tight_layout()
    figure.savefig(output/'helm_control.pdf')
    figure.savefig(output/'helm_control.png', dpi=180)
    plt.close(figure)
    pure = [face for face in faces if face['default_MESA']['frac_HELM'] >= 1-1e-12]
    dump(output/'review_summary.json', {'case_comparison': case_rows,
        'pure_HELM_face_ids': [face['face_id'] for face in pure],
        'pure_HELM_T_scaled_jumps': [face['total_jump']['Cv_t']['signed_scaled'] for face in pure if face['direction']=='T'],
        'summaries': result['summaries']})


if __name__ == '__main__':
    main()
