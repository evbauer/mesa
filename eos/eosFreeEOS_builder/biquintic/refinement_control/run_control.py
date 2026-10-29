"""Same-query archived-native refinement audit; no formal order or source-accuracy claim."""

import argparse
import csv
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np

from refinement import (
    BIQUINTIC, FIELDS, STRIDES, build_grids, common_queries, error_summary,
    evaluate_candidate, paired_ratios, table_reference,
)
from freeeos_potential import read_table

WORKSPACE = BIQUINTIC.parents[3]
OUTPUT = WORKSPACE / 'freeeos_biquintic_progress/experiments/021_refinement_control'
TABLE = WORKSPACE / 'mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data'
REGION_SOURCE = WORKSPACE / 'freeeos_biquintic_progress/experiments/015_smooth_validation/results.json'
TABLE_SHA = 'e6ca7ce27bee8197bb257a8c4886a06366e6951e58958fbed3735c877ee68502'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def regions(sample_q, sample_t):
    saved = json.loads(REGION_SOURCE.read_text())
    masks, bounds = {}, {}
    for name, patch in saved['patches'].items():
        box = patch['actual_bounds']
        bounds[name] = box
        masks[name] = ((sample_q >= box['logQ'][0] - 1e-12) & (sample_q <= box['logQ'][1] + 1e-12)
                       & (sample_t >= box['logT'][0] - 1e-12) & (sample_t <= box['logT'][1] + 1e-12))
    bounds['highT_geometry'] = {'logQ': [-8, 5], 'logT': [7.0, 7.8]}
    masks['highT_geometry'] = sample_t >= 7.0
    return masks, bounds


def plot_pairs(errors):
    os.environ['MPLCONFIGDIR'] = str(OUTPUT / '.mplconfig')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(2, 2, figsize=(9, 8), constrained_layout=True)
    for row, (coarse, fine) in enumerate([(8, 4), (4, 2)]):
        for column, field in enumerate(['Cv', 'chiRho']):
            axis = axes[row, column]
            coarse_error = abs(errors[f'BQ_split_stride{coarse}'][field])
            fine_error = abs(errors[f'BQ_split_stride{fine}'][field])
            valid = np.isfinite(coarse_error) & np.isfinite(fine_error)
            display_coarse, display_fine = np.maximum(coarse_error[valid], 1e-16), np.maximum(fine_error[valid], 1e-16)
            axis.scatter(display_coarse, display_fine, s=7, alpha=.25, rasterized=True)
            lower = min(display_coarse.min(), display_fine.min())
            upper = max(display_coarse.max(), display_fine.max())
            axis.plot([lower, upper], [lower, upper], color='black', linewidth=.8)
            axis.set_xscale('log')
            axis.set_yscale('log')
            axis.set_xlabel(f'|archived-field relative error|, stride {coarse}')
            axis.set_ylabel(f'|archived-field relative error|, stride {fine}')
            axis.set_title(f'{field}: same individual targets, {coarse}→{fine}')
            axis.grid(alpha=.2)
    figure.suptitle('BQ reference split: paired native holdouts, fixed sample at every resolution\nBelow diagonal = refinement improves; display-only floor 1e-16; not source accuracy')
    figure.savefig(OUTPUT / 'paired_refinement.png', dpi=160)
    plt.close(figure)


def run(args):
    if digest(TABLE) != TABLE_SHA:
        raise ValueError('Unexpected source-table SHA-256')
    _, logq, logt, rows = read_table(TABLE)
    grids = build_grids(rows, logq, logt)
    density_index, temperature_index, common, coverage = common_queries(rows, logq, logt, grids)
    choices = np.flatnonzero(common)
    if not choices.size:
        raise ValueError('No common heldouts across all retained grids/candidates')
    if choices.size > args.max_queries:
        choices = np.sort(np.random.default_rng(20261028).choice(choices, args.max_queries, replace=False))
    selected_iq, selected_it = density_index[choices], temperature_index[choices]
    assert np.all(selected_iq % 2 == 1) and np.all(selected_it % 2 == 1)
    sample_q, sample_t = logq[selected_iq], logt[selected_it]
    reference = table_reference(rows[selected_iq, selected_it])
    region_masks, region_bounds = regions(sample_q, sample_t)
    errors, results, products = {}, {}, {'iq': selected_iq, 'it': selected_it, 'logQ': sample_q, 'logT': sample_t,
                                        'printed_logRho': rows[selected_iq, selected_it, 18]}
    for field in FIELDS:
        products['reference_' + field] = reference[field]
    for stride, grid in grids.items():
        for candidate_name in grid['candidates']:
            label = f'{candidate_name}_stride{stride}'
            state = evaluate_candidate(candidate_name, grid, sample_q, sample_t)
            errors[label] = {field: (state[field] - reference[field]) / abs(reference[field]) for field in FIELDS}
            positivity = {}
            for field in FIELDS:
                failure = ~np.isfinite(state[field]) | (state[field] <= 0)
                positivity[field] = {'count': int(failure.sum()), 'examples_saved_indices':
                                     np.column_stack((selected_iq[failure], selected_it[failure]))[:8].tolist()}
            results[label] = {
                'queries': int(len(sample_q)), 'coarse_shape': list(grid['rows'].shape[:2]),
                'spacing_dex': [float(grid['q'][1] - grid['q'][0]), float(grid['t'][1] - grid['t'][0])],
                'relative_archived_field_errors': {field: error_summary(error) for field, error in errors[label].items()},
                'positivity_or_positive_response': positivity,
                'regions': {name: {'queries': int(mask.sum()), 'errors': {field: error_summary(error[mask])
                             for field, error in errors[label].items()}} for name, mask in region_masks.items()},
            }
            for field in FIELDS:
                products[label + '_' + field + '_relative_error'] = errors[label][field]
    ratio_results = {}
    for candidate_name in grids[2]['candidates']:
        for coarse, fine in [(8, 4), (4, 2)]:
            label = f'{candidate_name}_{coarse}_over_{fine}'
            ratio_results[label] = {}
            for field in FIELDS:
                ratios, stats = paired_ratios(errors[f'{candidate_name}_stride{coarse}'][field],
                                               errors[f'{candidate_name}_stride{fine}'][field])
                ratio_results[label][field] = stats
                products[label + '_' + field + '_individual_ratio'] = ratios
    output = {
        'claim': 'Common native archived-field holdout/refinement evidence; NOT independent source accuracy or formal convergence order.',
        'table_sha256': TABLE_SHA, 'region_source_sha256': digest(REGION_SOURCE),
        'evaluator_sha256': digest(BIQUINTIC / 'freeeos_potential.py'), 'command': sys.argv,
        'nested_strides': list(STRIDES), 'stride1': 'Not a holdout: every printed native target would be retained; excluded from measured error ratios.',
        'coverage': coverage, 'eligible_odd_odd_native_material_own_queries': int(len(common)),
        'all_candidate_grid_common_queries': int(common.sum()), 'selected_queries': int(len(sample_q)),
        'selection': 'Both original indices odd; every target withheld at strides 2/4/8; one deterministic common sample across every candidate/grid.',
        'training': 'Constraints and missing twists recomputed only from retained coarse rows; validity metadata may inspect original flags/coordinates. No fine twists are reused.',
        'NAK4_control': 'Independent raw log/response columns on each four-node per-axis native-only not-a-knot window; exactly a local tensor cubic. Not a global spline or MESA runtime cache.',
        'donors': 'All derivative/window donors native/finite and clear of exclusions. Material/own-use are query and containing-cell restrictions, not a claim donors belong to an active MESA blend.',
        'region_bounds_from_015': region_bounds, 'region_counts': {name: int(mask.sum()) for name, mask in region_masks.items()},
        'region_note': 'Overlapping geometry bins; active_chiT is the 015 patch label, NOT a measured active fraction.',
        'results': results, 'paired_individual_error_ratios': ratio_results,
        'C3': 'Optional control omitted here to preserve a bounded primary-BQ/common-sample study; no cross-sample C3 convergence claim.',
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    paths = [OUTPUT / name for name in ['results.json', 'paired_queries.npz', 'error_summary.csv', 'SUMMARY.md', 'paired_refinement.png']]
    if not args.overwrite and any(path.exists() for path in paths):
        raise ValueError('Products exist; use --overwrite to repeat this experiment')
    paths[0].write_text(json.dumps(output, indent=2, allow_nan=False) + '\n')
    np.savez_compressed(paths[1], **products)
    with paths[2].open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['candidate', 'field', 'count', 'nonfinite', 'median_abs_relative_error', 'p99', 'max'])
        for label, candidate in results.items():
            for field, stats in candidate['relative_archived_field_errors'].items():
                writer.writerow([label, field, stats['count'], stats['nonfinite'], stats.get('median'), stats.get('p99'), stats.get('max')])
    lines = ['# Same-query native refinement control', '', output['claim'], '',
             f"{len(common)} eligible odd/odd nodes; {common.sum()} common to all grids/candidates; {len(sample_q)} selected.",
             'Measured strides 2/4/8 only; stride 1 retains these native targets and cannot be called heldout.', '',
             '| Candidate | Pgas p99 / max | E p99 / max | S p99 / max | Cv p99 / max | chiRho p99 / max | chiT p99 / max |',
             '| --- | --- | --- | --- | --- | --- | --- |']
    for label, candidate in results.items():
        cells = []
        for field in ['Pgas', 'E', 'S', 'Cv', 'chiRho', 'chiT']:
            stats = candidate['relative_archived_field_errors'][field]
            cells.append(f"{stats.get('p99', float('nan')):.3e} / {stats.get('max', float('nan')):.3e}")
        lines.append('| ' + label + ' | ' + ' | '.join(cells) + ' |')
    lines.extend(['', 'Region sample counts: ' + ', '.join(f'{name}={mask.sum()}' for name, mask in region_masks.items()), '',
                  'Individual coarse/fine error ratios, zero-denominator counts, and improvement fractions are in results.json/paired_queries.npz.',
                  'NAK4 is a local four-node not-a-knot column control, not an executed MESA cache or globally fitted fallback-contaminated spline.',
                  'Shrinking-mesh trends are empirical archived-field reproduction; no formal order is inferred.'])
    paths[3].write_text('\n'.join(lines) + '\n')
    (OUTPUT / 'README.md').write_text('\n'.join(lines) + '\n')
    plot_pairs(errors)
    print('\n'.join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-queries', type=int, default=4000)
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    if not 1 <= args.max_queries <= 10000:
        parser.error('Require 1..10000 common queries')
    run(args)


if __name__ == '__main__':
    main()
