"""Recorded-only final mask/resource/core-margin verification; no source or jet rebuild."""

from datetime import datetime
import json
import time

import numpy as np

from run_control import FIELDS, HERE, OUTPUT, WORKSPACE, checksum, domain_masks, new_cell_mask, support_guard


def coordinate_population_pair(first, second):
    populations = []
    for coordinates in (first, second):
        coordinates = np.asarray(coordinates)
        if coordinates.ndim != 2 or coordinates.shape[1] != 2 or coordinates.dtype.kind not in 'iu' or np.any(coordinates < 0):
            raise ValueError('Require nonnegative integer original cell coordinate pairs')
        population = {tuple(coordinate): index for index, coordinate in enumerate(coordinates)}
        if len(population) != len(coordinates):
            raise ValueError('Duplicate original cell coordinates')
        populations.append(population)
    first_map, second_map = populations
    common = sorted(first_map.keys() & second_map.keys())
    return (np.array([first_map[key] for key in common], dtype=int),
            np.array([second_map[key] for key in common], dtype=int),
            np.array([first_map[key] for key in sorted(first_map.keys() - second_map.keys())], dtype=int),
            np.array([second_map[key] for key in sorted(second_map.keys() - first_map.keys())], dtype=int))


def compare_archives(results):
    archives = {}
    for floor in ('3.2', '3.1'):
        with np.load(OUTPUT / f'floor_{floor}/cell_bounds_and_core_queries.npz', allow_pickle=False) as archive:
            archives[floor] = {key: archive[key] for key in archive.files}
    first, second = archives['3.2'], archives['3.1']
    coordinates = {floor: np.column_stack((archive['density_cell_index'], archive['temperature_cell_index']))
                   for floor, archive in archives.items()}
    common_first, common_second, only_first, only_second = coordinate_population_pair(coordinates['3.2'], coordinates['3.1'])
    np.testing.assert_array_equal(coordinates['3.2'][common_first], coordinates['3.1'][common_second])
    for key in ('logQ_axis', 'logT_axis'):
        np.testing.assert_array_equal(first[key], second[key])
    for key in ('sample_logQ', 'sample_logT'):
        np.testing.assert_array_equal(first[key][common_first], second[key][common_second])
    payload = {}
    summary = {'method': 'exact original (density_cell_index, temperature_cell_index) coordinate-set intersection/difference',
               'recorded_only': True, 'actual_MESA_fraction_measured_or_inferred': False,
               'overlap_jets_assumed_equal': False,
               'source_floor_caution': 'Changing source floor can change overlap donor support and jets; paired archived fields are compared, not assumed equal.',
               'status_codes': {'unsupported': -1, 'unresolved': 0, 'positive': 1, 'negative_witness': 2},
               'input_archive_sha256': {floor: result['archive_sha256'] for floor, result in results['candidates'].items()},
               'populations': {}, 'overlap_field_comparisons': {}}
    groups = {'overlap': {'3.2': common_first, '3.1': common_second},
              'only_floor_3.2': {'3.2': only_first}, 'only_floor_3.1': {'3.1': only_second}}
    for name, members in groups.items():
        reference_floor = next(iter(members))
        indices = members[reference_floor]
        cell_coordinates = coordinates[reference_floor][indices]
        archive = archives[reference_floor]
        payload[name + '_cell_coordinates'] = cell_coordinates
        population = {'cells': len(indices), 'supported_floors': list(members), 'fields_by_floor': {}}
        population['closed_cell_logT_min'] = float(archive['logT_axis'][cell_coordinates[:, 1]].min()) if len(indices) else None
        population['closed_cell_logT_max'] = float(archive['logT_axis'][cell_coordinates[:, 1] + 1].max()) if len(indices) else None
        for floor in ('3.2', '3.1'):
            prefix = name + '_floor_' + floor + '_'
            if floor not in members:
                population['fields_by_floor'][floor] = {'supported': False, 'unsupported_cells': len(indices)}
                for field in FIELDS:
                    payload[prefix + field + '_status'] = np.full(len(indices), -1, dtype=np.int8)
                continue
            selected = members[floor]
            current = archives[floor]
            count = len(current['density_cell_index'])
            payload[prefix + 'archive_index'] = selected
            for key, values in current.items():
                if (key.startswith(tuple(field + '_' for field in FIELDS)) or key in
                    ('sample_logQ', 'sample_logT', 'source_corner_valid', 'donor_corner_support',
                     'full_one_branch_own_use', 'material_corner_screen')):
                    if values.ndim and values.shape[0] == count:
                        payload[prefix + key] = values[selected]
            field_summary = {'supported': True, 'full_one_branch_own_use_cells': int(current['full_one_branch_own_use'][selected].sum())}
            for field in FIELDS:
                status = current[field + '_status'][selected]
                field_summary[field] = {label: int((status == code).sum()) for label, code in
                                        (('positive', 1), ('negative_witness', 2), ('unresolved', 0))}
            field_summary['both_positive'] = int(np.logical_and.reduce([current[field + '_status'][selected] == 1 for field in FIELDS]).sum())
            population['fields_by_floor'][floor] = field_summary
        summary['populations'][name] = population
    for field in FIELDS:
        difference = second[field + '_core_value'][common_second] - first[field + '_core_value'][common_first]
        scale = np.maximum(1., np.maximum(abs(first[field + '_core_value'][common_first]), abs(second[field + '_core_value'][common_second])))
        payload['overlap_' + field + '_core_difference_floor_3.1_minus_3.2'] = difference
        payload['overlap_' + field + '_core_difference_pair_scale'] = scale
        summary['overlap_field_comparisons'][field] = {
            'changed_status_cells': int((first[field + '_status'][common_first] != second[field + '_status'][common_second]).sum()),
            'changed_sample_values': int(np.count_nonzero(difference)),
            'max_abs_core_difference': float(abs(difference).max()),
            'max_pair_scaled_core_difference': float((abs(difference) / scale).max()),
            'pair_scale': 'max(1, abs(core_floor_3.2), abs(core_floor_3.1)); descriptive, not physical error'}
    np.savez_compressed(OUTPUT / 'floor_comparison.npz', **payload)
    summary['paired_archive_sha256'] = checksum(OUTPUT / 'floor_comparison.npz')
    (OUTPUT / 'floor_comparison.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    results['cross_floor_comparison'] = summary
    (OUTPUT / 'results.json').write_text(json.dumps(results, indent=2, allow_nan=False) + '\n')
    unique = summary['populations']['only_floor_3.1']
    measured = unique['fields_by_floor']['3.1']
    marker = '## Cross-floor coordinate-set classification'
    readme = (OUTPUT / 'README.md').read_text().split(marker)[0].rstrip()
    readme += f'\n\n{marker}\n\n'
    readme += f"Exact original-coordinate sets give {len(common_first):,} overlap cells, {len(only_first):,} floor3.2-only cells, and **{len(only_second):,} floor3.1-only cells**.\n\n"
    readme += '| Population | Closed-cell logT range | Cv positive / witness / unresolved | Prho positive / witness / unresolved | Both positive | 024 geometric own-use |\n|---|---|---|---|---:|---:|\n'
    counts = lambda field: ' / '.join(str(measured[field][key]) for key in ('positive', 'negative_witness', 'unresolved'))
    readme += f"| Floor3.1-only | {unique['closed_cell_logT_min']:.8g} to {unique['closed_cell_logT_max']:.8g} | {counts('Cv_numerator')} | {counts('P_rho')} | {measured['both_positive']} | {measured['full_one_branch_own_use_cells']} |\n\n"
    readme += 'The opposite-floor status is **unsupported (-1)**, not unresolved; no unavailable numeric values are invented. The unchanged024 predicate excludes this low-temperature population; no actual MESA fraction is measured or inferred. Lowering source floor can change overlap donor support and jets: equality is not assumed. `floor_comparison.npz` pairs all archived per-cell field values, bounds, guards, statuses, resources and four-query coordinates on overlap, and retains these for unique populations; `floor_comparison.json` records coordinate-set counts and observed overlap differences. This finalizer only reads completed archives; no new bound run, source/jet rebuild, or cache access.\n'
    (OUTPUT / 'README.md').write_text(readme)
    return summary


def main():
    started = time.monotonic()
    source_files = list(HERE.glob('*.py')) + [HERE / 'README.md']
    before = {str(path.relative_to(WORKSPACE)): checksum(path) for path in source_files}
    results = json.loads((OUTPUT / 'results.json').read_text())
    if results['status'] != 'complete' or not results['read_only_inputs_unchanged']:
        raise AssertionError('Require complete guarded audit and unchanged read-only inputs')
    checks = {}
    for floor, result in results['candidates'].items():
        path = OUTPUT / f'floor_{floor}/cell_bounds_and_core_queries.npz'
        if checksum(path) != result['archive_sha256'] or not result['calibration_passed']:
            raise AssertionError('Archive fingerprint/calibration mismatch')
        with np.load(path, allow_pickle=False) as archive:
            selected = new_cell_mask(archive['candidate_complete_node_mask'], archive['candidate_source_valid_node_mask'],
                                     archive['baseline_centered_complete_node_mask'], results['limits']['max_cells'])
            np.testing.assert_array_equal(selected, archive['selected_new_cell_mask'])
            density_indices, temperature_indices = np.nonzero(selected)
            np.testing.assert_array_equal(density_indices, archive['density_cell_index'])
            np.testing.assert_array_equal(temperature_indices, archive['temperature_cell_index'])
            own, _, _, _ = domain_masks(archive['logQ_axis'], archive['logT_axis'],
                                        np.ones((len(archive['logQ_axis']), len(archive['logT_axis']))))
            np.testing.assert_array_equal(own[selected], archive['full_one_branch_own_use'])
            if not archive['source_corner_valid'].all():
                raise AssertionError('Four-corner native guard failed')
            minimum = support_guard(archive['donor_corner_support'], archive['logT_axis'], float(floor))
            if minimum != result['minimum_selected_donor_logT']:
                raise AssertionError('Donor floor metadata differs')
            count = len(density_indices)
            for field in FIELDS:
                core = archive[field + '_core_value']
                prediction = archive[field + '_polynomial_plus_radiation_value']
                margin = archive[field + '_effective_margin']
                ratio = archive[field + '_difference_over_effective_margin']
                if core.shape != (count, 4) or not np.isfinite(core).all() or not np.isfinite(prediction).all() or not (margin > 0).all():
                    raise AssertionError('Missing/nonfinite four-query calibration or margin')
                np.testing.assert_allclose(abs(prediction - core) / margin[:, None], ratio, rtol=1e-13, atol=1e-15)
                if (ratio > 1).any() or not archive[field + '_processed'].all():
                    raise AssertionError('Out-of-margin/unprocessed guarded cell')
                status = archive[field + '_status']
                if np.any((status == 1)[:, None] & (core <= 0)):
                    raise AssertionError('Positive certificate contradicts restored-core sample')
                if archive[field + '_visited_nodes'].max() > results['limits']['max_nodes'] or archive[field + '_maximum_depth'].max() > results['limits']['max_depth']:
                    raise AssertionError('Adaptive resource cap exceeded')
                witnesses = result['negative_witnesses'][field]
                if len(witnesses) != int((status == 2).sum()):
                    raise AssertionError('Explicit witness/status counts differ')
                for evidence in witnesses:
                    if not np.isfinite(evidence['witness']['core_value']) or evidence['witness']['core_value'] >= 0:
                        raise AssertionError('Not a finite negative restored-core witness')
            checks[floor] = {'new_cells': count, 'four_queries_per_cell': True, 'new_native_masks_and_floor': True,
                             'geometric_own_use_identity': True, 'guard_margin_agreement': True, 'resource_caps': True,
                             'finite_negative_witnesses_only': True}
    comparison = compare_archives(results)
    after = {str(path.relative_to(WORKSPACE)): checksum(path) for path in source_files}
    if before != after:
        raise AssertionError('Code changed during final recorded-only verification')
    report = {'status': 'passed', 'timestamp_local': datetime.now().astimezone().isoformat(),
              'elapsed_seconds': time.monotonic() - started, 'source_or_jet_rebuild': False, 'MESA_or_cache_access': False,
              'code_unchanged_during_verification': True, 'frozen_source_hashes': before, 'candidates': checks,
              'cross_floor_comparison': comparison}
    (OUTPUT / 'validation_results.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(f'Recorded-only verification passed for both floors; {report["elapsed_seconds"]:.2f}s; code frozen.')


if __name__ == '__main__':
    main()
