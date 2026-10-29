"""Nested retained-grid controls on one common, never-retained source holdout set."""

import sys
from pathlib import Path

import numpy as np

SIDECAR = Path(__file__).resolve().parent
BIQUINTIC = SIDECAR.parent
sys.path.insert(0, str(BIQUINTIC))

from freeeos_potential import (
    CRAD, LN10, LOW_ORDERS, construct_twists, evaluate, query_cells,
    source_constraints, thermodynamics, valid_cells,
)


STRIDES = (2, 4, 8)
FIELDS = ('Pgas', 'E', 'S', 'Cv', 'Cp', 'chiRho', 'chiT')


def own_use(logq, logt):
    return (logt >= 3.2) & ((logt >= 6.5) | (logq <= -.5) | (logq + 2 * logt - 12 <= -6.5))


def native_mask(rows, logq, logt):
    valid = (rows[..., 17] == 0) & (logt[None, :] >= 3.2)
    valid &= np.all(np.isfinite(rows), axis=-1)
    expected_rho = logq[:, None] + 2 * logt[None, :] - 12
    valid &= abs(rows[..., 18] - expected_rho) < 1e-10
    valid &= np.all(rows[..., [4, 5, 6, 7]] > 0, axis=-1)
    if rows.shape[0] > 520:
        valid[520, :] = False
    return valid


def material_mask(rows, logq, logt):
    ratio = CRAD * (10.0**logt[None, :])**4 / (3 * 10.0**rows[..., 1])
    return (ratio < 1e-3) & (logq[:, None] >= -8) & (logq[:, None] <= 5) & (logt[None, :] <= 7.8)


def window_cells(valid):
    cells = np.zeros((valid.shape[0] - 1, valid.shape[1] - 1), dtype=bool)
    interior = np.ones((valid.shape[0] - 3, valid.shape[1] - 3), dtype=bool)
    for density_offset in range(4):
        for temperature_offset in range(4):
            interior &= valid[density_offset:density_offset + interior.shape[0],
                              temperature_offset:temperature_offset + interior.shape[1]]
    cells[1:-1, 1:-1] = interior
    return cells


def lagrange_four(position):
    locations = np.array([-1., 0., 1., 2.])
    weights = np.ones((len(position), 4))
    for node_index, node in enumerate(locations):
        for other_index, other in enumerate(locations):
            if node_index != other_index:
                weights[:, node_index] *= (position - other) / (node - other)
    return weights


def not_a_knot_columns(rows, logq, logt, sample_q, sample_t):
    density_index, temperature_index, inside = query_cells(logq, logt, sample_q, sample_t)
    if not np.all(inside & window_cells(np.ones(rows.shape[:2], dtype=bool))[density_index, temperature_index]):
        raise ValueError('Four-node not-a-knot window outside coarse grid')
    density_weights = lagrange_four((sample_q - logq[density_index]) / (logq[1] - logq[0]))
    temperature_weights = lagrange_four((sample_t - logt[temperature_index]) / (logt[1] - logt[0]))
    columns = np.zeros((len(sample_q), 22))
    for density_offset in range(4):
        for temperature_offset in range(4):
            columns += (density_weights[:, density_offset] * temperature_weights[:, temperature_offset])[:, None] * rows[
                density_index + density_offset - 1, temperature_index + temperature_offset - 1]
    with np.errstate(over='ignore', invalid='ignore'):
        return {'Pgas': 10.0**columns[:, 1], 'E': 10.0**columns[:, 2], 'S': 10.0**columns[:, 3],
                'Cv': columns[:, 7], 'Cp': columns[:, 6], 'chiRho': columns[:, 4], 'chiT': columns[:, 5]}


def build_grids(rows, logq, logt):
    grids = {}
    full_valid = native_mask(rows, logq, logt)
    for stride in STRIDES:
        coarse_q, coarse_t = logq[::stride], logt[::stride]
        coarse_rows = rows[::stride, ::stride].copy()
        valid = full_valid[::stride, ::stride]
        candidates = {}
        for split in (False, True):
            low, source_valid, _ = source_constraints(coarse_rows, coarse_q, coarse_t, split_radiation=split)
            source_valid &= valid
            nodes, complete, _ = construct_twists(low, source_valid,
                (coarse_q[1] - coarse_q[0]) * LN10, (coarse_t[1] - coarse_t[0]) * LN10,
                radius=2, route='mean')
            candidates['BQ_split' if split else 'BQ_total'] = {
                'nodes': nodes, 'complete': complete, 'cells': valid_cells(complete), 'split': split}
        candidates['NAK4_columns'] = {'cells': window_cells(valid)}
        grids[stride] = {'rows': coarse_rows, 'q': coarse_q, 't': coarse_t, 'candidates': candidates,
                         'material_cells': valid_cells(material_mask(coarse_rows, coarse_q, coarse_t)
                             & own_use(coarse_q[:, None], coarse_t[None, :]))}
    return grids


def common_queries(rows, logq, logt, grids):
    density_indices, temperature_indices = np.indices(rows.shape[:2])
    eligible = (density_indices % 2 == 1) & (temperature_indices % 2 == 1)
    eligible &= native_mask(rows, logq, logt) & material_mask(rows, logq, logt)
    eligible &= own_use(logq[:, None], logt[None, :])
    density_indices, temperature_indices = np.nonzero(eligible)
    sample_q, sample_t = logq[density_indices], logt[temperature_indices]
    common = np.ones(len(sample_q), dtype=bool)
    coverage = {}
    for stride, grid in grids.items():
        density_index, temperature_index, inside = query_cells(grid['q'], grid['t'], sample_q, sample_t)
        for name, candidate in grid['candidates'].items():
            accepted = inside & grid['material_cells'][density_index, temperature_index]
            accepted &= candidate['cells'][density_index, temperature_index]
            coverage[f'{name}_stride{stride}'] = {'eligible_queries': int(len(common)),
                                               'valid_queries': int(accepted.sum())}
            common &= accepted
    return density_indices, temperature_indices, common, coverage


def evaluate_candidate(name, grid, sample_q, sample_t):
    if name == 'NAK4_columns':
        return not_a_knot_columns(grid['rows'], grid['q'], grid['t'], sample_q, sample_t)
    candidate = grid['candidates'][name]
    native = evaluate(candidate['nodes'], candidate['complete'], grid['q'], grid['t'], sample_q, sample_t,
                      orders=LOW_ORDERS, add_radiation=candidate['split'])
    return thermodynamics(native, sample_q, sample_t)


def table_reference(rows):
    return {'Pgas': 10.0**rows[..., 1], 'E': 10.0**rows[..., 2], 'S': 10.0**rows[..., 3],
            'Cv': rows[..., 7], 'Cp': rows[..., 6], 'chiRho': rows[..., 4], 'chiT': rows[..., 5]}


def error_summary(error):
    finite = np.asarray(error)[np.isfinite(error)]
    result = {'count': int(np.size(error)), 'finite': int(finite.size), 'nonfinite': int(np.size(error) - finite.size)}
    if finite.size:
        result.update(zip(['median', 'p99', 'max'], map(float, np.quantile(abs(finite), [.5, .99, 1]))))
    return result


def paired_ratios(coarse_error, fine_error):
    coarse, fine = abs(coarse_error), abs(fine_error)
    valid = np.isfinite(coarse) & np.isfinite(fine)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratios = coarse / fine
    finite_ratio = valid & np.isfinite(ratios)
    return ratios, {
        'finite_pair_count': int(valid.sum()), 'finite_ratio_count': int(finite_ratio.sum()),
        'fine_zero_coarse_positive': int((valid & (fine == 0) & (coarse > 0)).sum()),
        'both_zero_undefined': int((valid & (fine == 0) & (coarse == 0)).sum()),
        'fraction_improved_with_refinement': float((valid & (fine < coarse)).sum() / max(1, valid.sum())),
        'fraction_worsened_with_refinement': float((valid & (fine > coarse)).sum() / max(1, valid.sum())),
        'coarse_over_fine_abs_error_ratio': error_summary(ratios[finite_ratio]),
    }
