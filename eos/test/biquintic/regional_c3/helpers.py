"""Exact archived-query and retained-source helpers for the regional C3 sidecar."""

from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
BIQUINTIC = HERE.parent
WORKSPACE = BIQUINTIC.parents[3]
sys.path.insert(0, str(BIQUINTIC))

from freeeos_potential import CRAD, LOW_ORDERS, source_constraints

FIELDS = ('P', 'Pgas', 'E', 'S', 'Cv', 'Cp', 'chiRho', 'chiT', 'gamma1', 'gamma3', 'grad_ad')
PATCHES = ('ionization', 'helium_ionization', 'dense_tail')
IMPORTANCES = (.01, .1)
LAMBDAS = (.001, .01, .1)


def parse_selection(selection, source_shape):
    if len(selection) != 4 or any(isinstance(value, bool) or not isinstance(value, (int, np.integer)) for value in selection):
        raise ValueError('015 selection must contain four integer indices/counts')
    start_q, start_t, count_q, count_t = map(int, selection)
    if min(start_q, start_t) < 0 or start_q % 2 or start_t % 2:
        raise ValueError('Require nonnegative globally even patch starts')
    if min(count_q, count_t) < 7 or count_q % 2 == 0 or count_t % 2 == 0:
        raise ValueError('Require odd patch dimensions with at least four retained nodes per axis')
    if start_q + count_q > source_shape[0] or start_t + count_t > source_shape[1]:
        raise ValueError('Patch outside source grid')
    return start_q, start_t, count_q, count_t


def paired_indices(selection, source_shape):
    start_q, start_t, count_q, count_t = parse_selection(selection, source_shape)
    retained_q = start_q + np.arange(0, count_q, 2)
    retained_t = start_t + np.arange(0, count_t, 2)
    withheld_q, withheld_t = np.meshgrid(start_q + np.arange(1, count_q, 2),
                                        start_t + np.arange(1, count_t, 2), indexing='ij')
    return retained_q, retained_t, withheld_q.ravel(), withheld_t.ravel()


def validate_query_identity(archive, logq, logt, selection):
    start_q, start_t, count_q, count_t = parse_selection(selection, (len(logq), len(logt)))
    retained_q, retained_t, withheld_q, withheld_t = paired_indices(selection, (len(logq), len(logt)))
    np.testing.assert_array_equal(archive['patch_logQ'], logq[start_q:start_q + count_q])
    np.testing.assert_array_equal(archive['patch_logT'], logt[start_t:start_t + count_t])
    np.testing.assert_array_equal(archive['withheld_logQ'], logq[withheld_q])
    np.testing.assert_array_equal(archive['withheld_logT'], logt[withheld_t])
    if archive['metrics'].tolist()[:11] != list(FIELDS):
        raise ValueError('Archived state field order differs from 015 contract')
    if archive['local_coarse_common'].shape != withheld_q.shape or archive['local_coarse_common'].dtype != bool:
        raise ValueError('Invalid archived local-donor subset flags')
    if np.intersect1d(retained_q, withheld_q).size or np.intersect1d(retained_t, withheld_t).size:
        raise ValueError('Retained/withheld source index leak')
    return retained_q, retained_t, withheld_q, withheld_t


def retained_constraints(rows, logq, logt, selection):
    retained_q, retained_t, withheld_q, withheld_t = paired_indices(selection, rows.shape[:2])
    retained_rows = rows[retained_q[:, None], retained_t[None, :]].copy()
    coarse_q, coarse_t = logq[retained_q], logt[retained_t]
    low, valid, _ = source_constraints(retained_rows, coarse_q, coarse_t, minimum_logT=3.2)
    if not valid.all() or any(not np.isfinite(low[orders]).all() for orders in LOW_ORDERS):
        raise ValueError('Retained native constraints must all be finite/valid')
    return low, coarse_q, coarse_t, retained_q, retained_t


def printed_reference(rows):
    rows = np.asarray(rows)
    if rows.shape[-1] != 22 or not np.isfinite(rows).all():
        raise ValueError('Expected finite 22-column printed source rows')
    reference = {'Pgas': 10.0**rows[..., 1], 'E': 10.0**rows[..., 2], 'S': 10.0**rows[..., 3],
                 'chiRho': rows[..., 4], 'chiT': rows[..., 5], 'Cp': rows[..., 6], 'Cv': rows[..., 7],
                 'gamma1': rows[..., 13], 'gamma3': rows[..., 14], 'grad_ad': rows[..., 15]}
    reference['P'] = reference['Pgas'] + CRAD * (10.0**rows[..., 0])**4 / 3
    if any(not np.isfinite(reference[field]).all() or np.any(reference[field] <= 0) for field in FIELDS):
        raise ValueError('Printed reference fields must be finite/positive for relative comparisons')
    return reference


def value_array(state):
    return np.stack([state[field] for field in FIELDS])


def summary(values):
    values = np.asarray(values)
    finite = np.isfinite(values)
    result = {'count': int(values.size), 'nonfinite': int((~finite).sum())}
    if finite.any():
        result.update(zip(('median', 'p99', 'max'), map(float, np.quantile(abs(values[finite]), [.5, .99, 1]))))
    return result


def errors(candidate, reference, sample_q, sample_t, original_q=None, original_t=None):
    result = {}
    for field in FIELDS:
        normalized = (candidate[field] - reference[field]) / abs(reference[field])
        index = int(np.argmax(abs(normalized)))
        worst = {'query_index_0based': index, 'logQ': float(sample_q[index]), 'logT': float(sample_t[index]),
                 'candidate': float(candidate[field][index]), 'reference': float(reference[field][index]),
                 'reference_scale': float(abs(reference[field][index])), 'signed_relative_difference': float(normalized[index])}
        if original_q is not None:
            worst.update(original_iq=int(original_q[index]), original_it=int(original_t[index]))
        result[field] = {**summary(normalized), 'worst': worst}
    return result


def positivity(state):
    return {field: {'count': int(np.asarray(state[field]).size), 'nonfinite': int((~np.isfinite(state[field])).sum()),
                    'nonpositive_finite': int((np.isfinite(state[field]) & (state[field] <= 0)).sum())} for field in FIELDS}


def paired_improvement(candidate, baseline, reference, sample_q, sample_t):
    result, individual = {}, {}
    for field in FIELDS:
        candidate_error = abs((candidate[field] - reference[field]) / reference[field])
        baseline_error = abs((baseline[field] - reference[field]) / reference[field])
        difference = baseline_error - candidate_error
        regression = int(np.argmin(difference))
        improvement = int(np.argmax(difference))
        result[field] = {'strictly_improved': int((difference > 0).sum()), 'strictly_worse': int((difference < 0).sum()),
                         'exact_ties': int((difference == 0).sum()), 'count': int(difference.size),
                         'improved_fraction': float(np.mean(difference > 0)), 'median_signed_error_reduction': float(np.median(difference)),
                         'largest_improvement': {'query_index_0based': improvement, 'logQ': float(sample_q[improvement]),
                             'logT': float(sample_t[improvement]), 'relative_error_reduction': float(difference[improvement])},
                         'largest_regression': {'query_index_0based': regression, 'logQ': float(sample_q[regression]),
                             'logT': float(sample_t[regression]), 'relative_error_reduction': float(difference[regression]),
                             'C3': float(candidate[field][regression]), 'C4': float(baseline[field][regression]),
                             'printed_reference': float(reference[field][regression])}}
        individual[field] = difference
    return result, np.stack([individual[field] for field in FIELDS])
