"""Parsing, exact saved-query identity, reference columns, and retained-only leak checks."""

import unittest

import numpy as np

from helpers import (FIELDS, paired_improvement, parse_selection, printed_reference,
                     retained_constraints, validate_query_identity)


def synthetic_rows(logq, logt):
    rows = np.zeros((len(logq), len(logt), 22))
    rows[..., 0] = logt[None, :]
    rows[..., 1:4] = [5., 12., 8.]
    rows[..., 4:8] = [1., 1., 2e7, 1e7]
    rows[..., 13:16] = [1.5, 1.4, .2]
    rows[..., 18] = logq[:, None] + 2 * logt[None, :] - 12
    return rows


class RegionalTests(unittest.TestCase):
    def test_objective_normalization_separates_importance_squared_and_lambda(self):
        from run_control import objective_normalization
        solve = {'source_rows': 120, 'potential_scale': 1e12, 'relative_source_floors': {'E': 2.},
                 'data_objective': 3e-8, 'roughness_objective': 4e-5, 'roughness': {'formula': 'declared API formula'}}
        low = objective_normalization(solve, .01, .001)
        high = objective_normalization(solve, .1, .001)
        self.assertAlmostEqual(high['source_response_importance_squared'] / low['source_response_importance_squared'], 100.)
        self.assertAlmostEqual(low['source_operator_and_target_multiplier'], 1 / np.sqrt(120))
        self.assertAlmostEqual(low['total_augmented_objective'], 7e-8)
        self.assertEqual(low['source_mean_divisor'], high['source_mean_divisor'])
        self.assertEqual(low['roughness_operator_multiplier'], high['roughness_operator_multiplier'])

    def test_selection_parsing_requires_original_even_starts_and_odd_dimensions(self):
        self.assertEqual(parse_selection([2, 4, 9, 7], (15, 15)), (2, 4, 9, 7))
        for selection in [[True, 0, 9, 9], [1, 0, 9, 9], [0, 0, 8, 9], [0, 0, 5, 9], [8, 0, 9, 9]]:
            with self.assertRaises(ValueError):
                parse_selection(selection, (15, 15))

    def test_exact_query_identity_rejects_reordered_or_displaced_holdouts(self):
        logq, logt = np.arange(11) * .03, 3.5 + np.arange(11) * .02
        withheld_q, withheld_t = np.meshgrid(np.arange(1, 9, 2), np.arange(1, 9, 2), indexing='ij')
        archive = {'patch_logQ': logq[:9], 'patch_logT': logt[:9], 'withheld_logQ': logq[withheld_q.ravel()],
                   'withheld_logT': logt[withheld_t.ravel()], 'metrics': np.array(FIELDS),
                   'local_coarse_common': np.zeros(16, dtype=bool)}
        retained_q, retained_t, actual_q, actual_t = validate_query_identity(archive, logq, logt, [0, 0, 9, 9])
        np.testing.assert_array_equal(actual_q, withheld_q.ravel())
        np.testing.assert_array_equal(actual_t, withheld_t.ravel())
        np.testing.assert_array_equal(retained_q, np.arange(0, 9, 2))
        displaced = dict(archive, withheld_logT=archive['withheld_logT'][::-1])
        with self.assertRaises(AssertionError):
            validate_query_identity(displaced, logq, logt, [0, 0, 9, 9])

    def test_retained_only_zero_leak_from_withheld_cross_rows_and_outside_patch(self):
        logq, logt = -3 + np.arange(11) * .03, 3.5 + np.arange(11) * .02
        rows = synthetic_rows(logq, logt)
        original = retained_constraints(rows, logq, logt, [0, 0, 9, 9])
        keep = np.zeros(rows.shape[:2], dtype=bool)
        keep[np.arange(0, 9, 2)[:, None], np.arange(0, 9, 2)[None, :]] = True
        perturbed = rows.copy()
        perturbed[~keep] = np.random.default_rng(32).normal(size=((~keep).sum(), 22)) * 1e20
        changed = retained_constraints(perturbed, logq, logt, [0, 0, 9, 9])
        for actual, expected in zip(changed, original):
            np.testing.assert_array_equal(actual, expected)

    def test_printed_reference_columns_use_13_14_15_and_total_pressure(self):
        rows = synthetic_rows(np.array([-3.]), np.array([3.5]))[0]
        rows[:, 8:11] = [100, 200, 300]
        reference = printed_reference(rows)
        self.assertEqual(reference['gamma1'][0], 1.5)
        self.assertEqual(reference['gamma3'][0], 1.4)
        self.assertEqual(reference['grad_ad'][0], .2)
        self.assertGreater(reference['P'][0], reference['Pgas'][0])

    def test_individual_improvement_is_paired_to_same_reference_and_order(self):
        reference = {field: np.array([10., 10.]) for field in FIELDS}
        baseline = {field: np.array([11., 12.]) for field in FIELDS}
        candidate = {field: np.array([10.5, 13.]) for field in FIELDS}
        result, individual = paired_improvement(candidate, baseline, reference, [0., 1.], [4., 4.])
        self.assertEqual(result['Cv']['strictly_improved'], 1)
        self.assertEqual(result['Cv']['strictly_worse'], 1)
        np.testing.assert_allclose(individual, np.tile([.05, -.1], (11, 1)))


if __name__ == '__main__':
    unittest.main()
