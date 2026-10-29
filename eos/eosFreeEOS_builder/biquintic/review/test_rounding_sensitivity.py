"""Focused tests for printed-decimal intervals and perturbation isolation."""

from decimal import Decimal
import unittest

import numpy as np

from rounding_sensitivity import FIELDS, compare_face, perturb, printed_interval, stencil_support, variability


class RoundingSensitivityTests(unittest.TestCase):
    def test_printed_exponent_and_inward_bounds(self):
        for token, expected_half in (("1.23456789012345E+15", 5.),
                                     ("-1.23456789012345E-02", 5e-17),
                                     ("0.00000000000000E+00", 5e-15),
                                     ("2.00", .005)):
            half, lower, upper = printed_interval(token)
            self.assertEqual(half, expected_half)
            center = Decimal(token)
            decimal_half = Decimal(5).scaleb(center.as_tuple().exponent - 1)
            self.assertGreaterEqual(Decimal.from_float(float(lower)), center - decimal_half)
            self.assertLessEqual(Decimal.from_float(float(upper)), center + decimal_half)

    def test_deterministic_shared_draws_and_preservation(self):
        rows = np.ones((2, 3, 22))
        indices = np.array([1, 4])
        half = np.full((2, len(FIELDS)), 5e-15)
        lower, upper = 1 - half, 1 + half
        narrow = perturb(rows, indices, lower, upper, (2, 3), 18)
        broad = perturb(rows, indices, lower, upper, tuple(FIELDS), 18)
        np.testing.assert_array_equal(narrow[:, :, [2, 3]], broad[:, :, [2, 3]])
        np.testing.assert_array_equal(narrow[:, :, [0, 17, 18, 19, 20, 21]], rows[:, :, [0, 17, 18, 19, 20, 21]])
        untouched = np.array([0, 2, 3, 5])
        np.testing.assert_array_equal(broad.reshape(-1, 22)[untouched], rows.reshape(-1, 22)[untouched])
        np.testing.assert_array_equal(narrow, perturb(rows, indices, lower, upper, (2, 3), 18))

    def test_face_comparison_does_not_hide_zero_or_small_jumps(self):
        baseline = {"left": {"Cv": np.array([1e8]), "Cv_t": np.array([2.])},
                    "right": {"Cv": np.array([1e8]), "Cv_t": np.array([1.])}}
        changed = {"left": {"Cv": np.array([1e8]), "Cv_t": np.array([2.1])},
                   "right": {"Cv": np.array([1e8]), "Cv_t": np.array([1.])}}
        result = compare_face([changed], baseline, "Cv_t", [18])
        self.assertAlmostEqual(result["max_jump_noise_over_abs_baseline_jump"], .1)
        self.assertFalse(result["sampled_rounding_negligible_vs_jump"])
        self.assertEqual(result["worst_jump_change"]["seed"], 18)
        baseline["right"]["Cv_t"] = np.array([2.])
        result = compare_face([changed], baseline, "Cv_t", [18])
        self.assertIsNone(result["max_jump_noise_over_abs_baseline_jump"])
        self.assertFalse(result["sampled_rounding_negligible_vs_jump"])

    def test_support_includes_two_stage_radius_and_excludes_invalid(self):
        selected = np.zeros((13, 13), dtype=bool)
        selected[6, 6] = True
        valid = np.ones_like(selected)
        valid[6, 7] = False
        support = stencil_support(selected, valid)
        self.assertTrue(support[2, 2])
        self.assertTrue(support[10, 10])
        self.assertFalse(support[1, 6])
        self.assertFalse(support[6, 7])
        self.assertEqual(support.sum(), 80)

    def test_distinct_absolute_and_scaled_extrema(self):
        result = variability(np.array([[102., 1.5]]), np.array([100., 1.]), np.array([100., 1.]),
                             [18], np.array([2., 3.]), np.array([6., 7.]))
        self.assertEqual(result["worst_absolute_change"]["logQ"], 2.)
        self.assertEqual(result["worst_field_scaled_change"]["logQ"], 3.)
        self.assertEqual(result["field_scaled_change"]["max"], .5)


if __name__ == "__main__":
    unittest.main()
