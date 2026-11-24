"""Synthetic analytic-radiation and polynomial-plus-radiation calibrations."""

import sys
import unittest

sys.dont_write_bytecode = True

import numpy as np

from bounds import DTYPE, POWER_TO_BERNSTEIN
from bound_split import adaptive_split, radiation_bounds, radiation_cv_numerator, radiation_guard
from freeeos_potential import CRAD, radiation_potential, thermodynamics


class SplitTests(unittest.TestCase):
    def test_analytic_radiation_derivatives_and_signs(self):
        sample_q, sample_t = np.array([-4., -1., 2.]), np.array([4., 6., 7.])
        native = radiation_potential(sample_q, sample_t)
        state = thermodynamics(native, sample_q, sample_t)
        expected = radiation_cv_numerator(sample_q, sample_t)
        np.testing.assert_allclose(state["Cv"]*state["T"], expected, rtol=3e-15)
        np.testing.assert_array_equal(native[1, 0]+native[2, 0], np.zeros(3))
        np.testing.assert_array_equal(state["P_rho"], np.zeros(3))
        self.assertTrue(np.all(expected > 0))

    def test_monotonic_radiation_rectangle_enclosure(self):
        logbox = tuple(map(DTYPE, (-3., -2.7, 4., 4.2)))
        lower, upper = radiation_bounds(logbox)
        self.assertEqual(lower, radiation_cv_numerator(logbox[1], logbox[2]))
        self.assertEqual(upper, radiation_cv_numerator(logbox[0], logbox[3]))
        for density_phase in np.linspace(0, 1, 13):
            for temperature_phase in np.linspace(0, 1, 13):
                value = radiation_cv_numerator(logbox[0]+(logbox[1]-logbox[0])*DTYPE(density_phase),
                                               logbox[2]+(logbox[3]-logbox[2])*DTYPE(temperature_phase))
                self.assertLessEqual(lower, value)
                self.assertGreaterEqual(upper, value)

    def test_negative_remainder_positive_after_radiation(self):
        logbox = tuple(map(DTYPE, (-3., -2.9, 4., 4.1)))
        lower, _ = radiation_bounds(logbox)
        coefficients = np.full((6, 6), -.5*lower, dtype=DTYPE)
        result = adaptive_split(coefficients, lower*DTYPE("1e-12"), logbox, "Cv_numerator")
        self.assertEqual(result["status"], "positive")
        self.assertGreater(result["lower_bound"], 0)
        self.assertEqual(result["visited_nodes"], 1)

    def test_combined_polynomial_and_radiation_enclosure(self):
        generator = np.random.default_rng(2801)
        coefficients = generator.normal(size=(6, 6)).astype(DTYPE)
        logbox = tuple(map(DTYPE, (0., .01, 0., .1)))
        lower, upper = radiation_bounds(logbox)
        from bounds import bernstein_value
        for density_phase, temperature_phase in generator.uniform(size=(31, 2)):
            value = bernstein_value(coefficients, DTYPE(density_phase), DTYPE(temperature_phase))
            value += radiation_cv_numerator(logbox[0]+(logbox[1]-logbox[0])*DTYPE(density_phase),
                                            logbox[2]+(logbox[3]-logbox[2])*DTYPE(temperature_phase))
            self.assertLessEqual(coefficients.min()+lower, value)
            self.assertGreaterEqual(coefficients.max()+upper, value)

    def test_correlated_negative_polynomial_needs_analytic_leaf_bounds(self):
        logbox = tuple(map(DTYPE, (0., .01, 0., .1)))
        radiation_min, _ = radiation_bounds(logbox)
        power = np.zeros((6, 6), dtype=DTYPE)
        power[0, 1] = -DTYPE("1.2")*radiation_min
        coefficients = POWER_TO_BERNSTEIN @ power @ POWER_TO_BERNSTEIN.T
        self.assertLess(coefficients.min()+radiation_min, 0)
        result = adaptive_split(coefficients, 1e-12, logbox, "Cv_numerator")
        self.assertEqual(result["status"], "positive")
        self.assertGreater(result["visited_nodes"], 1)
        self.assertEqual(result["unresolved_leaves"], 0)

    def test_prho_has_no_analytic_radiation_sign_term(self):
        logbox = tuple(map(DTYPE, (-3., -2.9, 4., 4.1)))
        result = adaptive_split(-np.ones((6, 6), dtype=DTYPE), 1e-12, logbox, "P_rho", restored_core=False, max_depth=0)
        self.assertEqual(result["status"], "unresolved")
        self.assertLess(result["lower_bound"], 0)

    def test_poly_plus_radiation_adaptive_cover(self):
        logbox = tuple(map(DTYPE, (0., 1., 0., 1.)))
        power = np.zeros((6, 6), dtype=DTYPE)
        power[0, 0], power[1, 0], power[2, 0] = DTYPE(".26"), -1, 1
        coefficients = POWER_TO_BERNSTEIN @ power @ POWER_TO_BERNSTEIN.T
        result = adaptive_split(coefficients, 1e-12, logbox, "Cv_numerator", restored_core=False)
        self.assertEqual(result["status"], "positive")
        self.assertGreater(result["visited_nodes"], 1)
        self.assertEqual(result["unresolved_leaves"], 0)

    def test_radiation_restoration_guard_is_explicit(self):
        _, upper = radiation_bounds(tuple(map(DTYPE, (-10., -9.9, 7., 7.1))))
        self.assertEqual(radiation_guard("P_rho", upper, False), 0)
        guard = radiation_guard("P_rho", upper, True)
        self.assertGreater(guard, 1)
        coefficients = np.full((6, 6), guard/2, dtype=DTYPE)
        logbox = tuple(map(DTYPE, (-10., -9.9, 7., 7.1)))
        representation = adaptive_split(coefficients, 1e-12, logbox, "P_rho", restored_core=False, max_depth=0)
        restored = adaptive_split(coefficients, 1e-12, logbox, "P_rho", restored_core=True, max_depth=0)
        self.assertEqual(representation["status"], "positive")
        self.assertEqual(restored["status"], "unresolved")

    def test_depth_and_nodes_remain_bounded(self):
        result = adaptive_split(np.zeros((6, 6), dtype=DTYPE), 1e-12,
                                tuple(map(DTYPE, (0., 1., 0., 1.))), "P_rho", max_depth=4, max_nodes=5)
        self.assertEqual(result["status"], "unresolved")
        self.assertEqual(result["visited_nodes"], 5)
        self.assertEqual(result["unresolved_leaves"], 4)


if __name__ == "__main__":
    unittest.main()
