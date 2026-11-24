"""Analytic thermodynamic invariance under constant and linear-T gauges."""

import unittest

import numpy as np

from freeeos_potential import EVAL_ORDERS, to_physical
from gauge_control import apply_gauge, gauge_partials


class GaugeControlTests(unittest.TestCase):
    def test_analytic_gauge_physical_derivatives_and_invariants(self):
        sample_q = np.array([-1., 0., 1.])
        sample_t = np.array([3.5, 4., 5.])
        constant, coefficient = 1e12, 1e8
        native = gauge_partials(sample_q, sample_t, constant, coefficient)
        physical = to_physical(native, EVAL_ORDERS)
        for orders in EVAL_ORDERS:
            if orders[0]:
                np.testing.assert_array_equal(physical[orders], 0.)
        np.testing.assert_allclose(physical[0, 0] - physical[0, 1], constant)
        np.testing.assert_allclose(-physical[0, 1] / 10**sample_t, -coefficient)
        np.testing.assert_array_equal(physical[0, 1] - physical[0, 2], 0.)
        np.testing.assert_array_equal(2 * physical[0, 2] - physical[0, 1] - physical[0, 3], 0.)

    def test_grid_gauge_addition_leaves_density_and_mixed_jets_unchanged(self):
        low = np.zeros((3, 3, 4, 5))
        shifted = apply_gauge(low, np.zeros((4, 1)), np.ones((1, 5)), 3., 7.)
        np.testing.assert_array_equal(low, 0.)
        np.testing.assert_array_equal(shifted[1:, :], 0.)
        np.testing.assert_array_equal(shifted[0, 0], 73.)
        np.testing.assert_array_equal(shifted[0, 1], 70.)
        np.testing.assert_array_equal(shifted[0, 2], 70.)


if __name__ == "__main__":
    unittest.main()
