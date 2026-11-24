"""Analytic/cardinal and invalid-donor checks for the optional edge control."""

import unittest

import numpy as np

from edge_control import build_jets, polynomial_derivative
from freeeos_potential import LOW_ORDERS, construct_twists, local_derivative


class EdgeControlTests(unittest.TestCase):
    def test_quartic_derivative_both_axes_including_endpoints(self):
        coordinate = np.arange(13) * .1
        valid = np.ones((13, 13), bool)
        for axis in (0, 1):
            values = np.broadcast_to((coordinate**4)[:, None] if axis == 0 else (coordinate**4)[None, :], valid.shape)
            expected = np.broadcast_to((4 * coordinate**3)[:, None] if axis == 0 else (4 * coordinate**3)[None, :], valid.shape)
            actual, accepted, metadata = polynomial_derivative(values, valid, .1, axis)
            self.assertTrue(accepted.all())
            np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
            line = metadata["position"][:, 0] if axis == 0 else metadata["position"][0]
            np.testing.assert_array_equal(line[:3], [0, 1, 2])
            np.testing.assert_array_equal(line[-3:], [2, 3, 4])

    def test_centered_only_matches_core_bitwise(self):
        values = np.random.default_rng(18).normal(size=(25, 25))
        valid = np.ones_like(values, bool)
        valid[12, 12] = False
        for axis in (0, 1):
            actual, accepted, _ = polynomial_derivative(values, valid, .1, axis, shifted=False)
            expected, expected_mask = local_derivative(values, valid, .1, axis, radius=2)
            np.testing.assert_array_equal(accepted, expected_mask)
            np.testing.assert_array_equal(actual[accepted], expected[accepted])

    def test_quartic_tensor_nine_jets_at_all_native_nodes(self):
        first, second = np.linspace(-.4, .6, 13), np.linspace(-.3, .7, 13)
        coefficients = np.random.default_rng(19).normal(size=(5, 5))
        mesh_first, mesh_second = np.meshgrid(first, second, indexing="ij")
        low = np.full((3, 3, 13, 13), np.nan)
        expected = np.empty_like(low)
        for density_order in range(3):
            for temperature_order in range(3):
                derivative = np.polynomial.polynomial.polyder(coefficients, density_order, axis=0)
                derivative = np.polynomial.polynomial.polyder(derivative, temperature_order, axis=1)
                expected[density_order, temperature_order] = np.polynomial.polynomial.polyval2d(mesh_first, mesh_second, derivative)
        for orders in LOW_ORDERS:
            low[orders] = expected[orders]
        nodes, complete, _ = build_jets(low, np.ones((13, 13), bool), first[1] - first[0], second[1] - second[0])
        self.assertTrue(complete.all())
        np.testing.assert_allclose(nodes, expected, rtol=2e-10, atol=2e-10)

    def test_invalid_poison_isolation_and_coverage_superset(self):
        low = np.random.default_rng(20).normal(size=(3, 3, 25, 25))
        valid = np.ones((25, 25), bool)
        valid[:, :5] = False
        valid[12, 12] = False
        changed = low.copy()
        changed[:, :, ~valid] = 1e200
        nodes, complete, metadata = build_jets(low, valid, .1, .1)
        altered, altered_complete, _ = build_jets(changed, valid, .1, .1)
        core_nodes, core_complete, _ = construct_twists(low, valid, .1, .1)
        np.testing.assert_array_equal(complete, altered_complete)
        np.testing.assert_array_equal(nodes[:, :, complete], altered[:, :, complete])
        np.testing.assert_array_equal(nodes[:, :, core_complete], core_nodes[:, :, core_complete])
        self.assertFalse(np.any(core_complete & ~complete))
        self.assertFalse(complete[:, :5].any())
        self.assertFalse(complete[12, 12])
        self.assertGreater(complete.sum(), core_complete.sum())
        self.assertGreaterEqual(metadata["support"][..., 2][complete].min(), 5)

    def test_short_run_cannot_bridge_invalid_gap(self):
        values = np.ones((11, 11))
        valid = np.ones_like(values, bool)
        valid[:, 4:7] = False
        derivative, accepted, _ = polynomial_derivative(values, valid, .1, 1)
        self.assertFalse(accepted.any())
        self.assertTrue(np.isnan(derivative).all())


if __name__ == "__main__":
    unittest.main()
