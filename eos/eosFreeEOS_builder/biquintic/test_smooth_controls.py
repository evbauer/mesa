"""Analytic checks for restricted smooth-potential feasibility controls."""

import unittest

import numpy as np

from freeeos_potential import EVAL_ORDERS, LOW_ORDERS
from smooth_curvature import normal_third_jump, reconcile_curvatures
from smooth_patch import evaluate_fit, fit_constraints, fit_thermodynamic_constraints, spline_basis


class SmoothControlTests(unittest.TestCase):
    def test_c3_reconciliation_quintic(self):
        coordinate = np.linspace(-.8, .9, 23)
        polynomial = np.array([.2, .9, -.7, .3, .6, -.1])
        values = np.polynomial.polynomial.polyval(coordinate, polynomial)
        slopes = np.polynomial.polynomial.polyval(coordinate, np.polynomial.polynomial.polyder(polynomial))
        curvature = np.polynomial.polynomial.polyval(coordinate, np.polynomial.polynomial.polyder(polynomial, 2))
        restored = reconcile_curvatures(values, slopes, curvature, coordinate[1] - coordinate[0])
        np.testing.assert_allclose(restored, curvature, rtol=5e-12, atol=5e-12)
        self.assertLess(np.max(abs(normal_third_jump(values, slopes, restored, coordinate[1] - coordinate[0]))), 5e-11)

    def test_basis_partition_and_polynomial_derivatives(self):
        grid = np.linspace(-.3, .5, 9)
        query = np.linspace(-.3, .5, 23)
        np.testing.assert_allclose(spline_basis(grid, query).sum(axis=1), 1, atol=2e-15)
        for order in (1, 2, 3):
            self.assertLess(np.max(abs(spline_basis(grid, query, order).sum(axis=1))), 2e-9)

    def test_soft_constraint_fit_tensor_polynomial(self):
        first = np.linspace(-.2, .4, 9)
        second = np.linspace(-.3, .5, 8)
        coefficients = np.zeros((6, 6))
        coefficients[0, 0], coefficients[1, 0], coefficients[0, 1] = 7., 1.2, -2.
        coefficients[2, 2], coefficients[4, 1] = .9, .3
        low = np.full((3, 3, len(first), len(second)), np.nan)
        mesh_first, mesh_second = np.meshgrid(first, second, indexing="ij")
        for orders in LOW_ORDERS:
            polynomial = np.polynomial.polynomial.polyder(coefficients, orders[0], axis=0)
            polynomial = np.polynomial.polynomial.polyder(polynomial, orders[1], axis=1)
            low[orders] = np.polynomial.polynomial.polyval2d(mesh_first, mesh_second, polynomial)
        generator = np.random.default_rng(17)
        query_first = generator.uniform(first[0], first[-1], 20)
        query_second = generator.uniform(second[0], second[-1], 20)
        for fitter in (fit_constraints, fit_thermodynamic_constraints):
            fitted, _ = fitter(low, first, second, 1.)
            actual = evaluate_fit(fitted, first, second, query_first, query_second)
            for orders in EVAL_ORDERS:
                polynomial = np.polynomial.polynomial.polyder(coefficients, orders[0], axis=0)
                polynomial = np.polynomial.polynomial.polyder(polynomial, orders[1], axis=1)
                expected = np.polynomial.polynomial.polyval2d(query_first, query_second, polynomial)
                np.testing.assert_allclose(actual[orders], expected, rtol=3e-6, atol=3e-6)


if __name__ == "__main__":
    unittest.main()
