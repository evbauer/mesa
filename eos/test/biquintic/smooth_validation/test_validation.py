"""Small analytic checks for the read-only smoother validation sidecar."""

import unittest

import numpy as np

from validate_multislice import fourth_partials
from freeeos_potential import EVAL_ORDERS, LOW_ORDERS
from smooth_patch import evaluate_fit, fit_constraints, fit_thermodynamic_constraints, spline_basis


def manufactured():
    grid_u, grid_tau = np.linspace(-.2, .4, 9), np.linspace(-.3, .5, 8)
    polynomial = np.zeros((6, 6))
    polynomial[0, 0], polynomial[1, 0], polynomial[0, 1] = 7., 1.2, -2.
    polynomial[2, 2], polynomial[4, 1] = .9, .3
    mesh_u, mesh_tau = np.meshgrid(grid_u, grid_tau, indexing="ij")
    low = np.full((3, 3, len(grid_u), len(grid_tau)), np.nan)
    for orders in LOW_ORDERS:
        derivative = np.polynomial.polynomial.polyder(polynomial, orders[0], axis=0)
        derivative = np.polynomial.polynomial.polyder(derivative, orders[1], axis=1)
        low[orders] = np.polynomial.polynomial.polyval2d(mesh_u, mesh_tau, derivative)
    return grid_u, grid_tau, polynomial, low


class SmoothValidationTests(unittest.TestCase):
    def test_basis_partition_through_fourth_derivative(self):
        grid, query = np.linspace(-.3, .5, 9), np.linspace(-.3, .5, 23)
        np.testing.assert_allclose(spline_basis(grid, query).sum(axis=1), 1, atol=2e-15)
        for order in range(1, 5):
            self.assertLess(np.max(abs(spline_basis(grid, query, order).sum(axis=1))), 2e-7)

    def check_polynomial(self, fitter, weight):
        grid_u, grid_tau, polynomial, low = manufactured()
        fitted, _ = fitter(low, grid_u, grid_tau, weight)
        generator = np.random.default_rng(17)
        sample_u, sample_tau = generator.uniform(grid_u[0], grid_u[-1], 20), generator.uniform(grid_tau[0], grid_tau[-1], 20)
        actual = evaluate_fit(fitted, grid_u, grid_tau, sample_u, sample_tau)
        for orders in EVAL_ORDERS:
            derivative = np.polynomial.polynomial.polyder(polynomial, orders[0], axis=0)
            derivative = np.polynomial.polynomial.polyder(derivative, orders[1], axis=1)
            expected = np.polynomial.polynomial.polyval2d(sample_u, sample_tau, derivative)
            np.testing.assert_allclose(actual[orders], expected, rtol=3e-6, atol=3e-6)
        return fitted, grid_u, grid_tau, polynomial

    def test_native_control_exact_polynomial(self):
        self.check_polynomial(fit_constraints, 1.)

    def test_thermodynamic_weights_exact_polynomial(self):
        for weight in (.01, .1, 1.):
            self.check_polynomial(fit_thermodynamic_constraints, weight)

    def test_fourth_partials_exact_polynomial(self):
        fitted, grid_u, grid_tau, polynomial = self.check_polynomial(fit_thermodynamic_constraints, 1.)
        sample_u, sample_tau = np.array([.1, .2]), np.array([.1, .2])
        actual = fourth_partials(fitted, grid_u, grid_tau, sample_u, sample_tau)
        for density_order in range(5):
            derivative = np.polynomial.polynomial.polyder(polynomial, density_order, axis=0)
            derivative = np.polynomial.polynomial.polyder(derivative, 4 - density_order, axis=1)
            expected = np.polynomial.polynomial.polyval2d(sample_u, sample_tau, derivative)
            np.testing.assert_allclose(actual[str((density_order, 4 - density_order))], expected, rtol=2e-5, atol=2e-5)

    def test_simple_knots_fourth_derivative_gap_converges(self):
        grid = np.linspace(-.4, .6, 9)
        coefficients = np.random.default_rng(18).normal(size=(13, 13))
        face, tangent = np.array([grid[4]]), np.array([.123])
        differences = []
        for step in (1e-4, 1e-5):
            left = fourth_partials(coefficients, grid, grid, face - step, tangent)
            right = fourth_partials(coefficients, grid, grid, face + step, tangent)
            differences.append(max(abs(left[field][0] - right[field][0]) for field in left))
        self.assertLess(differences[1], .12 * differences[0])


if __name__ == "__main__":
    unittest.main()
