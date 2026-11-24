"""Exact-quadrature, null-mode selection and regularized analytic fit checks."""

import unittest

import numpy as np

from regularized import (
    evaluate_fit, insert_double_knots, null_witness, polynomial_coefficients, roughness_operator, solve_regularized,
)
from spline_control import LOW_ORDERS


class RegularizationTests(unittest.TestCase):
    def test_knot_insertion_preserves_original_C4_potential(self):
        grid_u, grid_tau = np.linspace(-.3, .5, 9), np.linspace(-.4, .6, 8)
        coefficients = np.random.default_rng(31).normal(size=(len(grid_u) + 4, len(grid_tau) + 4))
        inserted = insert_double_knots(coefficients, grid_u, grid_tau)
        sample_u, sample_tau = np.array([-.13, .23]), np.array([.13, .23])
        old = evaluate_fit(coefficients, grid_u, grid_tau, sample_u, sample_tau, multiplicity=1)
        new = evaluate_fit(inserted, grid_u, grid_tau, sample_u, sample_tau)
        for orders in [(0, 0), (3, 0), (2, 1), (1, 2), (0, 3)]:
            np.testing.assert_allclose(old[orders], new[orders], atol=1e-10, rtol=1e-11)

    def test_exact_normal_fourth_polynomial_quadrature(self):
        grid_u, grid_tau = np.linspace(-.3, .5, 9), np.linspace(-.4, .6, 8)
        polynomial = np.zeros((6, 6))
        polynomial[5, 2] = .7
        fitted = polynomial_coefficients(polynomial, grid_u, grid_tau)
        rough, meta = roughness_operator(grid_u, grid_tau)
        actual = np.sum((rough @ fitted.ravel(order='F'))**2)
        spacing = np.diff(grid_u)[0]
        area = (grid_u[-1] - grid_u[0]) * (grid_tau[-1] - grid_tau[0])
        expected = spacing**8 / area * (.7 * 120)**2 * (grid_u[-1]**3 - grid_u[0]**3) / 3 * (grid_tau[-1]**5 - grid_tau[0]**5) / 5
        np.testing.assert_allclose(actual, expected, rtol=1e-9, atol=1e-13)

    def test_roughness_eliminates_nonzero_source_null_witness(self):
        grid_u, grid_tau = np.linspace(-.3, .5, 9), np.linspace(-.4, .6, 8)
        witness, _ = null_witness(grid_u, grid_tau)
        rough, _ = roughness_operator(grid_u, grid_tau)
        self.assertGreater(np.sum((rough @ witness.ravel(order='F'))**2), 1e-4)

    def test_bicubic_exact_zero_roughness_full_rank_and_stationarity(self):
        grid_u, grid_tau = np.linspace(-.3, .5, 9), np.linspace(-.4, .6, 8)
        polynomial = np.zeros((4, 4))
        polynomial[0, 0], polynomial[1, 0], polynomial[0, 1], polynomial[2, 2], polynomial[3, 1] = 7., 1.2, -2., .9, .3
        mesh_u, mesh_tau = np.meshgrid(grid_u, grid_tau, indexing='ij')
        low = np.full((3, 3, len(grid_u), len(grid_tau)), np.nan)
        for first, second in LOW_ORDERS:
            derivative = np.polynomial.polynomial.polyder(polynomial, first, axis=0)
            derivative = np.polynomial.polynomial.polyder(derivative, second, axis=1)
            low[first, second] = np.polynomial.polynomial.polyval2d(mesh_u, mesh_tau, derivative)
        fitted, diagnostics = solve_regularized(low, grid_u, grid_tau, .01)
        self.assertEqual(diagnostics['rank'], diagnostics['unknowns'])
        self.assertLess(diagnostics['roughness_objective'], 1e-16)
        self.assertLess(diagnostics['augmented_stationarity_relative'], 1e-10)
        sample_u, sample_tau = np.array([-.13, .23]), np.array([.13, .23])
        actual = evaluate_fit(fitted, grid_u, grid_tau, sample_u, sample_tau)
        for first, second in [(0, 0), (3, 0), (2, 1), (1, 2), (0, 3)]:
            derivative = np.polynomial.polynomial.polyder(polynomial, first, axis=0)
            derivative = np.polynomial.polynomial.polyder(derivative, second, axis=1)
            expected = np.polynomial.polynomial.polyval2d(sample_u, sample_tau, derivative)
            np.testing.assert_allclose(actual[first, second], expected, rtol=2e-5, atol=2e-5)


if __name__ == '__main__':
    unittest.main()
