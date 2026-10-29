"""Analytic controls for double-knot C3 tensor quintic basis and soft fit."""

import unittest

import numpy as np

from spline_control import (
    EVAL_ORDERS, LOW_ORDERS, RankDeficientFit, evaluate_fit, fit_thermodynamic_constraints,
    knots_for, null_witness, polynomial_coefficients, spline_basis, thermodynamic_combinations,
)
from freeeos_potential import thermodynamic_partials


def polynomial_data():
    grid_u, grid_tau = np.linspace(-.3, .5, 9), np.linspace(-.4, .6, 8)
    coefficients = np.zeros((6, 6))
    coefficients[0, 0], coefficients[1, 0], coefficients[0, 1] = 7., 1.2, -2.
    coefficients[2, 2], coefficients[4, 1], coefficients[5, 3] = .9, .3, .02
    mesh_u, mesh_tau = np.meshgrid(grid_u, grid_tau, indexing='ij')
    low = np.full((3, 3, len(grid_u), len(grid_tau)), np.nan)
    for first, second in LOW_ORDERS:
        derivative = np.polynomial.polynomial.polyder(coefficients, first, axis=0)
        derivative = np.polynomial.polynomial.polyder(derivative, second, axis=1)
        low[first, second] = np.polynomial.polynomial.polyval2d(mesh_u, mesh_tau, derivative)
    return grid_u, grid_tau, coefficients, low


class SplineTests(unittest.TestCase):
    def test_double_knots_partition_and_derivative_partition(self):
        grid = np.linspace(-.4, .6, 11)
        knots = knots_for(grid)
        self.assertEqual(len(knots) - 6, 2 * len(grid) + 2)
        self.assertTrue(all(np.count_nonzero(knots == knot) == 2 for knot in grid[1:-1]))
        query = np.linspace(grid[0], grid[-1], 41)
        np.testing.assert_allclose(spline_basis(grid, query).sum(axis=1), 1, atol=2e-15)
        for order in range(1, 6):
            np.testing.assert_allclose(spline_basis(grid, query, order).sum(axis=1) * .1**order, 0, atol=2e-12)

    def test_polynomial_representation_not_unique_six_constraint_fit(self):
        grid_u, grid_tau, polynomial, low = polynomial_data()
        sample_u = np.random.default_rng(29).uniform(grid_u[0], grid_u[-1], 24)
        sample_tau = np.random.default_rng(30).uniform(grid_tau[0], grid_tau[-1], 24)
        fitted = polynomial_coefficients(polynomial, grid_u, grid_tau)
        actual = evaluate_fit(fitted, grid_u, grid_tau, sample_u, sample_tau)
        for first, second in EVAL_ORDERS:
            derivative = np.polynomial.polynomial.polyder(polynomial, first, axis=0)
            derivative = np.polynomial.polynomial.polyder(derivative, second, axis=1)
            expected = np.polynomial.polynomial.polyval2d(sample_u, sample_tau, derivative)
            np.testing.assert_allclose(actual[first, second], expected, rtol=2e-5, atol=2e-5)
        for weight in (.01, .1):
            with self.assertRaises(RankDeficientFit) as failure:
                fit_thermodynamic_constraints(low, grid_u, grid_tau, weight)
            self.assertEqual(failure.exception.diagnostics['rank'], 322)
            self.assertEqual(failure.exception.diagnostics['unknowns'], 360)

    def test_third_continuity_but_not_fourth(self):
        grid = np.linspace(-.4, .6, 11)
        face = grid[5]
        left, right = np.nextafter(face, -np.inf), np.nextafter(face, np.inf)
        for order in range(4):
            np.testing.assert_allclose(spline_basis(grid, [left], order), spline_basis(grid, [right], order),
                                       atol=3e-10 / .1**order, rtol=1e-11)
        fourth_gap = spline_basis(grid, [right], 4) - spline_basis(grid, [left], 4)
        self.assertGreater(np.max(abs(fourth_gap)) * .1**4, 1.)

    def test_physical_response_one_sided_continuity(self):
        grid_u, grid_tau, polynomial, _ = polynomial_data()
        fitted = polynomial_coefficients(polynomial, grid_u, grid_tau)
        ln10 = np.log(10)
        face = grid_u[4]
        query_t = np.array([.123])
        traces = []
        for coordinate in [np.nextafter(face, -np.inf), np.nextafter(face, np.inf)]:
            native = evaluate_fit(fitted, grid_u, grid_tau, np.array([coordinate]), query_t)
            traces.append(thermodynamic_partials(native, np.array([coordinate / ln10]), query_t / ln10))
        for field in ['P', 'E', 'S', 'Cv', 'Cp', 'chiRho', 'chiT']:
            for suffix in ['', '_x', '_t']:
                values = [trace[field + suffix] for trace in traces]
                np.testing.assert_allclose(values[0], values[1], atol=1e-10, rtol=1e-10)

    def test_six_combinations_and_unknown_bound(self):
        grid_u, grid_tau, _, low = polynomial_data()
        combinations = thermodynamic_combinations(low)
        self.assertEqual(len(combinations), 6)
        np.testing.assert_allclose(combinations['TCv'], low[0, 1] - 2 * low[1, 0] - low[0, 2]
                                  + 4 * low[1, 1] - 4 * low[2, 0], atol=1e-14)
        with self.assertRaisesRegex(ValueError, 'QR bound'):
            fit_thermodynamic_constraints(low, grid_u, grid_tau, .1, max_unknowns=10)

    def test_nonzero_tensor_null_mode_changes_offnode_not_six_node_constraints(self):
        grid_u, grid_tau, _, _ = polynomial_data()
        witness, dimensions = null_witness(grid_u, grid_tau)
        self.assertEqual(dimensions, [2, 2])
        mesh_u, mesh_tau = np.meshgrid(grid_u, grid_tau, indexing='ij')
        nodal = evaluate_fit(witness, grid_u, grid_tau, mesh_u.ravel(), mesh_tau.ravel())
        spacing_u, spacing_tau = np.diff(grid_u)[0], np.diff(grid_tau)[0]
        for first, second in LOW_ORDERS:
            np.testing.assert_allclose(nodal[first, second] * spacing_u**first * spacing_tau**second, 0, atol=1e-12)
        offnode = evaluate_fit(witness, grid_u, grid_tau, np.array([-.27]), np.array([-.35]))
        self.assertGreater(abs(offnode[0, 0, 0]), 1e-8)


if __name__ == '__main__':
    unittest.main()
