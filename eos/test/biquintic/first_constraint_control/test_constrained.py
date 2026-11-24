"""Exact equality, explicit roughness stationarity and nonzero soft-objective controls."""

import unittest

import numpy as np
from scipy.sparse import csr_matrix

from constrained import HardConstraintFailure, solve_constrained, solve_hard_system
from spline_control import evaluate_fit, thermodynamic_combinations
from freeeos_potential import LOW_ORDERS


def analytic_low():
    grid_u, grid_tau = np.linspace(-.3, .5, 9), np.linspace(-.4, .6, 8)
    mesh_u, mesh_tau = np.meshgrid(grid_u, grid_tau, indexing='ij')
    polynomial = np.zeros((4, 4))
    polynomial[0, 0], polynomial[1, 0], polynomial[0, 1], polynomial[2, 2], polynomial[3, 1] = 7., 1.2, -2., .9, .3
    low = np.full((3, 3, len(grid_u), len(grid_tau)), np.nan)
    for first, second in LOW_ORDERS:
        derivative = np.polynomial.polynomial.polyder(polynomial, first, axis=0)
        derivative = np.polynomial.polynomial.polyder(derivative, second, axis=1)
        low[first, second] = np.polynomial.polynomial.polyval2d(mesh_u, mesh_tau, derivative)
    return low, grid_u, grid_tau


class ConstrainedTests(unittest.TestCase):
    def test_analytic_exact_first_and_zero_bicubic_objective(self):
        low, grid_u, grid_tau = analytic_low()
        coefficients, solve = solve_constrained(low, grid_u, grid_tau, .01)
        self.assertEqual(solve['hard_rank'], 216)
        self.assertEqual(solve['reduced_rank'], solve['nullspace_columns'])
        self.assertLess(solve['total_objective'], 1e-16)
        mesh_u, mesh_tau = np.meshgrid(grid_u, grid_tau, indexing='ij')
        actual = thermodynamic_combinations(evaluate_fit(coefficients, grid_u, grid_tau, mesh_u.ravel(), mesh_tau.ravel()))
        targets = thermodynamic_combinations(low)
        for field in ['E', 'TS', 'P_over_rho']:
            np.testing.assert_allclose(actual[field], targets[field].ravel(), rtol=1e-10, atol=1e-10)

    def test_incompatible_soft_responses_leave_first_exact_nonzero_objective(self):
        low, grid_u, grid_tau = analytic_low()
        changed = low.copy()
        changed[2, 0] += np.random.default_rng(33).normal(size=low.shape[2:]) * .01
        coefficients, solve = solve_constrained(changed, grid_u, grid_tau, .01)
        self.assertGreater(solve['total_objective'], 1e-10)
        self.assertLess(solve['hard_relative_residual_max'], 1e-9)
        mesh_u, mesh_tau = np.meshgrid(grid_u, grid_tau, indexing='ij')
        actual = thermodynamic_combinations(evaluate_fit(coefficients, grid_u, grid_tau, mesh_u.ravel(), mesh_tau.ravel()))
        targets = thermodynamic_combinations(changed)
        for field in ['E', 'TS', 'P_over_rho']:
            np.testing.assert_allclose(actual[field], targets[field].ravel(), rtol=1e-10, atol=1e-10)

    def test_inconsistent_dependent_hard_rows_fail_not_minimum_norm(self):
        hard = np.array([[1., 0., 0.], [1., 0., 0.]])
        with self.assertRaises(HardConstraintFailure):
            solve_hard_system(hard, np.array([1., 2.]), csr_matrix(np.eye(3)), np.zeros(3))


if __name__ == '__main__':
    unittest.main()
