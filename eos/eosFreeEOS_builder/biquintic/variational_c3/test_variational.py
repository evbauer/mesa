"""Analytic minimization, validity isolation, source-jet and shared-face controls."""

import unittest

import numpy as np

from normal_jets import (
    JET_ORDERS, LOW_ORDERS, STIFFNESS, coarse_bridge_guard, construct_variational_jets,
    line_energy, line_system, solve_line, solve_normal,
)
from c3_potential import evaluate
from freeeos_potential import thermodynamic_partials


class VariationalTests(unittest.TestCase):
    def test_original_invalid_midpoints_break_coarse_runs(self):
        full_valid = np.ones((19, 17), dtype=bool)
        full_valid[7, :] = False
        full_valid[:, 11] = False
        coarse = full_valid[::2, ::2]
        guard, evidence = coarse_bridge_guard(full_valid, coarse)
        self.assertFalse(guard[3:5, :].any())
        self.assertFalse(guard[:, 5:7].any())
        self.assertEqual(evidence['axis0_invalid_original_midpoints'], 9)
        self.assertEqual(evidence['axis1_invalid_original_midpoints'], 10)

    def test_quadratic_and_cubic_reproduction(self):
        coordinate = np.linspace(-1, 2, 31)
        for cubic in [0., .7]:
            value = 2 + .8 * coordinate - .3 * coordinate**2 + cubic * coordinate**3
            first = .8 - .6 * coordinate + 3 * cubic * coordinate**2
            second = -.6 + 6 * cubic * coordinate
            third, stationarity = solve_line(value, first, second, .1)
            np.testing.assert_allclose(third, 6 * cubic, atol=2e-10)
            self.assertLess(stationarity, 1e-12)
            self.assertLess(line_energy(value, first, second, third, .1), 1e-16)

    def test_positive_definite_minimum_energy(self):
        generator = np.random.default_rng(25)
        value, first, second = generator.normal(size=(3, 19))
        spacing = .17
        third, stationarity = solve_line(value, first, second, spacing)
        diagonal, off_diagonal, rhs = line_system(value, first, second, spacing)
        matrix = np.diag(diagonal) + np.diag(off_diagonal, 1) + np.diag(off_diagonal, -1)
        self.assertGreater(np.linalg.eigvalsh(matrix).min(), 0)
        np.testing.assert_allclose(matrix @ (spacing**3 * third), rhs, rtol=1e-12, atol=1e-10)
        self.assertLess(stationarity, 1e-12)
        minimum = line_energy(value, first, second, third, spacing)
        for scale in [.1, 1., 10.]:
            competitor = third + scale * generator.normal(size=19) / spacing**3
            self.assertLessEqual(minimum, line_energy(value, first, second, competitor, spacing) * (1 + 1e-12))
        np.testing.assert_array_equal(STIFFNESS[np.ix_([3, 7], [3, 7])], [[16., 4.], [4., 16.]])

    def test_no_bridge_and_nonlocal_within_segment(self):
        coordinate = np.arange(25) * .1
        value = np.exp(coordinate)[:, None]
        valid = np.ones_like(value, dtype=bool)
        valid[12] = False
        original, mask, segments = solve_normal(value, value, value, valid, .1, 0, 2)
        changed = value.copy()
        changed[15] += 1
        changed[12] = np.nan
        perturbed, perturbed_mask, _ = solve_normal(changed, value, value, valid, .1, 0, 2)
        np.testing.assert_array_equal(original[:12], perturbed[:12])
        np.testing.assert_array_equal(mask, perturbed_mask)
        self.assertFalse(mask[12, 0])
        self.assertNotEqual(original[-1, 0], perturbed[-1, 0])
        self.assertEqual([(segment['start'], segment['stop_exclusive']) for segment in segments], [(0, 12), (13, 25)])

    def test_six_constraints_and_c3_faces(self):
        logq, logt = np.linspace(-2, 1, 35), np.linspace(3.5, 5, 33)
        du, dtau = np.diff(logq)[0] * np.log(10), np.diff(logt)[0] * np.log(10)
        value = np.exp(.3 * logq[:, None] * np.log(10) + .4 * logt[None, :] * np.log(10))
        low = np.full((3, 3, len(logq), len(logt)), np.nan)
        for first, second in LOW_ORDERS:
            low[first, second] = value * .3**first * .4**second
        nodes, complete, _ = construct_variational_jets(low, np.ones(value.shape, bool), du, dtau)
        for orders in LOW_ORDERS:
            np.testing.assert_array_equal(nodes[orders][complete], low[orders][complete])
        self.assertTrue(complete[16:19, 15:18].all())
        for axis in [0, 1]:
            sample_q = np.array([logq[17] if axis == 0 else (logq[16] + logq[17]) / 2])
            sample_t = np.array([(logt[16] + logt[17]) / 2 if axis == 0 else logt[17]])
            left = (np.array([16]), np.array([16]))
            right = (np.array([17 if axis == 0 else 16]), np.array([16 if axis == 0 else 17]))
            left_trace = evaluate(nodes, complete, logq, logt, sample_q, sample_t, JET_ORDERS, left)
            right_trace = evaluate(nodes, complete, logq, logt, sample_q, sample_t, JET_ORDERS, right)
            np.testing.assert_allclose(left_trace, right_trace, atol=1e-12, rtol=1e-12)

    def test_analytic_physical_jacobians_of_returned_potential(self):
        logq, logt = np.linspace(-4, 1, 39), np.linspace(3.5, 5, 37)
        ln10 = np.log(10)
        coordinate_q, coordinate_t = logq[:, None] * ln10, logt[None, :] * ln10
        thermal = np.broadcast_to(1e8 * np.exp(coordinate_t), (len(logq), len(logt)))
        low = np.full((3, 3, *thermal.shape), np.nan)
        low[0, 0] = thermal * (coordinate_q + .5 * coordinate_t - 30)
        low[1, 0], low[2, 0], low[1, 1] = thermal, 0., thermal
        low[0, 1], low[0, 2] = low[0, 0] + .5 * thermal, low[0, 0] + thermal
        nodes, complete, _ = construct_variational_jets(low, np.ones(thermal.shape, bool),
                                                       np.diff(logq)[0] * ln10, np.diff(logt)[0] * ln10)
        sample_q, sample_t = np.array([-1.43]), np.array([4.237])

        def state(query_q, query_t):
            return thermodynamic_partials(evaluate(nodes, complete, logq, logt, query_q, query_t), query_q, query_t)

        center = state(sample_q, sample_t)
        self.assertGreater(center['Cv'][0], 0)
        step = 1e-4
        for suffix, shift_q, shift_t in [('_x', step / ln10, 0), ('_t', -2 * step / ln10, step / ln10)]:
            plus, minus = state(sample_q + shift_q, sample_t + shift_t), state(sample_q - shift_q, sample_t - shift_t)
            for field in ['P', 'E', 'S', 'Cv', 'Cp', 'chiRho', 'chiT']:
                numerical = (plus[field] - minus[field]) / (2 * step)
                scale = np.maximum(1, np.maximum(abs(center[field]), abs(center[field + suffix])))
                np.testing.assert_allclose((numerical - center[field + suffix]) / scale, 0, atol=1e-6)


if __name__ == '__main__':
    unittest.main()
