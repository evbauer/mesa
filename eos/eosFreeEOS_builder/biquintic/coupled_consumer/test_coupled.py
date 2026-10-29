"""Analytic and independent-column derivative checks for the coupled solve."""

import sys
import unittest
from pathlib import Path

import numpy as np


sys.path.insert(0, str(Path(__file__).resolve().parent))

from models import CoupledModels, LEGACY_GRAD_AD_COLUMN, ideal_state
from solver import Problem, reference, residual, solve
from gravity_control import gravity_residual
from freeeos_potential import LN10, read_table
from thermal_consumer.models import SOURCE, legacy


class CoupledTests(unittest.TestCase):
    def test_ideal_manufactured_root_and_all_jacobian_entries(self):
        initial = np.array([np.log(1e-4), np.log(1e6)])
        target = initial + np.array([.003, -.004])
        target_state = ideal_state(*target)
        heating = -gravity_residual(target_state, *(target - initial))["value"]
        problem = Problem(lambda x, t, cell=None: ideal_state(x, t), *initial,
                          target_state["P"], heating, np.exp(initial[1]) * 1.5e8,
                          (initial[0] - .15, initial[0] + .15),
                          (initial[1] - .04, initial[1] + .04))
        values, jacobian, _ = residual(problem, *target)
        np.testing.assert_allclose(values, 0, atol=1e-14)
        for axis in range(2):
            for step in (1e-4, 1e-5, 1e-6):
                offset = np.zeros(2)
                offset[axis] = step
                finite = (residual(problem, *(target + offset))[0]
                          - residual(problem, *(target - offset))[0]) / (2 * step)
                np.testing.assert_allclose(finite, jacobian[:, axis], rtol=2e-8, atol=2e-9)
        independent = reference(problem)
        self.assertEqual(independent["status"], "bracketed")
        np.testing.assert_allclose(independent["root"], target, atol=1e-12, rtol=0)
        trial = solve(problem, initial)
        self.assertEqual(trial["status"], "converged")
        np.testing.assert_allclose(trial["root"], independent["root"], atol=1e-9, rtol=0)

    def test_legacy_total_pressure_and_gravity_partials_are_numerical(self):
        models = CoupledModels()
        density_log = .0968001524174138 * np.log(10)
        temperature_log = 6.368400076208707 * np.log(10)
        self.assertTrue(models.valid(density_log, temperature_log))
        state = models.at("legacy", density_log, temperature_log)
        for suffix, axis in (("_x", 0), ("_t", 1)):
            step = 1e-5
            offset = np.array([step, 0]) if axis == 0 else np.array([0, step])
            plus = models.at("legacy", density_log + offset[0], temperature_log + offset[1])
            minus = models.at("legacy", density_log - offset[0], temperature_log - offset[1])
            for field in ("P", "Cp", "chiRho", "chiT", "grad_ad"):
                finite = (plus[field] - minus[field]) / (2 * step)
                self.assertAlmostEqual(state[field + suffix] / max(1, abs(finite)),
                                       finite / max(1, abs(finite)), delta=2e-5)

    def test_legacy_gradient_is_independent_cache_column(self):
        models = CoupledModels()
        _, _, _, source_rows = read_table(SOURCE)
        self.assertEqual(models.tables.cache.shape[1], 27)
        self.assertEqual(LEGACY_GRAD_AD_COLUMN, 6)
        for density_index, temperature_index in ((20, 50), (50, 100), (100, 150)):
            self.assertEqual(models.tables.cache[0, LEGACY_GRAD_AD_COLUMN,
                                                density_index, temperature_index],
                             source_rows[density_index, temperature_index, 15])
        density_log = .0968001524174138 * LN10
        temperature_log = 6.368400076208707 * LN10
        sample_t = np.array([temperature_log / LN10])
        sample_q = np.array([density_log / LN10 - 2 * sample_t[0] + 12])
        state = models.at("legacy", density_log, temperature_log)
        values = legacy.evaluate_spline(models.tables.cache, models.tables.logq,
                                        models.tables.logt, sample_q, sample_t)
        native_q = legacy.evaluate_spline(models.tables.cache, models.tables.logq,
                                          models.tables.logt, sample_q, sample_t, (1, 0))
        native_t = legacy.evaluate_spline(models.tables.cache, models.tables.logq,
                                          models.tables.logt, sample_q, sample_t, (0, 1))
        column = LEGACY_GRAD_AD_COLUMN
        self.assertEqual(state["grad_ad"], values[0, column])
        self.assertAlmostEqual(state["grad_ad_x"], native_q[0, column] / LN10, places=12)
        self.assertAlmostEqual(state["grad_ad_t"],
                               (native_t[0, column] - 2 * native_q[0, column]) / LN10,
                               places=12)

    def test_reference_reports_unbracketed_pressure(self):
        initial = np.array([np.log(1e-4), np.log(1e6)])
        problem = Problem(lambda x, t, cell=None: ideal_state(x, t), *initial,
                          1e200, 0., np.exp(initial[1]) * 1.5e8,
                          (initial[0] - .1, initial[0] + .1),
                          (initial[1] - .01, initial[1] + .01))
        outcome = reference(problem)
        self.assertEqual(outcome["status"], "reference_failure")
        self.assertEqual(outcome["reason"], "pressure_not_bracketed")


if __name__ == "__main__":
    unittest.main()
