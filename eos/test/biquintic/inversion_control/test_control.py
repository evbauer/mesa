"""Small exact-function tests of inversion safeguards, not stellar EOS tests."""

import unittest

import numpy as np

from solve import solve


class InversionTests(unittest.TestCase):
    def test_analytic_newton_and_escape(self):
        def callback(temperature_log, indices):
            temperature = np.exp(temperature_log)
            return temperature, temperature, np.ones_like(temperature)
        target = np.array([1.0])
        lower, upper = np.array([-.7]), np.array([.18])
        raw = solve(callback, target, lower, upper, lower, False, rtol=1e-12)
        self.assertEqual(raw['status'][0], 'left_declared_bracket')
        safe = solve(callback, target, lower, upper, lower, True, rtol=1e-12)
        self.assertEqual(safe['status'][0], 'converged')
        self.assertLess(abs(safe['lnT'][0]), 1e-12)
        self.assertGreater(safe['safeguards'][0], 0)

    def test_unsafe_cv_bisection(self):
        def callback(temperature_log, indices):
            return np.exp(temperature_log), -np.ones_like(temperature_log), -np.ones_like(temperature_log)
        lower, upper = np.array([-1.]), np.array([1.])
        target = np.array([1.0])
        safe = solve(callback, target, lower, upper, lower, True)
        self.assertEqual(safe['status'][0], 'converged')
        self.assertGreater(safe['unsafe_Cv'][0], 0)
        raw = solve(callback, target, lower, upper, lower, False)
        self.assertEqual(raw['status'][0], 'unsafe_Cv')

    def test_invalid_bracket(self):
        def callback(temperature_log, indices):
            temperature = np.exp(temperature_log)
            return temperature, temperature, np.ones_like(temperature)
        result = solve(callback, np.array([10.]), np.array([0.]), np.array([1.]), np.array([.5]), True)
        self.assertEqual(result['status'][0], 'invalid_bracket')

    def test_batch_analytic_radiation_jacobian(self):
        rho = np.array([1e-4, 1e-2, 1.])
        truth = np.log(np.array([1e4, 2e4, 3e4]))
        def callback(temperature_log, indices):
            temperature = np.exp(temperature_log)
            radiation = 7.56573325028e-15 * temperature**4 / rho[indices]
            energy = 1.5e8 * temperature + radiation
            jacobian = 1.5e8 * temperature + 4 * radiation
            return energy, jacobian, jacobian / temperature
        indices = np.arange(3)
        target, jacobian, _ = callback(truth, indices)
        step = 1e-5
        numerical = (callback(truth + step, indices)[0] - callback(truth - step, indices)[0]) / (2 * step)
        np.testing.assert_allclose(jacobian, numerical, rtol=1e-8)
        result = solve(callback, target, truth - .7, truth + .2, truth + .2, True, rtol=1e-12)
        self.assertTrue(np.all(result['status'] == 'converged'))
        np.testing.assert_allclose(result['lnT'], truth, rtol=0, atol=1e-12)


if __name__ == '__main__':
    unittest.main()
