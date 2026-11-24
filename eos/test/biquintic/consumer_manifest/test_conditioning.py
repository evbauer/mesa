"""Cancellation identities and gauge dependence, without source-table truth."""

import unittest

import numpy as np

from conditioning import indicators


class ConditioningTests(unittest.TestCase):
    def test_ideal_gas_and_entropy_gauge(self):
        temperature = np.array([1e5])
        heat_capacity = 2e8
        entropy = 12 * heat_capacity
        physical = np.zeros((3, 3, 1))
        physical[0, 1] = -temperature * entropy
        physical[0, 2] = -temperature * (entropy + heat_capacity)
        physical[0, 0] = temperature * heat_capacity + physical[0, 1]
        result = indicators(physical, temperature)
        self.assertAlmostEqual(result["kappa_f"][0], 25.)
        self.assertAlmostEqual(result["kappa_f_over_T"][0], 1.)
        shifted = physical.copy()
        for order in (0, 1, 2):
            shifted[0, order] += 10 * heat_capacity * temperature
        shifted_result = indicators(shifted, temperature)
        self.assertAlmostEqual(shifted_result["Cv"][0], heat_capacity)
        self.assertAlmostEqual(shifted_result["kappa_f"][0], 5.)
        self.assertAlmostEqual(shifted_result["kappa_f_over_T"][0], 1.)

    def test_energy_gauge_changes_quotient_conditioning(self):
        temperature = np.array([1e5])
        physical = np.zeros((3, 3, 1))
        physical[0, 0] = temperature * 2e8
        physical[0, 2] = -temperature * 2e8
        before = indicators(physical, temperature)
        physical[0, 0] += 100 * temperature * 2e8
        after = indicators(physical, temperature)
        self.assertAlmostEqual(before["kappa_f"][0], after["kappa_f"][0])
        self.assertAlmostEqual(after["kappa_f_over_T"][0], 201.)
        self.assertAlmostEqual(before["Cv"][0], after["Cv"][0])

    def test_invalid_heat_capacity_rejected(self):
        with self.assertRaises(ValueError):
            indicators(np.zeros((3, 3, 1)), np.array([1e5]))


if __name__ == "__main__":
    unittest.main()
