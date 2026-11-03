"""Focused analytic controls for the 070 normalization and native endpoints."""

import math
import unittest

import numpy as np

from analyze import temperature_endpoints
from common import scaled_jump
from native import NativeHelm, REAL, SLOTS, decimal_cv


def polynomial_fields(density, temperature, density_power, temperature_power):
    result = []
    for density_order, temperature_order in SLOTS:
        if density_order > density_power or temperature_order > temperature_power:
            result.append(np.zeros((len(density), len(temperature))))
            continue
        factor = (math.factorial(density_power) / math.factorial(density_power - density_order)
                  * math.factorial(temperature_power) / math.factorial(temperature_power - temperature_order))
        result.append(np.asarray(factor * density[:, None]**(density_power - density_order)
                                 * temperature[None, :]**(temperature_power - temperature_order), dtype=float))
    return result


class MapTests(unittest.TestCase):
    def test_scale_invariant_and_all_zero(self):
        inputs = (np.array([2.0, 0]), np.array([3.0, 0]),
                  np.array([.4, 0]), np.array([-.2, 0]))
        original, denominator = scaled_jump(*inputs)
        self.assertEqual(original[1], 0)
        self.assertEqual(denominator[1], 0)
        for factor in (1e-20, 1e20):
            rescaled, _ = scaled_jump(*(factor * item for item in inputs))
            np.testing.assert_allclose(rescaled, original, rtol=2e-15)
        self.assertAlmostEqual(original[0], .2)

    def test_nonfinite_is_failure(self):
        with self.assertRaises(ValueError):
            scaled_jump(1, 1, float("nan"), 0)
        with self.assertRaises(ValueError):
            scaled_jump(1, float("inf"), 0, 0)

    def test_unequal_linear_width_endpoint(self):
        table = NativeHelm.__new__(NativeHelm)
        table.density = np.array([1, 1.3], dtype=REAL)
        table.temperature = np.array([1.89, 2, 2.37], dtype=REAL)
        table.fields = polynomial_fields(table.density, table.temperature, 3, 6)
        cv, left, right, difference = temperature_endpoints(table, np.array([0]), 1, .85)
        decimal_left = decimal_cv(table, 1.15, 2, .85, (0, 0))
        decimal_right = decimal_cv(table, 1.15, 2, .85, (0, 1))
        self.assertAlmostEqual(float(cv[0]), float(decimal_left[0]), delta=1e-10)
        self.assertAlmostEqual(float(left[0]), float(decimal_left[2]), delta=2e-9)
        self.assertAlmostEqual(float(right[0]), float(decimal_right[2]), delta=2e-9)
        self.assertAlmostEqual(float(difference[0]), float(decimal_left[2] - decimal_right[2]), delta=2e-9)

    def test_density_face_temperature_partial_is_continuous(self):
        table = NativeHelm.__new__(NativeHelm)
        table.density = np.array([1, 1.2, 1.5], dtype=REAL)
        table.temperature = np.array([2, 2.31], dtype=REAL)
        table.fields = polynomial_fields(table.density, table.temperature, 4, 6)
        left = table.cv(REAL(1.2), REAL(2.13), .85, (0, 0))
        right = table.cv(REAL(1.2), REAL(2.13), .85, (1, 0))
        np.testing.assert_allclose(left, right, rtol=2e-11, atol=1e-8)


if __name__ == "__main__":
    unittest.main()
