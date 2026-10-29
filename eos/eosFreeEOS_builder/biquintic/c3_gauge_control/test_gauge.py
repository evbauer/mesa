"""Exact gauge physics, input roundtrip and explicit noninvariant roughness checks."""

import unittest

import numpy as np

from gauge import (
    EVAL_ORDERS, LOW_ORDERS, analytic_gauge_roughness, gauge_native,
    remove_input_gauge, shifted_low, thermodynamic_partials,
)


class GaugeTests(unittest.TestCase):
    def test_known_a_bT_thermodynamic_transform(self):
        sample_q, sample_t = np.array([-2., 1.]), np.array([5.1, 6.7])
        ln10 = np.log(10)
        thermal = 1e8 * 10.0**sample_t
        value = thermal * (sample_q * ln10 + .5 * sample_t * ln10 - 40)
        native = np.full((4, 4, len(sample_q)), np.nan)
        for first, second in EVAL_ORDERS:
            if first == 0:
                native[first, second] = value + .5 * second * thermal
            elif first == 1:
                native[first, second] = thermal
            else:
                native[first, second] = 0
        shift_a, shift_b = 7e11, 3e8
        addition = gauge_native(sample_t, shift_a, shift_b)
        changed = native.copy()
        for orders in EVAL_ORDERS:
            changed[orders] += addition[orders]
        base = thermodynamic_partials(native, sample_q, sample_t)
        altered = thermodynamic_partials(changed, sample_q, sample_t)
        np.testing.assert_allclose(altered['E'], base['E'] + shift_a, rtol=1e-12)
        np.testing.assert_allclose(altered['S'], base['S'] - shift_b, rtol=1e-12)
        for field in ['P', 'Pgas', 'Cv', 'Cp', 'chiRho', 'chiT']:
            for suffix in ['', '_x', '_t']:
                np.testing.assert_allclose(altered[field + suffix], base[field + suffix], rtol=1e-11, atol=1e-8)

    def test_low_gauge_only_pure_temperature_orders_and_roundtrip(self):
        low = np.full((3, 3, 7, 9), np.nan)
        for orders in LOW_ORDERS:
            low[orders] = np.random.default_rng(sum(orders) + 34).uniform(1e10, 1e11, (7, 9))
        logt = np.linspace(3.5, 4.3, 9)
        shifted = shifted_low(low, logt, 1e10, 1e8)
        restored = remove_input_gauge(shifted, logt, 1e10, 1e8)
        for orders in LOW_ORDERS:
            np.testing.assert_allclose(restored[orders], low[orders], rtol=1e-13)
            if orders[0] > 0:
                np.testing.assert_array_equal(shifted[orders], low[orders])

    def test_normal_fourth_roughness_is_not_a_bT_invariant(self):
        logq, logt = np.linspace(.35, 1.37, 18), np.linspace(6.5, 6.98, 13)
        constant = analytic_gauge_roughness(logq, logt, 0., 1e15)
        thermal = analytic_gauge_roughness(logq, logt, 1e9, 1e15)
        hundredfold = analytic_gauge_roughness(logq, logt, 1e11, 1e15)
        self.assertEqual(constant, 0)
        self.assertGreater(thermal, 0)
        np.testing.assert_allclose(hundredfold, 10000 * thermal, rtol=1e-14)


if __name__ == '__main__':
    unittest.main()
