"""Analytic checks independent of archived stellar results."""

import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

from native import REAL, NativeHelm, basis, decimal_cv, polynomial_corners, read_cache, response_jump, tensor_jet, thermal_jacobian


class NativeTests(unittest.TestCase):
    def test_endpoint_nodal_constraints(self):
        for width in (.017, 21.0):
            for endpoint in (0, 1):
                for order in range(3):
                    expected = np.zeros((2, 3))
                    expected[endpoint, order] = 1
                    np.testing.assert_allclose(basis(endpoint, width, order), expected, atol=1e-15)

    def test_tensor_polynomials_on_unequal_linear_widths(self):
        density_axis, temperature_axis = (1.0, 1.23), (2.0, 2.78)
        density, temperature = 1.113, 2.34
        for density_power, temperature_power in ((0, 0), (1, 2), (2, 3), (5, 5)):
            corners = polynomial_corners(density_axis, temperature_axis, density_power, temperature_power)
            for density_order, temperature_order in ((0, 2), (0, 3), (1, 2)):
                expected = 0.0
                if density_power >= density_order and temperature_power >= temperature_order:
                    expected = math.factorial(density_power) / math.factorial(density_power - density_order) * math.factorial(
                        temperature_power) / math.factorial(temperature_power - temperature_order) * density**(
                            density_power - density_order) * temperature**(temperature_power - temperature_order)
                value = tensor_jet(corners, .23, .78, (density - 1) / .23,
                                   (temperature - 2) / .78, density_order, temperature_order)
                self.assertAlmostEqual(float(value), expected, delta=max(1e-9, abs(expected)*2e-11))

    def test_nonuniform_sextic_endpoint_jump(self):
        face, left_width, right_width = REAL(2), REAL('.11'), REAL('.37')
        traces = []
        for axis, fraction in (((face-left_width, face), 1), ((face, face+right_width), 0)):
            corners = polynomial_corners((1, 1.3), axis, 0, 6)
            traces.append(tensor_jet(corners, REAL('.3'), axis[1]-axis[0], REAL('.41'), fraction, 0, 3))
        self.assertAlmostEqual(float(traces[1] - traces[0]), float(6*(left_width**3 + right_width**3)), delta=1e-9)

    def test_density_faces_preserve_requested_partials(self):
        generator = np.random.default_rng(48)
        nodes = generator.normal(size=(9, 3, 2))
        for orders in ((0, 2), (0, 3), (1, 2)):
            left = tensor_jet(nodes[:, :2], .19, .31, 1, .37, *orders)
            right = tensor_jet(nodes[:, 1:], .27, .31, 0, .37, *orders)
            self.assertAlmostEqual(float(left), float(right), delta=1e-10)

    def test_temperature_faces_preserve_mixed_partial(self):
        generator = np.random.default_rng(49)
        nodes = generator.normal(size=(9, 2, 3))
        for orders in ((0, 2), (1, 2)):
            left = tensor_jet(nodes[:, :, :2], .19, .31, .37, 1, *orders)
            right = tensor_jet(nodes[:, :, 1:], .19, .41, .37, 0, *orders)
            self.assertAlmostEqual(float(left), float(right), delta=1e-10)

    def test_thermal_jacobian_by_differentiation(self):
        current, previous, step = .3, .299, 1e-5
        for coefficient in (2.0, 4.0):
            def residual(temperature_log):
                return -np.exp(temperature_log) * np.exp(coefficient*temperature_log) * (temperature_log - previous)
            cv = np.exp(coefficient*current)
            state = [cv, 0.0, coefficient*cv]
            result = thermal_jacobian(state, state, current - previous)
            numerical = (residual(current+step) - residual(current-step))/(2*step)
            self.assertAlmostEqual(result['left_over_T']*np.exp(current), numerical, delta=2e-8)
        self.assertEqual(response_jump([1, 0, 0], [1, 0, .5])['Cv_t']['signed_scaled'], -.5)

    def test_bad_cache_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'bad.bin'
            path.write_bytes(b'wrong record')
            with self.assertRaises(ValueError):
                read_cache(path, 2, 2)

    def test_physical_log_coordinate_chain(self):
        table = NativeHelm.__new__(NativeHelm)
        table.density = np.array([1, 1.5], dtype=REAL)
        table.temperature = np.array([2, 3], dtype=REAL)
        table.fields = -polynomial_corners(table.density, table.temperature, 1, 3)
        density, temperature, ye = 1.25, 2.4, .85
        cv = 6*ye*density*temperature**2
        expected = [cv, cv, 2*cv]
        np.testing.assert_allclose(table.cv(density, temperature, ye), expected, rtol=1e-12)
        np.testing.assert_allclose([float(value) for value in decimal_cv(table, density, temperature, ye, (0, 0))], expected, rtol=1e-14)


if __name__ == '__main__':
    unittest.main()
