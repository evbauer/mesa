"""Standalone analytic, coordinate, and stencil-isolation tests."""

import unittest
import tempfile
from pathlib import Path

import numpy as np

from freeeos_potential import (
    CRAD, EVAL_ORDERS, LN10, LOW_ORDERS, construct_twists, evaluate,
    identity_residuals, local_derivative, radiation_potential, read_table, source_constraints, thermodynamics,
    thermodynamic_partials, to_native, to_physical, valid_cells,
)


class PreprocessingTests(unittest.TestCase):
    def test_strict_source_header_validation(self):
        header = [51., .7, .02, 2., 3.2, 3.22, .02, 2., 0., .03, .03]
        rows = np.zeros((2, 2, 22))
        rows[:, :, 0] = [3.2, 3.22]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "source.data"

            def write_source(values):
                with path.open("w") as stream:
                    stream.write("test source\n" + " ".join(map(str, values)) + "\n")
                    np.savetxt(stream, rows.reshape(-1, 22), fmt="%.17g")

            write_source(header)
            actual, logq, logt, actual_rows = read_table(path)
            np.testing.assert_array_equal(actual, header)
            np.testing.assert_allclose(logq, [0., .03])
            np.testing.assert_allclose(logt, [3.2, 3.22])
            np.testing.assert_array_equal(actual_rows, rows)
            for index, value in ((0, 51.5), (3, 2.5), (6, 0.), (9, 2.03), (1, .99)):
                altered = header.copy()
                altered[index] = value
                write_source(altered)
                with self.assertRaises(ValueError):
                    read_table(path)
            altered = header.copy()
            altered[9] = 2.03
            write_source(altered)
            _, logq, _, _ = read_table(path, validate_endpoints=False)
            np.testing.assert_allclose(logq, [0., .03])

    def test_sheared_roundtrip(self):
        generator = np.random.default_rng(20261028)
        physical = np.full((3, 3, 9), np.nan)
        for orders in LOW_ORDERS:
            physical[orders] = generator.normal(size=9)
        restored = to_physical(to_native(physical))
        for orders in LOW_ORDERS:
            np.testing.assert_allclose(restored[orders], physical[orders], atol=2e-15)

    def test_local_derivative_polynomial(self):
        coordinate = np.arange(29) * .13
        valid = np.ones((29, 11), dtype=bool)
        for radius in (1, 2, 3):
            power = 2 * radius
            values = np.broadcast_to(coordinate[:, None]**power, valid.shape)
            actual, accepted = local_derivative(values, valid, .13, 0, radius)
            expected = np.broadcast_to(power * coordinate[:, None]**(power - 1), valid.shape)
            np.testing.assert_allclose(actual[accepted], expected[accepted], rtol=5e-13, atol=1e-11)
            self.assertFalse(np.any(accepted[:radius]))
            self.assertFalse(np.any(accepted[-radius:]))

    def test_fallback_and_low_temperature_isolation(self):
        generator = np.random.default_rng(7)
        low = generator.normal(size=(3, 3, 31, 31))
        valid = np.ones((31, 31), dtype=bool)
        valid[16, 16] = False
        valid[:, :5] = False
        contaminated = low.copy()
        contaminated[:, :, ~valid] = 1e200
        nodes, accepted, _ = construct_twists(low, valid, .1, .12, radius=2)
        altered, accepted_altered, _ = construct_twists(contaminated, valid, .1, .12, radius=2)
        np.testing.assert_array_equal(accepted, accepted_altered)
        np.testing.assert_array_equal(nodes[:, :, accepted], altered[:, :, accepted])
        self.assertFalse(accepted[16, 16])
        self.assertFalse(np.any(accepted[:, :9]))
        self.assertTrue(np.any(accepted))

    def test_polynomial_hermite_and_shear(self):
        logq, logt = np.array([-.5, -.17]), np.array([.1, .29])
        native = np.zeros((3, 3, 2, 2))
        generator = np.random.default_rng(3)
        coefficients = generator.normal(size=(6, 6))
        for density_order in range(3):
            for temperature_order in range(3):
                polynomial = np.polynomial.polynomial.polyder(coefficients, density_order, axis=0)
                polynomial = np.polynomial.polynomial.polyder(polynomial, temperature_order, axis=1)
                for density_end in range(2):
                    for temperature_end in range(2):
                        native[density_order, temperature_order, density_end, temperature_end] = np.polynomial.polynomial.polyval2d(
                            logq[density_end] * LN10, logt[temperature_end] * LN10, polynomial)
        sample_q = generator.uniform(logq[0], logq[-1], 20)
        sample_t = generator.uniform(logt[0], logt[-1], 20)
        actual = evaluate(native, np.ones((2, 2), bool), logq, logt, sample_q, sample_t)
        for orders in EVAL_ORDERS:
            polynomial = np.polynomial.polynomial.polyder(coefficients, orders[0], axis=0)
            polynomial = np.polynomial.polynomial.polyder(polynomial, orders[1], axis=1)
            expected = np.polynomial.polynomial.polyval2d(sample_q * LN10, sample_t * LN10, polynomial)
            np.testing.assert_allclose(actual[orders], expected, rtol=5e-10, atol=2e-10)
        physical = to_physical(actual, EVAL_ORDERS)
        np.testing.assert_allclose(physical[1, 1], actual[1, 1] - 2 * actual[2, 0])
        np.testing.assert_allclose(physical[0, 3], actual[0, 3] - 6 * actual[1, 2]
                                   + 12 * actual[2, 1] - 8 * actual[3, 0])

    def test_source_constraints_ideal_gas_and_radiation(self):
        logq, logt = np.linspace(-2, 0, 9), np.linspace(3.2, 4.2, 9)
        temperature = np.broadcast_to(10**logt[None, :], (9, 9))
        density = 10**(logq[:, None] + 2 * logt[None, :] - 12)
        logdensity, logtemperature = np.log(density), np.log(temperature)
        gas_constant, heat_factor = 1.2e8, 1.5
        radiation_free = -CRAD * temperature**4 / (3 * density)
        physical = np.full((4, 4, 9, 9), np.nan)
        for density_order, temperature_order in EVAL_ORDERS:
            if density_order == 0:
                gas = gas_constant * temperature * (logdensity - heat_factor * logtemperature - heat_factor * temperature_order)
            elif density_order == 1:
                gas = gas_constant * temperature
            else:
                gas = np.zeros((9, 9))
            physical[density_order, temperature_order] = gas + (-1)**density_order * 4**temperature_order * radiation_free
        state = thermodynamics(to_native(physical), logq[:, None], logt[None, :])
        rows = np.zeros((9, 9, 22))
        rows[:, :, 18] = np.log10(density)
        rows[:, :, 1] = np.log10(state["Pgas"])
        rows[:, :, 2] = np.log10(state["E"])
        rows[:, :, 3] = np.log10(state["S"])
        rows[:, :, 4] = state["chiRho"]
        rows[:, :, 5] = state["chiT"]
        rows[:, :, 7] = state["Cv"]
        actual, valid, _ = source_constraints(rows, logq, logt)
        self.assertTrue(np.all(valid))
        expected = to_native(physical)
        for orders in LOW_ORDERS:
            np.testing.assert_allclose(actual[orders], expected[orders], rtol=1e-12, atol=.1)
        split, split_valid, _ = source_constraints(rows, logq, logt, split_radiation=True)
        np.testing.assert_array_equal(split_valid, valid)
        radiation = radiation_potential(logq[:, None], logt[None, :], LOW_ORDERS)
        for orders in LOW_ORDERS:
            np.testing.assert_allclose(split[orders] + radiation[orders], expected[orders], rtol=1e-12, atol=.1)
        reconstructed = thermodynamics(actual, logq[:, None], logt[None, :])
        for name in ("P", "E", "S", "Cv", "chiRho", "chiT"):
            np.testing.assert_allclose(reconstructed[name], state[name], rtol=5e-13)
        for residual in identity_residuals(reconstructed).values():
            self.assertLess(np.max(abs(residual)), 2e-15)
        damaged = rows.copy()
        damaged[4, 4, 18] += 1e-4
        damaged[5, 5, 17] = 1
        damaged[6, 6, 7] = -1
        damaged[3, 3, 8] = np.nan
        _, damaged_valid, _ = source_constraints(damaged, logq, logt)
        self.assertFalse(damaged_valid[4, 4])
        self.assertFalse(damaged_valid[5, 5])
        self.assertFalse(damaged_valid[6, 6])
        self.assertFalse(damaged_valid[3, 3])
        shifted_logt = logt - .2
        _, cold_valid, _ = source_constraints(rows, logq, shifted_logt)
        self.assertFalse(np.any(cold_valid[:, :2]))

    def test_thermodynamic_jacobians(self):
        def analytic_state(logdensity, logtemperature):
            temperature = np.exp(logtemperature)
            density = np.exp(logdensity)
            radiation_free = -CRAD * temperature**4 / (3 * density)
            physical = np.full((4, 4, len(logdensity)), np.nan)
            for density_order, temperature_order in EVAL_ORDERS:
                if density_order == 0:
                    gas = 1.2e8 * temperature * (logdensity - 1.5 * logtemperature - 1.5 * temperature_order)
                elif density_order == 1:
                    gas = 1.2e8 * temperature
                else:
                    gas = np.zeros_like(temperature)
                physical[density_order, temperature_order] = gas + (-1)**density_order * 4**temperature_order * radiation_free
            native = to_native(physical, EVAL_ORDERS)
            logq = (logdensity - 2 * logtemperature) / LN10 + 12
            logt = logtemperature / LN10
            return thermodynamic_partials(native, logq, logt)

        logdensity = np.array([-12., -7., 0.])
        logtemperature = np.array([12., 14., 16.])
        state = analytic_state(logdensity, logtemperature)
        step = 1e-4
        for axis, suffix in ((0, "_x"), (1, "_t")):
            shifted_density, shifted_temperature = logdensity.copy(), logtemperature.copy()
            if axis == 0:
                shifted_density += step
            else:
                shifted_temperature += step
            plus = analytic_state(shifted_density, shifted_temperature)
            if axis == 0:
                shifted_density -= 2 * step
            else:
                shifted_temperature -= 2 * step
            minus = analytic_state(shifted_density, shifted_temperature)
            for field in ("P", "Pgas", "E", "S", "Cv", "Cp", "chiRho", "chiT", "gamma1", "gamma3", "grad_ad"):
                expected = state[field + suffix]
                actual = (plus[field] - minus[field]) / (2 * step)
                scale = np.maximum.reduce([np.ones_like(expected), abs(expected), abs(state[field])])
                self.assertLess(np.max(abs(actual - expected) / scale), 5e-8, (field, suffix))

    def test_query_rejection(self):
        nodes = np.ones((3, 3, 3, 3))
        complete = np.ones((3, 3), dtype=bool)
        grid = np.arange(3.)
        for query in (np.nan, -.1, 2.1):
            with self.assertRaises(ValueError):
                evaluate(nodes, complete, grid, grid, np.array([query]), np.array([.5]))
        complete[1, 1] = False
        self.assertFalse(np.any(valid_cells(complete)))
        with self.assertRaises(ValueError):
            evaluate(nodes, complete, grid, grid, np.array([.5]), np.array([.5]))

    def test_analytic_radiation_addition(self):
        logq, logt = np.array([-2., -1.9]), np.array([4., 4.1])
        nodes = np.zeros((3, 3, 2, 2))
        query_q, query_t = np.array([-1.96]), np.array([4.04])
        actual = evaluate(nodes, np.ones((2, 2), bool), logq, logt, query_q, query_t, add_radiation=True)
        radiation = radiation_potential(query_q, query_t)
        for orders in EVAL_ORDERS:
            np.testing.assert_array_equal(actual[orders], radiation[orders])
        physical = to_physical(actual, EVAL_ORDERS)
        for density_order, temperature_order in EVAL_ORDERS:
            np.testing.assert_allclose(physical[density_order, temperature_order],
                                       (-1)**density_order * 4**temperature_order * actual[0, 0])

    def test_zero_temperature_pressure_balance(self):
        state = {"P": np.array([1e8]), "rho": np.array([1.]), "T": np.array([1.]),
                 "P_T": np.array([0.]), "E_rho": np.array([1e8]), "E_T": np.array([1e8]),
                 "S_rho": np.array([0.]), "S_T": np.array([1e8])}
        self.assertTrue(np.isnan(identity_residuals(state)["dsp"][0]))
        for residual in identity_residuals(state, balanced=True).values():
            np.testing.assert_array_equal(residual, 0)


if __name__ == "__main__":
    unittest.main()
