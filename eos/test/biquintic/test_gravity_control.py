"""Analytic residual-Jacobian controls with fixed previous-state coordinates."""

import unittest

import numpy as np

from gravity_control import gravity_residual


class GravityControlTests(unittest.TestCase):
    def test_nonconstant_coefficients_chain_rule_in_both_coordinates(self):
        def coefficients(density_coordinate, temperature_coordinate):
            state = {"T": np.exp(temperature_coordinate)}
            for field, scale, density_slope, temperature_slope in (
                    ("Cp", 2e8, .2, .3), ("grad_ad", .3, -.1, .2),
                    ("chiRho", 1.2, .4, -.2), ("chiT", .8, -.3, .1)):
                state[field] = scale * np.exp(density_slope * density_coordinate + temperature_slope * temperature_coordinate)
                state[field + "_x"] = density_slope * state[field]
                state[field + "_t"] = temperature_slope * state[field]
            return state

        density_coordinate, temperature_coordinate = np.array([.1, -.2]), np.array([.4, -.1])
        delta_x, delta_t, step = .03, -.02, 1e-6
        reference = gravity_residual(coefficients(density_coordinate, temperature_coordinate), delta_x, delta_t)
        for suffix in ("_x", "_t"):
            plus_x = density_coordinate + step if suffix == "_x" else density_coordinate
            minus_x = density_coordinate - step if suffix == "_x" else density_coordinate
            plus_t = temperature_coordinate + step if suffix == "_t" else temperature_coordinate
            minus_t = temperature_coordinate - step if suffix == "_t" else temperature_coordinate
            plus = gravity_residual(coefficients(plus_x, plus_t), delta_x + (step if suffix == "_x" else 0), delta_t + (step if suffix == "_t" else 0))
            minus = gravity_residual(coefficients(minus_x, minus_t), delta_x - (step if suffix == "_x" else 0), delta_t - (step if suffix == "_t" else 0))
            np.testing.assert_allclose((plus["value"] - minus["value"]) / (2 * step), reference[suffix], rtol=2e-10)

    def test_ideal_gas_simplification_and_fixed_start_derivatives(self):
        cp, cv, gas_constant = 2.5e8, 1.5e8, 1e8
        temperature = np.array([1e4, 1e6])
        state = {"T": temperature, "Cp": np.full(2, cp), "grad_ad": np.full(2, .4),
                 "chiRho": np.ones(2), "chiT": np.ones(2)}
        for field in ("Cp", "grad_ad", "chiRho", "chiT"):
            state[field + "_x"] = np.zeros(2)
            state[field + "_t"] = np.zeros(2)
        delta_x, delta_t = .01, .03
        result = gravity_residual(state, delta_x, delta_t)
        np.testing.assert_allclose(result["value"], temperature * (gas_constant * delta_x - cv * delta_t))
        np.testing.assert_allclose(result["_x"], temperature * gas_constant)
        np.testing.assert_allclose(result["_t"], result["value"] - temperature * cv)
        step = 1e-6
        plus, minus = state.copy(), state.copy()
        plus["T"], minus["T"] = temperature * np.exp(step), temperature * np.exp(-step)
        finite_difference = (gravity_residual(plus, delta_x, delta_t + step)["value"] -
                             gravity_residual(minus, delta_x, delta_t - step)["value"]) / (2 * step)
        np.testing.assert_allclose(finite_difference, result["_t"], rtol=1e-10)


if __name__ == "__main__":
    unittest.main()
