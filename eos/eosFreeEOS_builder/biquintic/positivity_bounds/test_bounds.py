"""Self-contained analytic, Hermite, subdivision, domain and core calibrations."""

import sys
from pathlib import Path
import unittest

sys.dont_write_bytecode = True

import numpy as np

from bounds import (
    DTYPE, HERMITE, POWER_TO_BERNSTEIN, adaptive_bound, bernstein_bounds,
    bernstein_value, cell_power, derivative, domain_masks, normalize_closed_cell, polynomial_value,
    sign_powers, subdivide,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from freeeos_potential import LN10, evaluate, quintic_basis, thermodynamics


def polynomial_jets(power, du, dtau):
    jets = np.empty((*power.shape[:-2], 2, 2, 3, 3), dtype=DTYPE)
    for density_end in range(2):
        for temperature_end in range(2):
            for density_order in range(3):
                for temperature_order in range(3):
                    jets[..., density_end, temperature_end, density_order, temperature_order] = polynomial_value(
                        derivative(power, density_order, temperature_order, du, dtau), density_end, temperature_end)
    return jets


class BoundsTests(unittest.TestCase):
    def test_hermite_matches_core_basis(self):
        positions = np.array([0., .13, .5, .87, 1.])
        for order in range(3):
            expected = quintic_basis(positions, .37, order)
            actual = np.zeros_like(expected)
            for endpoint in range(2):
                for nodal_derivative in range(3):
                    coefficients = np.polynomial.polynomial.polyder(HERMITE[3*endpoint+nodal_derivative], order)
                    actual[:, endpoint, nodal_derivative] = np.polynomial.polynomial.polyval(positions, coefficients) * DTYPE(.37)**(nodal_derivative-order)
            np.testing.assert_allclose(actual, expected, atol=2e-13, rtol=2e-13)

    def test_full_biquintic_reproduction_and_nine_jets(self):
        generator = np.random.default_rng(2401)
        expected = generator.normal(size=(7, 6, 6)).astype(DTYPE)
        jets = polynomial_jets(expected, .25, .125)
        actual, magnitude = cell_power(jets, .25, .125)
        np.testing.assert_allclose(actual, expected, atol=3e-13, rtol=3e-13)
        self.assertTrue(np.all(magnitude >= abs(actual)))
        np.testing.assert_allclose(polynomial_jets(actual, .25, .125), jets, atol=1e-10, rtol=1e-12)

    def test_derivative_signs_analytic_natural_coordinates(self):
        du, dtau = DTYPE(.3), DTYPE(.2)
        power = np.zeros((6, 6), dtype=DTYPE)
        power[0, 0], power[1, 0], power[0, 1] = 1, 2*du, 3*dtau
        power[2, 0], power[1, 1], power[0, 2] = 4*du**2, 5*du*dtau, 6*dtau**2
        powers, _ = sign_powers(power, abs(power), du, dtau)
        density_phase, temperature_phase = DTYPE(.37), DTYPE(.63)
        native_u, native_tau = density_phase*du, temperature_phase*dtau
        expected_cv = (3+5*native_u+12*native_tau)-2*(2+8*native_u+5*native_tau)-12+4*5-4*8
        expected_prho = 2+8*native_u+5*native_tau+8
        self.assertAlmostEqual(float(polynomial_value(powers["Cv_numerator"], density_phase, temperature_phase)), float(expected_cv), places=13)
        self.assertAlmostEqual(float(polynomial_value(powers["P_rho"], density_phase, temperature_phase)), float(expected_prho), places=13)

    def test_degree_five_bernstein_elevation(self):
        generator = np.random.default_rng(2402)
        power = generator.normal(size=(6, 6)).astype(DTYPE)
        coefficients = POWER_TO_BERNSTEIN @ power @ POWER_TO_BERNSTEIN.T
        for density_phase, temperature_phase in generator.uniform(size=(20, 2)):
            expected = polynomial_value(power, density_phase, temperature_phase)
            actual = bernstein_value(coefficients, DTYPE(density_phase), DTYPE(temperature_phase))
            self.assertAlmostEqual(float(actual), float(expected), places=13)
            self.assertLessEqual(coefficients.min(), actual)
            self.assertGreaterEqual(coefficients.max(), actual)

    def test_subdivision_maps_exact_subrectangles(self):
        generator = np.random.default_rng(2403)
        coefficients = generator.normal(size=(6, 6)).astype(DTYPE)
        children = subdivide(coefficients)
        for child_index, child in enumerate(children):
            density_side, temperature_side = divmod(child_index, 2)
            for density_phase, temperature_phase in ((.1, .7), (.5, .5), (.91, .13)):
                expected = bernstein_value(coefficients, DTYPE((density_side+density_phase)/2), DTYPE((temperature_side+temperature_phase)/2))
                actual = bernstein_value(child, DTYPE(density_phase), DTYPE(temperature_phase))
                self.assertAlmostEqual(float(actual), float(expected), places=14)

    def test_closed_cell_width_rescaling(self):
        generator = np.random.default_rng(2405)
        power = generator.normal(size=(2, 6, 6)).astype(DTYPE)
        density_stretch = np.array([1.00001, .99999], dtype=DTYPE)
        temperature_stretch = np.array([.99998, 1.00002], dtype=DTYPE)
        normalized = normalize_closed_cell(power, density_stretch, temperature_stretch)
        for density_phase, temperature_phase in ((0., 1.), (.31, .79), (1., 1.)):
            expected = polynomial_value(power, density_stretch*density_phase, temperature_stretch*temperature_phase)
            actual = polynomial_value(normalized, density_phase, temperature_phase)
            np.testing.assert_allclose(actual, expected, atol=1e-16, rtol=1e-16)

    def test_adaptive_positive_not_just_samples(self):
        power = np.zeros((6, 6), dtype=DTYPE)
        power[0, 0], power[1, 0], power[2, 0] = DTYPE(".26"), -1, 1
        coefficients, margin, *_ = bernstein_bounds(power, abs(power))
        self.assertLess(coefficients.min(), 0)
        result = adaptive_bound(coefficients, margin)
        self.assertEqual(result["status"], "positive")
        self.assertGreater(result["lower_bound"], 0)
        self.assertGreater(result["visited_nodes"], 1)
        self.assertEqual(result["unresolved_leaves"], 0)

    def test_zero_minimum_ambiguous_not_failure(self):
        coefficients = np.zeros((6, 6), dtype=DTYPE)
        result = adaptive_bound(coefficients, 1e-20, max_depth=2, max_nodes=21)
        self.assertEqual(result["status"], "unresolved")
        self.assertIsNone(result["witness"])
        self.assertEqual(result["visited_nodes"], 21)
        self.assertEqual(result["unresolved_leaves"], 16)

    def test_positive_minimum_below_margin_unresolved(self):
        result = adaptive_bound(np.full((6, 6), DTYPE("1e-12")), DTYPE("2e-12"), max_depth=0)
        self.assertEqual(result["status"], "unresolved")

    def test_negative_requires_core_witness_callback(self):
        coefficients = -np.ones((6, 6), dtype=DTYPE)
        self.assertEqual(adaptive_bound(coefficients, 1e-12, max_depth=0)["status"], "unresolved")
        result = adaptive_bound(coefficients, 1e-12,
                                witness=lambda density_phase, temperature_phase: {"value": -1., "density_phase": density_phase, "temperature_phase": temperature_phase})
        self.assertEqual(result["status"], "negative_witness")
        self.assertEqual(result["witness_core_queries"], 1)

    def test_node_budget_preserves_unresolved_cover(self):
        result = adaptive_bound(np.zeros((6, 6), dtype=DTYPE), 1e-20, max_depth=4, max_nodes=5)
        self.assertEqual(result["visited_nodes"], 5)
        self.assertEqual(result["unresolved_leaves"], 4)
        self.assertEqual(result["status"], "unresolved")
        with self.assertRaises(ValueError):
            adaptive_bound(np.zeros((6, 6)), 0, max_depth=5)

    def test_full_cell_single_branch_not_corner_union(self):
        logq, logt = np.array([-.51, -.49]), np.array([6.49, 6.51])
        own, material, branches, _ = domain_masks(logq, logt, np.zeros((2, 2)))
        self.assertFalse(own[0, 0])
        self.assertTrue(material[0, 0])
        self.assertFalse(any(branch[0, 0] for branch in branches.values()))

    def test_low_density_maximum_and_material_corner_bounds(self):
        logq, logt = np.array([-2., -1.9]), np.array([3.5, 3.6])
        own, material, _, maximum = domain_masks(logq, logt, np.full((2, 2), 9e-4))
        self.assertTrue(own[0, 0])
        self.assertTrue(material[0, 0])
        ratio = np.full((2, 2), 9e-4)
        ratio[1, 1] = 1e-3
        _, material, _, maximum = domain_masks(logq, logt, ratio)
        self.assertFalse(material[0, 0])
        self.assertEqual(maximum[0, 0], 1e-3)
        _, material, _, _ = domain_masks(np.array([4.99, 5.01]), logt, ratio*0)
        self.assertFalse(material[0, 0])

    def test_core_sample_agreement_with_conservative_margin(self):
        logq, logt = np.array([-3., -2.95]), np.array([4., 4.025])
        du, dtau = np.diff(logq)[0]*LN10, np.diff(logt)[0]*LN10
        generator = np.random.default_rng(2404)
        jets = generator.normal(size=(1, 2, 2, 3, 3)).astype(DTYPE)*1e12
        power, magnitude = cell_power(jets, du, dtau)
        powers, scales = sign_powers(power, magnitude, du, dtau)
        nodes = jets[0].transpose(2, 3, 0, 1).astype(float)
        complete = np.ones((2, 2), dtype=bool)
        for density_phase, temperature_phase in ((.01, .99), (.211, .337), (.5, .5), (.789, .663)):
            sample_q = np.array([logq[0]+np.diff(logq)[0]*density_phase])
            sample_t = np.array([logt[0]+np.diff(logt)[0]*temperature_phase])
            native = evaluate(nodes, complete, logq, logt, sample_q, sample_t)
            state = thermodynamics(native, sample_q, sample_t)
            actual_density = (sample_q-logq[0])/np.diff(logq)[0]
            actual_temperature = (sample_t-logt[0])/np.diff(logt)[0]
            for field in powers:
                _, margin, *_ = bernstein_bounds(powers[field], scales[field])
                value = polynomial_value(powers[field], actual_density, actual_temperature)
                expected = state["Cv"]*state["T"] if field == "Cv_numerator" else state["P_rho"]
                self.assertLess(float(abs(value-expected)[0]), float(margin[0]))


if __name__ == "__main__":
    unittest.main()
