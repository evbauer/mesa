"""Independent polynomial, estimator, continuity, and donor-isolation checks."""

import math
import sys
import unittest
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

from c3_potential import (
    COEFFICIENTS, EXACT_COEFFICIENTS, ESTIMATED_ORDERS, JET_ORDERS, ROUTE_GRAPH,
    construct_jets, evaluate, septic_basis,
)
from freeeos_potential import EVAL_ORDERS, LN10, LOW_ORDERS, source_constraints, thermodynamic_partials


def polynomial_jets(coefficients, logq, logt):
    first, second = np.meshgrid(logq * LN10, logt * LN10, indexing="ij")
    jets = np.empty((4, 4, len(logq), len(logt)))
    for orders in JET_ORDERS:
        differentiated = np.polynomial.polynomial.polyder(coefficients, orders[0], axis=0)
        differentiated = np.polynomial.polynomial.polyder(differentiated, orders[1], axis=1)
        jets[orders] = np.polynomial.polynomial.polyval2d(first, second, differentiated)
    return jets


def synthetic_rows(logq, logt):
    temperature = np.broadcast_to(10.0**logt, (len(logq), len(logt)))
    density = 10.0**(logq[:, None] + 2 * logt - 12)
    rows = np.zeros((*density.shape, 22))
    rows[..., 0] = logt
    rows[..., 1] = np.log10(density * temperature * 2e8)
    rows[..., 2] = np.log10(temperature * 3e8)
    rows[..., 3] = np.log10(np.full_like(density, 1e10))
    rows[..., 4:6] = 1
    rows[..., 6:8] = 3e8
    rows[..., 18] = np.log10(density)
    return rows


class C3Tests(unittest.TestCase):
    def test_exact_polynomial_matrix_independent_cardinality(self):
        matrix = np.zeros((8, 8))
        for endpoint in range(2):
            for derivative in range(4):
                for power in range(derivative, 8):
                    matrix[4 * endpoint + derivative, power] = (
                        math.prod(range(power - derivative + 1, power + 1))
                        * endpoint**(power - derivative))
        np.testing.assert_allclose(matrix @ COEFFICIENTS.T, np.eye(8), atol=2e-13)
        for constraint, coefficients in enumerate(EXACT_COEFFICIENTS):
            for endpoint in range(2):
                for derivative in range(4):
                    actual = sum(coefficients[power] * math.prod(range(power - derivative + 1, power + 1))
                                 * endpoint**(power - derivative) for power in range(derivative, 8))
                    self.assertEqual(actual, int(constraint == 4 * endpoint + derivative))
        for spacing in (.13, 1.7):
            for derivative in range(4):
                basis = septic_basis(np.array([0., 1.]), spacing, derivative)
                expected = np.zeros((2, 2, 4))
                expected[0, 0, derivative] = expected[1, 1, derivative] = 1
                np.testing.assert_array_equal(basis, expected)

    def test_tensor_degree7_all16_partials_and_nodes(self):
        generator = np.random.default_rng(8011)
        coefficients = generator.normal(size=(8, 8))
        logq, logt = np.array([-.21, .18]), np.array([-.12, .22])
        jets = polynomial_jets(coefficients, logq, logt)
        sample_q = np.r_[generator.uniform(*logq, 100), logq[[0, 0, 1, 1]]]
        sample_t = np.r_[generator.uniform(*logt, 100), logt[[0, 1, 0, 1]]]
        actual = evaluate(jets, np.ones((2, 2), bool), logq, logt, sample_q, sample_t, JET_ORDERS)
        for orders in JET_ORDERS:
            differentiated = np.polynomial.polynomial.polyder(coefficients, orders[0], axis=0)
            differentiated = np.polynomial.polynomial.polyder(differentiated, orders[1], axis=1)
            expected = np.polynomial.polynomial.polyval2d(sample_q * LN10, sample_t * LN10, differentiated)
            np.testing.assert_allclose(actual[orders], expected, atol=3e-8, rtol=2e-8)
            np.testing.assert_array_equal(actual[orders][-4:], jets[orders].ravel())

    def test_full_estimator_polynomial_and_explicit_routes(self):
        generator = np.random.default_rng(19)
        logq, logt = np.linspace(-.3, .3, 45), np.linspace(-.25, .25, 43)
        for radius, degree in ((2, 4), (3, 6)):
            coefficients = generator.normal(size=(degree + 1, degree + 1))
            exact = polynomial_jets(coefficients, logq, logt)
            low = exact[:3, :3].copy()
            for orders in ((2, 1), (1, 2), (2, 2)):
                low[orders] = 1e190
            nodes, complete, diagnostics = construct_jets(low, np.ones(low.shape[2:], bool),
                (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10, radius)
            self.assertTrue(complete.any())
            for orders in JET_ORDERS:
                np.testing.assert_allclose(nodes[orders][complete], exact[orders][complete], atol=4e-7, rtol=2e-7)
            for orders in ESTIMATED_ORDERS:
                route = diagnostics["routes"][orders]
                self.assertEqual(route["parents"], ROUTE_GRAPH[orders])
                expected = np.mean(route["values"], axis=0)
                np.testing.assert_array_equal(nodes[orders][complete], expected[complete])
            self.assertEqual(len(diagnostics["routes"][3, 0]["values"]), 1)
            self.assertEqual(len(diagnostics["routes"][0, 3]["values"]), 1)
            self.assertEqual(len(diagnostics["routes"][3, 3]["values"]), 2)

    def test_every_stage_propagates_invalid_and_zero_weight_center(self):
        generator = np.random.default_rng(8)
        low = generator.normal(size=(3, 3, 65, 61))
        valid = np.ones((65, 61), bool)
        valid[31, 30] = False
        valid[:, :5] = False
        for radius in (2, 3):
            reference, complete, diagnostics = construct_jets(low, valid, .12, .09, radius)
            for poison in (1e200, np.nan, np.inf):
                altered_low = low.copy()
                altered_low[:, :, ~valid] = poison
                altered, accepted, altered_diagnostics = construct_jets(altered_low, valid, .12, .09, radius)
                np.testing.assert_array_equal(accepted, complete)
                np.testing.assert_array_equal(altered[:, :, complete], reference[:, :, complete])
                for orders in JET_ORDERS:
                    mask = diagnostics["jet_masks"][orders]
                    np.testing.assert_array_equal(mask, altered_diagnostics["jet_masks"][orders])
                    self.assertFalse(mask[31, 30])
                self.assertFalse(complete[:, :5 + 3 * radius].any())
            for orders in ESTIMATED_ORDERS:
                expected_mask = np.ones(valid.shape, bool)
                for parent, axis in ROUTE_GRAPH[orders]:
                    parent_mask = diagnostics["jet_masks"][parent]
                    donor_mask = np.zeros_like(valid)
                    for first in range(valid.shape[0]):
                        for second in range(valid.shape[1]):
                            center = (first, second)[axis]
                            if radius <= center < valid.shape[axis] - radius:
                                donor_mask[first, second] = all(parent_mask[
                                    first + offset if axis == 0 else first,
                                    second + offset if axis == 1 else second]
                                    for offset in range(-radius, radius + 1))
                    expected_mask &= donor_mask
                np.testing.assert_array_equal(expected_mask, diagnostics["jet_masks"][orders])

    def test_source_fallback_cold_iq520_and_nonfinite_poisoning(self):
        logq, logt = np.linspace(-4, -1, 65), np.linspace(3, 4.5, 65)
        rows = synthetic_rows(logq, logt)
        rows[30, 30, 17] = 520
        rows[35, 35, 17] = 1
        rows[40, 40, 9] = np.nan
        rows[42, 42, 7] = -1
        low, valid, _ = source_constraints(rows, logq, logt)
        self.assertFalse(valid[30, 30] or valid[35, 35] or valid[40, 40] or valid[42, 42])
        self.assertFalse(valid[:, logt < 3.2 - 1e-12].any())
        for radius in (2, 3):
            baseline, complete, _ = construct_jets(low, valid, .1, .08, radius)
            poisoned_rows = rows.copy()
            poisoned_rows[..., 1:11][~valid] = 150
            poisoned_rows[..., 17][~valid] = 520
            poisoned, poisoned_valid, _ = source_constraints(poisoned_rows, logq, logt)
            np.testing.assert_array_equal(poisoned_valid, valid)
            altered, accepted, _ = construct_jets(poisoned, poisoned_valid, .1, .08, radius)
            np.testing.assert_array_equal(accepted, complete)
            np.testing.assert_array_equal(altered[:, :, complete], baseline[:, :, complete])

    def test_exact_face_all16_and_response_jacobians(self):
        generator = np.random.default_rng(34)
        logq, logt = np.linspace(-.4, .4, 4), np.linspace(3.3, 3.7, 4)
        jets = generator.normal(size=(4, 4, 4, 4)) + 3
        complete = np.ones((4, 4), bool)
        for axis in (0, 1):
            for fraction in (0., .37, 1.):
                sample_q = np.array([logq[2] if axis == 0 else logq[1] + fraction * (logq[2] - logq[1])])
                sample_t = np.array([logt[1] + fraction * (logt[2] - logt[1]) if axis == 0 else logt[2]])
                left = (np.array([1]), np.array([1]))
                right = (np.array([2 if axis == 0 else 1]), np.array([1 if axis == 0 else 2]))
                first = evaluate(jets, complete, logq, logt, sample_q, sample_t, JET_ORDERS, left)
                second = evaluate(jets, complete, logq, logt, sample_q, sample_t, JET_ORDERS, right)
                np.testing.assert_array_equal(first, second)
                first_state = thermodynamic_partials(first, sample_q, sample_t)
                second_state = thermodynamic_partials(second, sample_q, sample_t)
                for field in first_state:
                    np.testing.assert_array_equal(first_state[field], second_state[field])

    def test_bounds_invalid_cells_and_unsupplied_purethird(self):
        logq, logt = np.arange(3.), np.arange(3.)
        jets = np.ones((4, 4, 3, 3))
        complete = np.ones((3, 3), bool)
        for sample in (-1., 3., np.nan):
            with self.assertRaises(ValueError):
                evaluate(jets, complete, logq, logt, np.array([sample]), np.array([.5]))
        complete[0, 0] = False
        with self.assertRaises(ValueError):
            evaluate(jets, complete, logq, logt, [.5], [.5])
        for radius in (1, 4):
            with self.assertRaises(ValueError):
                construct_jets(jets[:3, :3], complete, 1, 1, radius)
        self.assertEqual(len(EVAL_ORDERS), 10)
        self.assertEqual(len(JET_ORDERS), 16)
        self.assertNotIn((3, 0), LOW_ORDERS)
        self.assertNotIn((0, 3), LOW_ORDERS)

    def test_redundant_route_coordinate_exchange_symmetry(self):
        generator = np.random.default_rng(67)
        low = generator.normal(size=(3, 3, 49, 51))
        valid = np.ones((49, 51), bool)
        valid[22, 25] = False
        for radius in (2, 3):
            nodes, complete, _ = construct_jets(low, valid, .17, .09, radius)
            exchanged, exchanged_complete, _ = construct_jets(low.transpose(1, 0, 3, 2),
                valid.T, .09, .17, radius)
            np.testing.assert_array_equal(complete, exchanged_complete.T)
            np.testing.assert_array_equal(nodes, exchanged.transpose(1, 0, 3, 2))

    def test_exponential_full_estimator_and_interpolant_convergence(self):
        from run_control import analytic_convergence

        result = analytic_convergence()
        coarse, _, fine = result["grids"]
        for radius, minimum_improvement in ((2, 10), (3, 30)):
            name = f"estimated_c3_r{radius}"
            for orders in ESTIMATED_ORDERS:
                before = coarse["models"][name]["sixteen_estimated_nodal_jets"][str(orders)]["max"]
                after = fine["models"][name]["sixteen_estimated_nodal_jets"][str(orders)]["max"]
                self.assertGreater(before / after, minimum_improvement)
            before = max(coarse["models"][name]["ten_interpolated_partials"][str(orders)]["max"] for orders in EVAL_ORDERS)
            after = max(fine["models"][name]["ten_interpolated_partials"][str(orders)]["max"] for orders in EVAL_ORDERS)
            self.assertGreater(before / after, minimum_improvement)

    def test_response_jacobians_against_physical_coordinate_fd(self):
        logq, logt = np.linspace(-.3, .3, 7), np.linspace(-.2, .2, 7)
        first, second = np.meshgrid(logq * LN10, logt * LN10, indexing="ij")
        jets = np.empty((4, 4, 7, 7))
        for orders in JET_ORDERS:
            jets[orders] = 1.3**orders[0] * .9**orders[1] * np.exp(1.3 * first + .9 * second)
        complete = np.ones((7, 7), bool)
        sample_q, sample_t = np.array([.033, -.071]), np.array([.021, -.095])

        def state(query_q, query_t):
            return thermodynamic_partials(evaluate(jets, complete, logq, logt, query_q, query_t), query_q, query_t)

        center = state(sample_q, sample_t)
        for step in (2e-5, 1e-5):
            for suffix, direction in (("_x", (1, 0)), ("_t", (-2, 1))):
                delta_q, delta_t = np.asarray(direction) * step / LN10
                plus, minus = state(sample_q + delta_q, sample_t + delta_t), state(sample_q - delta_q, sample_t - delta_t)
                for field in ("P", "E", "S", "Cv", "Cp", "chiRho", "chiT", "gamma1", "gamma3", "grad_ad"):
                    finite_difference = (plus[field] - minus[field]) / (2 * step)
                    scale = np.maximum.reduce([abs(center[field]), abs(center[field + suffix]), np.ones(len(sample_q))])
                    np.testing.assert_allclose((finite_difference - center[field + suffix]) / scale, 0, atol=2e-7)

    def test_stage_overflow_is_invalid_not_a_complete_jet(self):
        generator = np.random.default_rng(81)
        low = generator.uniform(-1, 1, size=(3, 3, 41, 41)) * 1e308
        with np.errstate(over="ignore", invalid="ignore"):
            nodes, complete, diagnostics = construct_jets(low, np.ones((41, 41), bool), 1e-100, 1e-100, 2)
        self.assertFalse(complete.any())
        self.assertTrue(np.isnan(nodes).all())
        for orders in ESTIMATED_ORDERS:
            self.assertFalse(diagnostics["jet_masks"][orders].any())

    def test_coarse_only_estimator_ignores_all_withheld_rows(self):
        from run_control import build_models

        logq, logt = np.linspace(-4, -1, 81), np.linspace(3.2, 4.5, 81)
        rows = synthetic_rows(logq, logt)
        baseline, _, common = build_models(logq[::2], logt[::2], rows[::2, ::2])
        self.assertTrue(common.any())
        altered = rows.copy()
        withheld = np.ones(rows.shape[:2], bool)
        withheld[::2, ::2] = False
        altered[withheld] = np.inf
        poisoned, _, poisoned_common = build_models(logq[::2], logt[::2], altered[::2, ::2])
        np.testing.assert_array_equal(common, poisoned_common)
        for name in baseline:
            np.testing.assert_array_equal(baseline[name]["complete"], poisoned[name]["complete"])
            np.testing.assert_array_equal(baseline[name]["nodes"], poisoned[name]["nodes"])


if __name__ == "__main__":
    unittest.main()
