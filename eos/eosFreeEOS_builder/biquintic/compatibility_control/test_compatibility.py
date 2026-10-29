"""Analytic checks for the Appendix A diagnostic and physical shear."""

import unittest

import numpy as np

from compatibility import face_state, line_compatibility, quintic_third_jump
from freeeos_potential import LN10, quintic_basis
from run_control import validate_composition


class CompatibilityTests(unittest.TestCase):
    def test_composition_guard_rejects_swapped_metadata(self):
        header = np.zeros(11)
        header[1:3] = [.70, .02]
        validate_composition(header, {"candidate": {"X": .70, "Z": .02}})
        with self.assertRaisesRegex(ValueError, "manifest composition"):
            validate_composition(header, {"candidate": {"X": .02, "Z": .70}})
        header[1:3] = [.02, .70]
        with self.assertRaisesRegex(ValueError, "source-table composition"):
            validate_composition(header, {"candidate": {"X": .02, "Z": .70}})

    def test_quintic_has_zero_jump_and_sextic_has_expected_jump(self):
        for spacing in (.13, .7, 2.0):
            locations = np.array([-spacing, 0.0, spacing])
            for degree in range(6):
                values = locations**degree
                first = degree * locations**(degree - 1) if degree else np.zeros(3)
                second = degree * (degree - 1) * locations**(degree - 2) if degree > 1 else np.zeros(3)
                self.assertAlmostEqual(quintic_third_jump(values, first, second, spacing), 0.0, delta=2e-10)
            values = locations**6
            first = 6 * locations**5
            second = 30 * locations**4
            self.assertAlmostEqual(quintic_third_jump(values, first, second, spacing), 12 * spacing**3,
                                   delta=2e-10)

    def test_formula_matches_endpoint_hermite_derivatives(self):
        spacing = .27
        values = np.array([1.3, -.2, .9])
        first = np.array([.6, -1.1, .4])
        second = np.array([-.7, .8, 2.1])
        left_basis = quintic_basis(np.array([1.0]), spacing, 3)[0]
        right_basis = quintic_basis(np.array([0.0]), spacing, 3)[0]
        left = sum(left_basis[end, order] * (values[end] if order == 0 else first[end] if order == 1 else second[end])
                   for end in range(2) for order in range(3))
        right = sum(right_basis[end, order] * (values[end + 1] if order == 0 else first[end + 1] if order == 1 else second[end + 1])
                    for end in range(2) for order in range(3))
        self.assertAlmostEqual(right - left, quintic_third_jump(values, first, second, spacing), places=10)

    def test_knot_line_uses_only_supplied_normal_derivatives(self):
        low = np.zeros((3, 3, 3, 2))
        normal = np.array([-1.0, 0.0, 1.0])
        low[0, 0, :, 0] = normal**6
        low[1, 0, :, 0] = 6 * normal**5
        low[2, 0, :, 0] = 30 * normal**4
        self.assertAlmostEqual(line_compatibility(low, "u", 1, 0, 1.0), 12.0)
        low[2, 1, :, 0] = 1e40
        low[1, 2, :, 0] = -1e40
        low[2, 2, :, 0] = 1e40
        self.assertAlmostEqual(line_compatibility(low, "u", 1, 0, 1.0), 12.0)

    def test_physical_cv_jump_shear_on_piecewise_cubic(self):
        logq = np.array([-1.0, 0.0, 1.0])
        logt = np.array([5.0, 6.0])
        nodes = np.zeros((3, 3, 3, 2))
        for density_index, sample_q in enumerate(logq):
            native_q = sample_q * LN10
            for temperature_index, sample_t in enumerate(logt):
                native_t = sample_t * LN10
                nodes[0, 0, density_index, temperature_index] = 1e9 * (native_q - native_t**2)
                nodes[1, 0, density_index, temperature_index] = 1e9
                nodes[0, 1, density_index, temperature_index] = -2e9 * native_t
                nodes[0, 2, density_index, temperature_index] = -2e9
                if density_index == 2:
                    nodes[0, 0, density_index, temperature_index] += 1e7 * native_q**3
                    nodes[1, 0, density_index, temperature_index] += 3e7 * native_q**2
                    nodes[2, 0, density_index, temperature_index] += 6e7 * native_q
        case = {"direction": "u", "logQ": 0.0, "logT": 5.5,
                "left_cell": [0, 0], "right_cell": [1, 0]}
        state = face_state(nodes, np.ones((3, 2), dtype=bool), logq, logt, case, False)
        self.assertAlmostEqual(state["normal_third_jump"], -6e7, delta=1e-3)
        self.assertAlmostEqual(state["Cv_t_jump"], state["shear_expected_Cv_t_jump"], delta=2e-4)
        self.assertAlmostEqual(state["Cv_x_jump"], state["shear_expected_Cv_x_jump"], delta=2e-4)
        self.assertAlmostEqual(state["Cv_t_jump"], -2 * state["Cv_x_jump"], delta=2e-4)

    def test_tau_face_has_cv_t_jump_and_continuous_cv_x(self):
        logq = np.array([-1.0, 0.0])
        logt = np.array([5.0, 6.0, 7.0])
        nodes = np.zeros((3, 3, 2, 3))
        face_t = logt[1] * LN10
        for density_index, sample_q in enumerate(logq):
            native_q = sample_q * LN10
            for temperature_index, sample_t in enumerate(logt):
                native_t = sample_t * LN10
                nodes[0, 0, density_index, temperature_index] = 1e9 * (native_q - native_t**2)
                nodes[1, 0, density_index, temperature_index] = 1e9
                nodes[0, 1, density_index, temperature_index] = -2e9 * native_t
                nodes[0, 2, density_index, temperature_index] = -2e9
                if temperature_index == 2:
                    offset = native_t - face_t
                    nodes[0, 0, density_index, temperature_index] += 1e7 * offset**3
                    nodes[0, 1, density_index, temperature_index] += 3e7 * offset**2
                    nodes[0, 2, density_index, temperature_index] += 6e7 * offset
        case = {"direction": "tau", "logQ": -.5, "logT": 6.0,
                "left_cell": [0, 0], "right_cell": [0, 1]}
        state = face_state(nodes, np.ones((2, 3), dtype=bool), logq, logt, case, False)
        self.assertAlmostEqual(state["normal_third_jump"], -6e7, delta=1e-3)
        self.assertAlmostEqual(state["Cv_t_jump"], state["shear_expected_Cv_t_jump"], delta=2e-4)
        self.assertAlmostEqual(state["Cv_x_jump"], 0.0, delta=2e-4)


if __name__ == "__main__":
    unittest.main()
