"""Small analytic guards for the Figure 2 diagnostic arithmetic."""

import unittest

import numpy as np

from run import compatibility, scaled_jump


class DiagnosticsTests(unittest.TestCase):
    def test_compatibility_quintic_and_sextic(self):
        spacing = 0.37
        positions = np.array([-spacing, 0.0, spacing])
        for degree, expected in ((5, 0.0), (6, 12 * spacing**3)):
            value = positions**degree
            slope = degree * positions**(degree - 1)
            curvature = degree * (degree - 1) * positions**(degree - 2)
            self.assertAlmostEqual(compatibility(value, slope, curvature, spacing),
                                   expected, delta=1e-10)

    def test_jump_scale_has_no_unit_floor(self):
        left = {"Cv": np.array([0.0, 1e-9]), "Cv_x": np.array([0.0, 1e-9])}
        right = {"Cv": np.array([0.0, 2e-9]), "Cv_x": np.array([0.0, 3e-9])}
        jump, scale = scaled_jump(left, right, "Cv_x")
        np.testing.assert_array_equal(jump, [0.0, 2 / 3])
        np.testing.assert_array_equal(scale, [0.0, 3e-9])

    def test_nonfinite_jump_rejected(self):
        left = {"Cv": np.array([np.nan]), "Cv_x": np.array([0.0])}
        right = {"Cv": np.array([1.0]), "Cv_x": np.array([0.0])}
        with self.assertRaisesRegex(ValueError, "Nonfinite"):
            scaled_jump(left, right, "Cv_x")


if __name__ == "__main__":
    unittest.main()
