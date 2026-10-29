"""Analytic controls for the standalone thermal residual and solver."""

import sys
import unittest
import json
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from solver import reference, residual, solve
from run_control import MANIFEST, verify_manifest


class SolverTests(unittest.TestCase):
    def test_constant_cv_jacobian_root_and_temporal_defect(self):
        cv = 1.5e8
        model = lambda temperature_log: (cv, 0.)
        initial = np.log(1e6)
        errors = []
        for step in (1e-2, 5e-3, 2.5e-3):
            heating = step * np.exp(initial) * cv
            lower, upper = initial, initial + 2 * step
            root, diagnostic = reference(model, initial, heating, lower, upper,
                                         np.linspace(lower, upper, 21))
            trial = solve(model, initial, heating, lower, upper, initial)
            self.assertEqual(trial.status, "converged")
            self.assertLess(abs(trial.root - root), 1e-10)
            self.assertTrue(diagnostic["sampled_positive_jacobian"])
            errors.append(abs(root - np.log(np.exp(initial) + heating / cv)))
        self.assertAlmostEqual(errors[0] / errors[1], 4, delta=.1)
        self.assertAlmostEqual(errors[1] / errors[2], 4, delta=.1)

    def test_residual_derivative_and_both_heating_signs(self):
        model = lambda temperature_log: (2e8 * np.exp(.2 * temperature_log),
                                         .2 * 2e8 * np.exp(.2 * temperature_log))
        initial = np.log(1e6)
        for heating in (-1e12, 1e12):
            value, jacobian = residual(model, initial, initial, heating)
            step = 1e-5
            finite = (residual(model, initial + step, initial, heating)[0]
                      - residual(model, initial - step, initial, heating)[0]) / (2 * step)
            self.assertAlmostEqual(jacobian, finite, delta=1e-9)
            self.assertEqual(np.sign(value), -np.sign(heating))

    def test_piecewise_cv_face_one_sided_and_invalid_cv(self):
        face = np.log(1e6)

        def model(temperature_log):
            offset = temperature_log - face
            slope = 1e8 if offset < 0 else 4e8
            return 2e8 + slope * offset, slope

        heating = 1e12
        left = residual(model, face - 1e-8, face - .002, heating)[1]
        right = residual(model, face + 1e-8, face - .002, heating)[1]
        self.assertGreater(abs(left - right), 1e-4)
        with self.assertRaises(ValueError):
            residual(lambda temperature_log: (0., 0.), face, face, heating)

    def test_manifest_rejects_swapped_composition_and_old_schema(self):
        manifest = json.loads(MANIFEST.read_text())
        manifest["candidate"]["X"], manifest["candidate"]["Z"] = .02, .7
        with self.assertRaisesRegex(ValueError, "X/Z metadata"):
            verify_manifest(None, manifest)
        manifest["candidate"]["X"], manifest["candidate"]["Z"] = .7, .02
        manifest["schema"] = "FREEEOS_CONSUMER_MANIFEST_V1"
        with self.assertRaisesRegex(ValueError, "schema/count"):
            verify_manifest(None, manifest)


if __name__ == "__main__":
    unittest.main()
