"""Focused controls for the full-face population and routing mask."""

import unittest

import numpy as np

from audit import classify_probe, face_support, native_faces, support_causes


class PureFaceAuditTests(unittest.TestCase):
    def test_complete_native_face_population_and_midpoints(self):
        logq = np.array([-1., 0., 1., 2.])
        logt = np.array([3.2, 3.3, 3.4, 3.5])
        for axis in (0, 1):
            density_index, temperature_index, face_q, face_t = native_faces(logq, logt, axis)
            self.assertEqual(len(face_q), 6)
            self.assertEqual(len(set(zip(density_index, temperature_index))), 6)
            if axis == 0:
                np.testing.assert_array_equal(face_q, logq[density_index + 1])
                np.testing.assert_allclose(face_t, (logt[temperature_index] + logt[temperature_index + 1]) / 2)
            else:
                np.testing.assert_allclose(face_q, (logq[density_index] + logq[density_index + 1]) / 2)
                np.testing.assert_array_equal(face_t, logt[temperature_index + 1])

    def test_fraction_sum_other_fractions_and_tolerance(self):
        probe = np.zeros((5, 14))
        probe[:, 5] = [1., 1. - 5e-13, 1. - 2e-12, 1., 1.]
        probe[1, 3] = 5e-13
        probe[2, 3] = 2e-12
        probe[3, 3] = 1e-8
        probe[4, 2] = 7
        probe[4, 3] = np.nan
        successful, pure, good = classify_probe(probe)
        np.testing.assert_array_equal(successful, [True, True, True, True, False])
        np.testing.assert_array_equal(pure, [True, True, False, False, False])
        np.testing.assert_array_equal(good[:4], [True, True, True, False])

    def test_two_adjoining_cells_and_separate_support_causes(self):
        logt = np.array([3.05, 3.15, 3.25, 3.35, 3.45])
        complete = np.ones((4, 5), dtype=bool)
        density_index = np.array([0, 0, 0, 0])
        temperature_index = np.array([0, 1, 2, 3])
        self.assertTrue(np.all(face_support(complete, density_index, temperature_index, 0)))
        complete[1, 3] = False
        rows = np.zeros((4, 5, 22))
        source_valid = np.ones((4, 5), dtype=bool)
        rows[1, 2, 17] = 1
        source_valid[1, 2] = False
        supported, cause = support_causes(rows, logt, source_valid, complete,
                                          density_index, temperature_index, 0)
        np.testing.assert_array_equal(supported, [True, True, False, False])
        np.testing.assert_array_equal(cause, [0, 0, 3, 4])
        complete[1, 0] = False
        complete[1, 1] = False
        supported, cause = support_causes(rows, logt, source_valid, complete,
                                          density_index, temperature_index, 0)
        np.testing.assert_array_equal(cause, [1, 2, 3, 4])


if __name__ == "__main__":
    unittest.main()
