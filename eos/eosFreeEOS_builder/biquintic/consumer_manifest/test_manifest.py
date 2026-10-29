"""Deterministic sampling and native-face geometry controls."""

import unittest

import numpy as np

from build_manifest import composition_metadata, face_indices, select_indices
from probe_paths import classify


class ManifestTests(unittest.TestCase):
    def test_source_header_composition_order(self):
        self.assertEqual(composition_metadata(np.array([51., .7, .02])), {"X": .7, "Z": .02})

    def test_balanced_distinct_deterministic_selection(self):
        jumps = np.arange(100, dtype=float)
        eligible = np.ones(100, dtype=bool)
        eligible[99] = False
        selected = select_indices(jumps, eligible)
        np.testing.assert_array_equal(selected["tail"], [98, 97, 96, 95, 94, 93])
        np.testing.assert_array_equal(selected["ordinary"], [20, 29, 39, 49, 59, 69])
        self.assertEqual(len(np.unique(np.concatenate(tuple(selected.values())))), 12)

    def test_insufficient_population_rejected(self):
        with self.assertRaises(ValueError):
            select_indices(np.arange(10), np.ones(10, dtype=bool))

    def test_native_faces(self):
        knots = np.arange(5, dtype=float)
        self.assertEqual(face_indices(knots, knots, 2., 1.4, "u"), ([1, 1], [2, 1]))
        self.assertEqual(face_indices(knots, knots, 1.4, 2., "tau"), ([1, 1], [1, 2]))
        with self.assertRaises(ValueError):
            face_indices(knots, knots, 1.4, 1.4, "u")

    def test_path_classification(self):
        self.assertEqual(classify(0, 1.), "pure_FreeEOS")
        self.assertEqual(classify(0, .5), "blend")
        self.assertEqual(classify(0, 0.), "inactive")
        self.assertEqual(classify(-1, 1.), "failed")
        self.assertEqual(classify(0, np.nan), "failed")
        self.assertEqual(classify(0, 2.), "failed")


if __name__ == "__main__":
    unittest.main()
