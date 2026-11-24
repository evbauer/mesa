"""Tests for face-classification and normalized jump reporting."""

import unittest

import numpy as np

from active_faces import classifications, jump_summary


class ActiveFaceTests(unittest.TestCase):
    def test_disjoint_actual_fraction_categories_and_failed_status(self):
        status = np.array([0, 0, 0, 4, 0])
        fraction = np.array([0., .4, 1., .5, np.nan])
        masks = classifications(status, fraction)
        np.testing.assert_array_equal(masks["active"], [False, True, True, False, False])
        np.testing.assert_array_equal(masks["pure_FreeEOS"], [False, False, True, False, False])
        np.testing.assert_array_equal(masks["blend"], [False, True, False, False, False])
        np.testing.assert_array_equal(masks["inactive"], [True, False, False, False, False])
        np.testing.assert_array_equal(masks["failed"], [False, False, False, True, True])

    def test_summary_selects_worst_only_within_category(self):
        result = jump_summary(np.array([.9, -.1, .2]), np.array([False, True, True]),
                              np.array([1., 2., 3.]), np.array([6., 7., 8.]), np.array([0., .5, 1.]))
        self.assertEqual(result["count"], 2)
        self.assertEqual(result["worst"]["logQ"], 3.)
        self.assertEqual(result["worst"]["logRho"], 7.)
        self.assertEqual(result["worst"]["fraction_FreeEOS"], 1.)
        self.assertEqual(jump_summary(np.ones(3), np.zeros(3, dtype=bool), np.ones(3), np.ones(3), np.ones(3)), {"count": 0})


if __name__ == "__main__":
    unittest.main()
