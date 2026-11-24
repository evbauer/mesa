"""Coordinate, fraction-label and no-probe analytic checks for the active slice."""

import unittest

import numpy as np

from run_slice import (
    LN10, REGION_NAMES, SLICE_LOGT, baseline, fraction_regions, legacy_state,
    region_spans, slice_queries,
)


class ActiveSliceTests(unittest.TestCase):
    def test_bounded_exact_query_grid_and_faces(self):
        sample_q, sample_t = slice_queries()
        self.assertEqual(len(sample_q), 1001)
        self.assertEqual(sample_q[0], .5)
        self.assertEqual(sample_q[-1], 1.1)
        np.testing.assert_array_equal(sample_t, np.full(1001, SLICE_LOGT))
        for face in (.77, .80):
            self.assertLess(np.min(abs(sample_q - face)), 1e-14)

    def test_actual_fractions_not_guessed_geometric_regions(self):
        values = np.array([1, .9, 0, .2, 1])
        np.testing.assert_array_equal(fraction_regions(values), [0, 1, 2, 1, 0])
        for values in (np.array([np.nan]), np.array([1.1]), np.array([-.1])):
            with self.assertRaises(ValueError):
                fraction_regions(values)

    def test_absent_inactive_region_is_not_invented(self):
        sample_q = np.linspace(.5, 1.1, 5)
        spans = region_spans(sample_q, fraction_regions(np.array([1, 1, .9, .5, .2])))
        self.assertEqual([span["label"] for span in spans], list(REGION_NAMES[:2]))
        self.assertEqual(sum(span["sample_count"] for span in spans), 5)
        self.assertEqual(spans[0]["shading_logQ_bounds_midpoint_estimate"][0], .5)
        self.assertEqual(spans[-1]["shading_logQ_bounds_midpoint_estimate"][-1], 1.1)

    def test_physical_temperature_chain_rule_for_independent_Cv(self):
        logq, logt = np.array([.5, .53]), np.array([6.74, 6.76])
        cache = np.zeros((4, max(baseline.FIELDS.values()) + 1, 2, 2))
        cache[0, baseline.FIELDS["Cv"]] = 10 + 3 * logq[:, None] + 8 * logt[None, :]
        cache[0, baseline.FIELDS["lnPgas"]] = 2
        sample_q, sample_t = np.array([.51, .52]), np.array([6.745, 6.755])
        state = legacy_state(cache, logq, logt, sample_q, sample_t)
        np.testing.assert_allclose(state["Cv"], 10 + 3 * sample_q + 8 * sample_t, rtol=1e-14)
        np.testing.assert_allclose(state["Cv_t"], (8 - 2 * 3) / LN10, rtol=1e-12)
        np.testing.assert_allclose(state["Pgas"], np.exp(2))


if __name__ == "__main__":
    unittest.main()
