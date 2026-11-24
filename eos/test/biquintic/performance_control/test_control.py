"""Small tests of orchestration math; analytic Fortran smoke is --smoke."""

import io
import unittest

import numpy as np

from run_benchmark import analytic_state, parse_report, saved_faces, summarize
from finalize_results import diagnostic_distribution, response_parent_comparison


class ControlTests(unittest.TestCase):
    def test_response_parent_scale_offsets_and_near_zero_derivative(self):
        real_values = np.ones((2, 43))
        ad_values = np.ones((2, 43))
        real_values[:, 10:21] = np.arange(11) + 1e8
        ad_values[:, 10:21] = real_values[:, 10:21]
        real_values[:, 26] = -3.
        ad_values[:, 26] = -3.25
        real_values[:, 37] = 0
        ad_values[:, 37] = .5
        real_values[:, 7] = -2e6
        ad_values[:, 7] = -2.1e6
        scales, differences = response_parent_comparison(real_values, ad_values)
        np.testing.assert_array_equal(scales[:, 21:32], real_values[:, 10:21])
        np.testing.assert_array_equal(scales[:, 32:43], real_values[:, 10:21])
        self.assertAlmostEqual(differences[0, 26], .25 / real_values[0, 15])
        self.assertAlmostEqual(differences[0, 37], .5 / real_values[0, 15])
        self.assertEqual(scales[0, 7], 2.1e6)

    def test_diagnostic_tails_are_not_hidden_by_median(self):
        summary = diagnostic_distribution([1e-10, 1e-9, 2e-8, .08])
        self.assertEqual(summary['counts_above_diagnostic_thresholds']['0.01'], 1)
        self.assertEqual(summary['p50_p90_p99_max'][-1], .08)

    def test_saved_coordinate_conversion_and_repeated_order(self):
        payload = io.BytesIO()
        probe = np.zeros((2, 14))
        probe[:, 5] = 1
        np.savez(payload, u_q=[1., 2.], u_t=[4., 5.], u_probe=probe,
                 tau_q=[3., 4.], tau_t=[6., 7.], tau_probe=probe)
        first, indices, metadata = saved_faces(payload.getvalue(), 3, 12)
        second, second_indices, unused = saved_faces(payload.getvalue(), 3, 12)
        np.testing.assert_array_equal(first, second)
        np.testing.assert_array_equal(indices, second_indices)
        expected_q = np.array([1., 2., 3., 4.])[indices]
        expected_t = np.array([4., 5., 6., 7.])[indices]
        np.testing.assert_allclose(first[:, 0] / np.log(10), expected_q + 2 * expected_t - 12)
        np.testing.assert_allclose(first[:, 1] / np.log(10), expected_t)
        self.assertEqual(metadata['saved_FreeEOS_positive_count'], 3)

    def test_timing_units_and_checksum_pass_normalization(self):
        rows = [{'process': process, 'queries': 1000, 'wall_seconds': .002,
                 'cpu_seconds': .002, 'checksum_per_pass': 5., 'failures': 0}
                for process in [1, 1, 2, 2]]
        summary = summarize(rows)
        self.assertAlmostEqual(summary['median_queries_per_second'], 500000.)
        self.assertAlmostEqual(summary['median_seconds_per_query'], 2e-6)
        self.assertEqual(summary['processes'], 2)
        self.assertEqual(summary['checksum_per_pass_min_max'], [5., 5.])

    def test_report_parser_and_gas_limit(self):
        records = parse_report(' phase_seconds cache_load 1.25\n timing real_eval 1 2 10 .1 .1 4 0\n')
        self.assertEqual(records['phase_seconds'][0], ['cache_load', '1.25'])
        self.assertEqual(records['timing'][0][0], 'real_eval')
        values = analytic_state(0, 0)
        self.assertAlmostEqual(values[0] / 2e7, 1.)
        self.assertAlmostEqual(values[4] / 3e7, 1.)
        self.assertAlmostEqual(values[8], 5 / 3)
        self.assertAlmostEqual(values[10], .4)


if __name__ == '__main__':
    unittest.main()
