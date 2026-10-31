"""Independent analytic control for the legacy logarithmic chain rules."""

import unittest

import numpy as np

from run_extended import LN10, baseline, legacy_state


class LegacyChainRuleTests(unittest.TestCase):
    def test_energy_entropy_and_both_cv_partials(self):
        logq, logt = np.array([.5, .53]), np.array([6.74, 6.76])
        cache = np.zeros((4, max(baseline.FIELDS.values()) + 1, 2, 2))
        cache[0, baseline.FIELDS["Cv"]] = 10 + 3 * logq[:, None] + 8 * logt[None, :]
        cache[0, baseline.FIELDS["lnS"]] = -2 + .4 * logq[:, None] + 1.1 * logt[None, :]
        cache[0, baseline.FIELDS["lnE"]] = 10 + .3 * logq[:, None] + .8 * logt[None, :]
        cache[0, baseline.FIELDS["dS_dT"]] = 9e99
        sample_q, sample_t = np.array([.51, .52]), np.array([6.745, 6.755])
        state = legacy_state(cache, logq, logt, sample_q, sample_t)
        cv = 10 + 3 * sample_q + 8 * sample_t
        entropy = np.exp(-2 + .4 * sample_q + 1.1 * sample_t)
        energy = np.exp(10 + .3 * sample_q + .8 * sample_t)
        expected_entropy_partial = entropy * (1.1 - 2 * .4) / LN10
        expected_energy_partial = energy / 10.0**sample_t * (.8 - 2 * .3) / LN10
        np.testing.assert_allclose(state["Cv_x"], 3 / LN10, rtol=1e-11)
        np.testing.assert_allclose(state["Cv_t"], (8 - 2 * 3) / LN10, rtol=1e-11)
        np.testing.assert_allclose(state["T_s_T"], expected_entropy_partial, rtol=1e-11)
        np.testing.assert_allclose(state["E_T"], expected_energy_partial, rtol=1e-11)
        np.testing.assert_allclose(state["dse_Cv"], expected_entropy_partial / cv - 1, rtol=1e-11)
        np.testing.assert_allclose(state["eT_Cv"], expected_energy_partial / cv - 1, rtol=1e-11)
        self.assertTrue(np.all(abs(state["dse_Cv"]) > .1))


if __name__ == "__main__":
    unittest.main()
