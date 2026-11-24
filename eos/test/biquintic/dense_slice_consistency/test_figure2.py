"""Synthetic controls for the returned entropy derivative."""

import unittest

import numpy as np

from run_figure2 import LN10, legacy, legacy_state


class Figure2ConsistencyTests(unittest.TestCase):
    def test_legacy_log_entropy_chain_rule_and_fixed_density_shear(self):
        logq = np.array([.5, .53])
        logt = np.array([6.74, 6.76])
        cache = np.zeros((4, max(legacy.FIELDS.values()) + 1, 2, 2))
        cache[0, legacy.FIELDS["Cv"]] = 7 + 2 * logq[:, None] + 3 * logt[None, :]
        cache[0, legacy.FIELDS["lnS"]] = -3 + .4 * logq[:, None] + 1.1 * logt[None, :]
        cache[0, legacy.FIELDS["dS_dT"]] = 1e80
        sample_q = np.array([.508, .518])
        sample_t = np.array([6.745, 6.755])
        state = legacy_state(cache, logq, logt, sample_q, sample_t)
        expected_entropy = np.exp(-3 + .4 * sample_q + 1.1 * sample_t)
        expected_partial = expected_entropy * (1.1 - .8) / LN10
        np.testing.assert_allclose(state["S"], expected_entropy, rtol=1e-12)
        np.testing.assert_allclose(state["T_s_T"], expected_partial, rtol=1e-11)
        np.testing.assert_allclose(state["dse_Cv"], expected_partial / state["Cv"] - 1, rtol=1e-11)
        step = 1e-5
        plus = legacy_state(cache, logq, logt, sample_q - 2 * step / LN10, sample_t + step / LN10)["S"]
        minus = legacy_state(cache, logq, logt, sample_q + 2 * step / LN10, sample_t - step / LN10)["S"]
        np.testing.assert_allclose((plus - minus) / (2 * step), expected_partial, rtol=2e-9)


if __name__ == "__main__":
    unittest.main()
