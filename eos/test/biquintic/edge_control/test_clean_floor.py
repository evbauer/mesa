"""Clean source-floor and recursive donor-support checks without runtime probes."""

import unittest

import numpy as np

from clean_floor import build_control, support_report
from edge_control import build_jets
from freeeos_potential import LOW_ORDERS, source_constraints


class CleanFloorTests(unittest.TestCase):
    def source_rows(self):
        logq, logt = np.linspace(-2., -1.4, 21), np.linspace(3.02, 3.42, 21)
        rows = np.zeros((21, 21, 21))
        rows[:, :, 1:4] = [5., 10., 7.]
        rows[:, :, 4:11] = 1.
        rows[:, :, 17] = 0.
        rows[:, :, 18] = logq[:, None] + 2 * logt[None, :] - 12
        return rows, logq, logt

    def test_floor3point1_changes_only_temperature_eligibility(self):
        rows, logq, logt = self.source_rows()
        default, default_valid, _ = source_constraints(rows, logq, logt)
        candidate, valid, _ = source_constraints(rows, logq, logt, minimum_logT=3.1)
        self.assertFalse(valid[:, logt < 3.1 - 1e-12].any())
        self.assertTrue(valid[:, logt >= 3.1 - 1e-12].all())
        for orders in LOW_ORDERS:
            np.testing.assert_array_equal(candidate[orders][default_valid], default[orders][default_valid])
        rows[7, 7, 17] = 1.
        rows[8, 8, 18] += .1
        rows[9, 9, 7] = -1.
        _, masked, _ = source_constraints(rows, logq, logt, minimum_logT=3.1)
        for index in (7, 8, 9):
            self.assertFalse(masked[index, index])

    def test_shifted_recursive_donors_stop_at_clean_floor(self):
        rows, logq, logt = self.source_rows()
        control = build_control(rows, logq, logt, 3.1, False)
        self.assertAlmostEqual(control["support"]["minimum_donor_logT"], 3.1)
        self.assertAlmostEqual(control["support"]["minimum_complete_node_logT"], 3.1)
        self.assertFalse(control["complete"][:, logt < 3.1 - 1e-12].any())
        altered = control["low"].copy()
        altered[:, :, ~control["valid"]] = 1e200
        nodes, complete, metadata = build_jets(altered, control["valid"], .03 * np.log(10), .02 * np.log(10))
        np.testing.assert_array_equal(complete, control["complete"])
        np.testing.assert_allclose(nodes[:, :, complete], control["nodes"][:, :, complete], rtol=1e-11)
        support_report(complete, metadata, logq, logt, 3.1)

    def test_invalid_clean_floors_rejected(self):
        rows, logq, logt = self.source_rows()
        for source_floor in (3.09, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                source_constraints(rows, logq, logt, minimum_logT=source_floor)


if __name__ == "__main__":
    unittest.main()
