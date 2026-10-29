"""Signed scaling and classification tests without running a MESA probe."""

import json
import unittest

import numpy as np

from active_quality import JACOBIANS, TABLE_FIELDS, VALUES, classified_groups, signed_comparison


class ActiveQualityTests(unittest.TestCase):
    def test_material_is_not_an_own_use_or_fraction_mask(self):
        coordinates = np.array([-1., .1, 4.52, 4.52])
        state = {field: np.ones(4) for field in set(VALUES + JACOBIANS + ("P_rho",))}
        reference = {field: np.ones(4) for field in TABLE_FIELDS}
        population = {"q": coordinates, "t": np.full(4, 3.56), "overlap": np.zeros(4, bool),
            "state": state, "reference": reference}
        probe = np.zeros((4, 14))
        probe[:, 5] = [1., .25, 0., np.nan]
        probe[-1, 2] = -1004
        result = classified_groups(population, probe)
        groups = result["new_coverage"]["categories"]
        self.assertEqual(groups["all_successful"]["count"], 3)
        self.assertEqual(groups["active"]["count"], 2)
        self.assertEqual(groups["pure_FreeEOS"]["count"], 1)
        self.assertEqual(groups["blend"]["count"], 1)
        self.assertEqual(groups["own_use"]["count"], 1)
        self.assertEqual(groups["active_own_use"]["count"], 1)
        self.assertEqual(groups["inactive"]["count"], 1)
        self.assertEqual(groups["failed"]["count"], 1)
        self.assertIsNone(groups["failed"]["metrics"]["Cv"]["worst"]["fraction_FreeEOS"])
        json.dumps(result, allow_nan=False)

    def test_holdout_signed_departure_retains_underestimate(self):
        candidate = {field: np.array([.2, 1.2]) for field in TABLE_FIELDS}
        reference = {field: np.ones(2) for field in TABLE_FIELDS}
        result = signed_comparison(candidate, reference, np.array([4.52, -1.]), np.full(2, 3.56), np.array([0., 1.]))
        self.assertAlmostEqual(result["Cv"]["max"], .8)
        self.assertAlmostEqual(result["Cv"]["worst"]["signed_scaled_difference"], -.8)
        self.assertEqual(result["Cv"]["worst"]["fraction_FreeEOS"], 0.)
        self.assertAlmostEqual(result["Cv"]["signed"]["median"], -.3)

    def test_face_jacobian_scales_include_response_not_unit_floor(self):
        left = {field: np.array([10.]) for field in VALUES + JACOBIANS}
        right = {field: np.array([10.]) for field in VALUES + JACOBIANS}
        left["Cv_t"], right["Cv_t"] = np.array([1.]), np.array([-2.])
        result = signed_comparison(left, right, np.array([-1.]), np.array([3.22]), np.array([1.]), faces=True)
        self.assertAlmostEqual(result["Cv_t"]["max"], .3)
        self.assertAlmostEqual(result["Cv_t"]["worst"]["signed_scaled_difference"], .3)
        left["chiT"], right["chiT"] = np.array([.01]), np.array([.01])
        left["chiT_t"], right["chiT_t"] = np.array([.002]), np.array([-.003])
        result = signed_comparison(left, right, np.array([-1.]), np.array([3.22]), np.array([1.]), faces=True)
        self.assertAlmostEqual(result["chiT_t"]["max"], .5)


if __name__ == "__main__":
    unittest.main()
