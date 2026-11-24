"""Guard common-forcing pairing against incomplete or mismatched controls."""

import copy
import unittest

from review_coupled import common_pairs


class CoupledReviewTests(unittest.TestCase):
    def setUp(self):
        record = {"family": "common", "case_id": "face", "sign": 1, "scale": .001,
                  "density_over_temperature_increment": -1, "initial": [1., 2.],
                  "nominal_target": [1.01, 1.99], "pressure_target": 3.,
                  "heating": 4., "energy_scale": 5., "box": [[0., 2.], [1., 3.]]}
        self.records = [{**copy.deepcopy(record), "model": model} for model in ("candidate", "legacy")]

    def test_matched_pair(self):
        self.assertEqual(len(common_pairs(self.records, expected_count=1)), 1)

    def test_incomplete_and_duplicate_models(self):
        for records in (self.records[:1], self.records + self.records[:1]):
            with self.assertRaises(ValueError):
                common_pairs(records, expected_count=1)

    def test_different_forcing_is_rejected(self):
        self.records[1]["heating"] += .01
        with self.assertRaisesRegex(ValueError, "heating"):
            common_pairs(self.records, expected_count=1)


if __name__ == "__main__":
    unittest.main()
