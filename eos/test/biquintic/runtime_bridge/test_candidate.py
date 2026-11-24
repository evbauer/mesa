"""Actual component/final API acceptance checks, not stellar acceptance."""

import csv
import json
import unittest

import numpy as np

from api_client import EOSClient, FIELDS, LN10, CRAD, MANIFEST_SHA
from runtime_environment import WORKSPACE, EXPERIMENT, installed_inventory, verify_private_cache, sha256


class CandidateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.output = EXPERIMENT / "candidate_tests"
        cls.before = installed_inventory()
        cls.client = EOSClient(cls.output)
        manifest = WORKSPACE / "freeeos_biquintic_progress/experiments/038_solver_manifest/queries_v2.json"
        if sha256(manifest) != MANIFEST_SHA:
            raise ValueError("Frozen manifest changed")
        cls.cases = json.loads(manifest.read_text())["cases"]
        cls.points = [(case["logRho"], case["logT"]) for case in cls.cases]
        with (WORKSPACE / "freeeos_biquintic_progress/experiments/041_consumer_path_selection/states/classified_paths.csv").open() as stream:
            cls.blends = sorted({(float(row["logRho"]), float(row["logT"]))
                                 for row in csv.DictReader(stream) if row["classification"] == "blend"})[:3]
        cls.metrics = {}

    def successful(self, scenario, route, point):
        result = self.client.query(scenario, route, *point)
        self.assertEqual(result[1], 0, (scenario, route, point))
        self.assertEqual(result[3], 1)
        self.assertTrue(np.all(np.isfinite(result[4:])))
        return result

    def test_component_matches_direct_adapter_and_preserves_hybrid(self):
        for point in self.points + self.blends:
            candidate = self.successful("candidate", "component", point)
            direct = self.successful("direct", "component", point)
            legacy = self.successful("off", "component", point)
            self.assertEqual(candidate[2], 1)
            np.testing.assert_array_equal(candidate[4:], direct[4:])
            np.testing.assert_array_equal(candidate[85:], legacy[85:])
            for block in (4, 31, 58):
                np.testing.assert_array_equal(candidate[block + 3:block + 6], legacy[block + 3:block + 6])

    def test_pure_final_matches_component_and_inactive_is_unchanged(self):
        for point in self.points:
            candidate = self.successful("candidate", "component", point)
            final = self.successful("candidate", "final", point)
            self.assertEqual(final[4 + FIELDS["frac_FreeEOS"]], 1)
            np.testing.assert_array_equal(candidate[4:], final[4:])
        legacy = self.successful("off", "final", (4., 6.))
        inactive = self.successful("candidate", "final", (4., 6.))
        self.assertEqual(inactive[2], 0)
        np.testing.assert_array_equal(inactive[4:], legacy[4:])

    def test_strict_adapter_errors_propagate(self):
        for route in ("component", "final"):
            for scenario, expected in (("composition_fault", 102), ("input_fault", 103)):
                result = self.client.query(scenario, route, *self.points[0])
                self.assertEqual(result[1], expected)
                self.assertEqual(result[2], 1)
        below_floor = (-10., 3.19)
        for route in ("component", "final"):
            result = self.client.query("candidate", route, *below_floor)
            legacy = self.client.query("off", route, *below_floor)
            self.assertEqual(legacy[1], 0)
            self.assertNotEqual(result[1], 0)
            self.assertEqual(result[2], 1)
        self.metrics["below_floor"] = {"logRho": -10., "logT": 3.19, "candidate_status": int(result[1])}

    def test_blend_values_follow_existing_fraction(self):
        fractions = []
        for point in self.blends:
            component = self.successful("candidate", "component", point)
            legacy_component = self.successful("off", "component", point)
            final = self.successful("candidate", "final", point)
            legacy_final = self.successful("off", "final", point)
            fraction = final[4 + FIELDS["frac_FreeEOS"]]
            self.assertGreater(fraction, 0)
            self.assertLess(fraction, 1)
            fractions.append(float(fraction))
            expected = legacy_final[4:31] + fraction * (component[4:31] - legacy_component[4:31])
            np.testing.assert_allclose(final[4:31], expected, rtol=5e-14, atol=1e-14)
        self.metrics["blend_fractions"] = fractions

    def test_off_node_identities_and_physical_partials(self):
        diagnostics, derivative_records = [], []
        for case in self.cases:
            face = np.array([case["logRho"], case["logT"]]) * LN10
            for side in (-1, 1):
                point = face.copy()
                point[0 if case["direction"] == "u" else 1] += side * .001
                center = self.successful("candidate", "component", point / LN10)
                values, density, thermal = center[4:31], center[31:58], center[58:85]
                rho, temperature = np.exp(point)
                gas, energy, entropy = np.exp(values[:3])
                radiation = CRAD * temperature**4 / 3
                pressure = gas + radiation
                pressure_t = gas * thermal[0] + 4 * radiation
                balances = [
                    (rho * energy * density[1] + pressure_t - pressure) / pressure,
                    temperature * entropy * thermal[2] / (energy * thermal[1]) - 1,
                    -rho * temperature * entropy * density[2] / pressure_t - 1,
                    energy * thermal[1] / (temperature * values[FIELDS["Cv"]]) - 1]
                diagnostics.append(balances)
                for route in ("component", "final"):
                    central = self.successful("candidate", route, point / LN10)
                    for axis, start in ((0, 31), (1, 58)):
                        analytic = central[start:start + 27]
                        analytic = analytic[:16]
                        errors = []
                        for step in (2e-4, 1e-4, 3e-5, 1e-5):
                            offset = np.zeros(2)
                            offset[axis] = step
                            plus = self.successful("candidate", route, (point + offset) / LN10)
                            minus = self.successful("candidate", route, (point - offset) / LN10)
                            outer_plus = self.successful("candidate", route, (point + 2 * offset) / LN10)
                            outer_minus = self.successful("candidate", route, (point - 2 * offset) / LN10)
                            finite = (8 * (plus[4:20] - minus[4:20]) - (outer_plus[4:20] - outer_minus[4:20])) / (12 * step)
                            scale = np.maximum.reduce([np.abs(analytic), np.abs(central[4:20]), np.ones(16)])
                            errors.append(np.abs(finite - analytic) / scale)
                        best = np.min(errors, axis=0)
                        derivative_records.append({"case_id": case["case_id"], "side": side, "route": route,
                            "axis": axis, "best_scaled_error": best.tolist(), "errors_by_step": np.asarray(errors).tolist()})
        diagnostics = np.abs(diagnostics)
        self.assertLess(np.max(diagnostics), 2e-11)
        self.metrics["identity_maxima"] = np.max(diagnostics, axis=0).tolist()
        self.metrics["identity_queries"] = len(diagnostics)
        self.metrics["finite_difference"] = derivative_records
        regular_columns = [column for column in range(16) if column != FIELDS["dE_dRho"]]
        self.assertLess(max(max(record["best_scaled_error"][column] for column in regular_columns)
                            for record in derivative_records), 3e-6)
        self.assertLess(max(record["best_scaled_error"][FIELDS["dE_dRho"]] for record in derivative_records), 2e-5)

    @classmethod
    def tearDownClass(cls):
        cls.client.close()
        verify_private_cache(cls.output)
        after = installed_inventory()
        if cls.before != after:
            raise AssertionError("Installed cache or EOS library changed")
        cls.metrics["installed_unchanged"] = True
        cls.metrics["scope"] = "actual component/final API; tagged legacy abundance fields; no stellar evolution"
        (cls.output / "metrics.json").write_text(json.dumps(cls.metrics, indent=2) + "\n")


if __name__ == "__main__":
    unittest.main()
