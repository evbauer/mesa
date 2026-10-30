"""Actual-MESA opt-out, hook transport, blend and strict-error regression."""

import csv
import io
import json
import subprocess
import unittest
from pathlib import Path

import numpy as np

from runtime_environment import (CODE, WORKSPACE, MESA, EXPERIMENT, MASSFILE, INLIST,
                                 installed_inventory, private_environment, sha256, verify_private_cache)


class HookTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.output = EXPERIMENT / "hook_tests"
        cls.output.mkdir(parents=True, exist_ok=True)
        cls.before = installed_inventory()
        cls.environment = private_environment(cls.output)
        manifest = WORKSPACE / "freeeos_biquintic_progress/experiments/038_solver_manifest/queries_v2.json"
        cases = json.loads(manifest.read_text())["cases"]
        cls.points = [(case["logRho"], case["logT"]) for case in cases]
        with (WORKSPACE / "freeeos_biquintic_progress/experiments/041_consumer_path_selection/states/classified_paths.csv").open() as stream:
            blend_points = sorted({(float(row["logRho"]), float(row["logT"]))
                                   for row in csv.DictReader(stream) if row["classification"] == "blend"})
        cls.points += blend_points[:3] + [(4., 6.)]
        cls.records = [(route, *point) for point in cls.points for route in ("component", "final")]
        cls.baseline = cls.run_driver("baseline_driver", [f"{route} {rho:.17g} {temp:.17g}" for route, rho, temp in cls.records])
        cls.off = cls.run_driver("hook_driver", [f"off {route} {rho:.17g} {temp:.17g}" for route, rho, temp in cls.records], "off")
        cls.passthrough = cls.run_driver("hook_driver", [f"passthrough {route} {rho:.17g} {temp:.17g}" for route, rho, temp in cls.records], "passthrough")

    @classmethod
    def run_driver(cls, executable, rows, name=None):
        name = name or executable
        command = [str(EXPERIMENT / "build/driver" / executable), str(MESA), str(INLIST), str(MASSFILE)]
        query = "\n".join(rows) + "\n"
        (cls.output / (name + ".queries.txt")).write_text(query)
        result = subprocess.run(command, input=query, text=True, capture_output=True,
                                env=cls.environment, cwd=cls.output, timeout=120)
        (cls.output / (name + ".stdout")).write_text(result.stdout)
        (cls.output / (name + ".stderr")).write_text(result.stderr)
        if result.returncode:
            raise RuntimeError("Driver failure: " + result.stderr[-2000:])
        data = np.loadtxt(io.StringIO(result.stdout), ndmin=2)
        if data.shape != (len(rows), 125):
            raise ValueError("Unexpected driver protocol shape")
        np.testing.assert_array_equal(data[:, 0], np.arange(1, len(rows) + 1))
        return data

    def test_opt_out_matches_unchanged_installed_eos(self):
        np.testing.assert_array_equal(self.baseline[:, 1], self.off[:, 1])
        success = self.baseline[:, 1] == 0
        expected = np.zeros(len(self.records))
        expected[-2] = -1
        np.testing.assert_array_equal(self.baseline[:, 1], expected)
        np.testing.assert_array_equal(self.baseline[success, 4:], self.off[success, 4:])
        np.testing.assert_array_equal(self.off[:, 2], 0)

    def test_passthrough_preserves_full_outputs_and_hybrid_rows(self):
        np.testing.assert_array_equal(self.off[:, 1], self.passthrough[:, 1])
        success = self.off[:, 1] == 0
        np.testing.assert_array_equal(self.off[success, 4:], self.passthrough[success, 4:])
        np.testing.assert_array_equal(self.passthrough[:, 3], 1)
        component = np.array([route == "component" for route, _, _ in self.records])
        np.testing.assert_array_equal(self.passthrough[component & success, 2], 1)
        np.testing.assert_array_equal(self.passthrough[~success, 2], 0)

    def test_errors_propagate_through_both_apis(self):
        rho, temp = self.points[0]
        rows = [f"{mode} {route} {rho:.17g} {temp:.17g}"
                for mode in ("missing", "unsafe", "failure") for route in ("component", "final")]
        data = self.run_driver("hook_driver", rows, "errors")
        np.testing.assert_array_equal(data[:, 1], [-7901, -7901, -7902, -7902, -7999, -7999])
        np.testing.assert_array_equal(data[:, 2], [0, 0, 0, 0, 1, 1])

    def test_hook_precedes_existing_component_blending(self):
        rows = [f"inject {route} {rho:.17g} {temp:.17g}" for route, rho, temp in self.records]
        injected = self.run_driver("hook_driver", rows, "inject")
        np.testing.assert_array_equal(injected[:, 1], self.off[:, 1])
        cv_index, fraction_index = 4 + 10, 4 + 21
        blend_count, inactive_count = 0, 0
        for offset in range(0, len(self.records), 2):
            component, final = offset, offset + 1
            fraction = self.off[final, fraction_index]
            if self.off[component, 1] != 0:
                self.assertEqual(fraction, 0)
                self.assertEqual(injected[final, 2], 0)
                np.testing.assert_array_equal(injected[final, 4:], self.off[final, 4:])
                inactive_count += 1
                continue
            self.assertEqual(injected[component, cv_index], 2 * self.off[component, cv_index])
            expected = self.off[final, cv_index] + fraction * self.off[component, cv_index]
            np.testing.assert_allclose(injected[final, cv_index], expected, rtol=5e-14)
            blend_count += 0 < fraction < 1
            inactive_count += fraction == 0
        self.assertGreater(blend_count, 0)
        self.assertGreater(inactive_count, 0)

    @classmethod
    def tearDownClass(cls):
        verify_private_cache(cls.output)
        after = installed_inventory()
        if cls.before != after:
            raise AssertionError("Installed cache or EOS library changed")
        record = {"scope": "Hook transport/configuration and original-EOS parity; no potential adapter or stellar test",
                  "query_pairs": len(cls.points), "before": cls.before, "after": after,
                  "installed_unchanged": True,
                  "script_sha256": sha256(__file__),
                  "hook_driver_sha256": sha256(EXPERIMENT / "build/driver/hook_driver"),
                  "baseline_driver_sha256": sha256(EXPERIMENT / "build/driver/baseline_driver")}
        (cls.output / "provenance.json").write_text(json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    unittest.main()
