"""Dense-audit masks, accounting, inherited diagnostics, and synthetic pipeline."""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

import dense_audit as audit
from freeeos_potential import CRAD, thermodynamic_partials, to_native, EVAL_ORDERS


def ideal_rows(logq, logt):
    temperature = np.broadcast_to(10.0**logt[None, :], (logq.size, logt.size))
    density = 10.0**(logq[:, None] + 2 * logt[None, :] - 12)
    gas_constant, heat_factor = 1.2e8, 1.5
    radiation_pressure = CRAD * temperature**4 / 3
    gas_pressure = gas_constant * density * temperature
    pressure = gas_pressure + radiation_pressure
    energy = heat_factor * gas_constant * temperature + 3 * radiation_pressure / density
    entropy = gas_constant * (heat_factor * np.log(temperature) + heat_factor - np.log(density))
    entropy += 4 * radiation_pressure / (density * temperature)
    cv = heat_factor * gas_constant + 12 * radiation_pressure / (density * temperature)
    rows = np.zeros((logq.size, logt.size, 22))
    rows[:, :, 0] = logt[None, :]
    rows[:, :, 1] = np.log10(gas_pressure)
    rows[:, :, 2] = np.log10(energy)
    rows[:, :, 3] = np.log10(entropy)
    rows[:, :, 4] = gas_pressure / pressure
    rows[:, :, 5] = (gas_pressure + 4 * radiation_pressure) / pressure
    rows[:, :, 7] = cv
    rows[:, :, 18] = np.log10(density)
    return rows


def saved_queries(logq, logt):
    density_index = np.array([6, 8, 10])
    temperature_index = np.array([6, 8, 10])
    phase_q = np.array([[.23], [.67]])
    phase_t = np.array([[.37], [.79]])
    sample_q = logq[density_index][None, :] + phase_q * (logq[1] - logq[0])
    sample_t = logt[temperature_index][None, :] + phase_t * (logt[1] - logt[0])
    shape = sample_q.shape
    data = {"logQ": sample_q, "logT": sample_t, "density_cell_index": density_index,
            "temperature_cell_index": temperature_index, "eligible": audit.broad_own_use(sample_q, sample_t),
            "material": np.array([[True, False, True], [True, True, False]]),
            "Cv_mismatch": np.full(shape, -2.0)}
    for route in ("mathematical", "physical"):
        for field in audit.IDENTITY_FIELDS:
            data[f"{route}_{field}"] = np.full(shape, .025 if route == "mathematical" else .0125)
    return data


class DenseAuditTests(unittest.TestCase):
    def test_default_output_is_progress_report_not_original_assessment(self):
        self.assertEqual(audit.DEFAULT_OUTPUT, audit.ROOT / "freeeos_biquintic_progress/experiments/004_dense_validation")
        self.assertNotIn("freeeos_report", audit.DEFAULT_OUTPUT.relative_to(audit.ROOT).parts)

    def setUp(self):
        build = Path(__file__).resolve().parent / "build"
        build.mkdir(exist_ok=True)
        self.temporary = tempfile.TemporaryDirectory(dir=build)
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.logq = np.linspace(-3, -.3, 19)
        self.logt = np.linspace(3.2, 5, 19)
        self.rows = ideal_rows(self.logq, self.logt)
        self.saved = saved_queries(self.logq, self.logt)
        self.archive = self.directory / "saved.npz"
        np.savez_compressed(self.archive, **self.saved)

    def test_distribution_counts_nonfinite_and_empty(self):
        result = audit.distribution([0, -2, 4, np.nan, np.inf])
        self.assertEqual((result["count"], result["finite"], result["nonfinite"]), (5, 3, 2))
        self.assertEqual(result["absolute"]["median"], 2)
        self.assertEqual(result["signed"]["min"], -2)
        self.assertIsNone(audit.distribution([])["absolute"])

    def test_saved_query_loading_and_strict_off_node_validation(self):
        data = audit.load_queries(self.archive, self.logq, self.logt)
        self.assertEqual(data["logQ"].shape, (2, 3))
        np.testing.assert_allclose(data["density_phase"][:, 0], [.23, .67])
        changed = {field: values.copy() for field, values in self.saved.items()}
        changed["logQ"][0, 0] = self.logq[changed["density_cell_index"][0]]
        changed["eligible"] = audit.broad_own_use(changed["logQ"], changed["logT"])
        np.savez_compressed(self.archive, **changed)
        with self.assertRaisesRegex(ValueError, "strictly off-node"):
            audit.load_queries(self.archive, self.logq, self.logt)

    def test_archive_rejects_bad_masks_indices_and_shapes(self):
        for field, replacement, expected in (
                ("eligible", ~self.saved["eligible"], "eligible flag"),
                ("material", self.saved["material"].astype(int), "boolean"),
                ("density_cell_index", self.saved["density_cell_index"] + 1, "indices disagree"),
                ("logT", self.saved["logT"][:1], "query shape")):
            with self.subTest(field=field):
                changed = dict(self.saved)
                changed[field] = replacement
                np.savez_compressed(self.archive, **changed)
                with self.assertRaisesRegex(ValueError, expected):
                    audit.load_queries(self.archive, self.logq, self.logt)
        changed = {field: values for field, values in self.saved.items() if field != "material"}
        np.savez_compressed(self.archive, **changed)
        with self.assertRaisesRegex(ValueError, "Missing"):
            audit.load_queries(self.archive, self.logq, self.logt)

    def test_domain_exclusions_account_for_all_queries_and_propagate_printed_mismatch(self):
        saved = audit.load_queries(self.archive, self.logq, self.logt)
        models = audit.build_models(self.rows, self.logq, self.logt)
        valid, result = audit.domain_summary(self.rows, self.logq, self.logt, saved, models)
        self.assertEqual(valid.sum(), 6)
        self.assertEqual(sum(result["exclusive_query_exclusions_in_order"].values()), 6)
        damaged = self.rows.copy()
        damaged[2, 6, 18] += 1e-4
        changed_models = audit.build_models(damaged, self.logq, self.logt)
        changed_valid, changed = audit.domain_summary(damaged, self.logq, self.logt, saved, changed_models)
        self.assertEqual(changed["source_printed_density_mismatch_nodes"], 1)
        self.assertEqual(changed["exclusive_query_exclusions_in_order"]["printed_density_corner_mismatch"], 0)
        self.assertLess(changed_valid.sum(), valid.sum())
        self.assertGreater(changed["exclusive_query_exclusions_in_order"]["outside_common_complete_twist_stencil"], 0)

    def test_chunks_are_identical_and_radius2_ideal_states_positive(self):
        saved = audit.load_queries(self.archive, self.logq, self.logt)
        models = audit.build_models(self.rows, self.logq, self.logt)
        valid, _ = audit.domain_summary(self.rows, self.logq, self.logt, saved, models)
        for model in models.values():
            first = audit.evaluate_chunks(model["nodes"], model["complete"], self.logq, self.logt,
                                          saved, valid, model["split_radiation"], 1)
            second = audit.evaluate_chunks(model["nodes"], model["complete"], self.logq, self.logt,
                                           saved, valid, model["split_radiation"], 5)
            for field in audit.COMPARISON_FIELDS:
                np.testing.assert_array_equal(first[0][field], second[0][field])
            for field in audit.POSITIVE_FIELDS:
                self.assertTrue(np.all(first[0][field] > 0))
            for residual in first[1]["returned_function"].values():
                self.assertLess(np.max(abs(residual)), 1e-12)
        with self.assertRaisesRegex(ValueError, "positive"):
            audit.evaluate_chunks(model["nodes"], model["complete"], self.logq, self.logt,
                                  saved, valid, model["split_radiation"], 0)

    def test_identity_zero_denominator_is_counted_not_silently_passed(self):
        sample_q, sample_t = np.array([-2.0]), np.array([4.0])
        density, temperature = 10**(sample_q + 2 * sample_t - 12), 10**sample_t
        physical = np.full((4, 4, 1), np.nan)
        for order_q, order_t in EVAL_ORDERS:
            if order_q == 0:
                gas = 1.2e8 * temperature * (np.log(density) - 1.5 * np.log(temperature) - 1.5 * order_t)
            elif order_q == 1:
                gas = 1.2e8 * temperature
            else:
                gas = np.zeros(1)
            physical[order_q, order_t] = gas
        state = thermodynamic_partials(to_native(physical, EVAL_ORDERS), sample_q, sample_t)
        for field in ("E_t", "S_t", "E_T", "S_T"):
            state[field][:] = 0
        diagnostics = audit.identity_diagnostics(state)
        summary = audit.summarize_identities(diagnostics, {"all_valid": np.array([True])})
        self.assertEqual(summary["all_valid"]["returned_function"]["dse"]["nonfinite"], 1)
        self.assertEqual(summary["all_valid"]["well_conditioned_returned_function"]["dse"]["conditioning_excluded"], 1)

    def test_cv_mismatch_never_becomes_candidate_departure(self):
        values = {field: np.full(6, 3.0) for field in audit.COMPARISON_FIELDS}
        masks = audit.subset_masks(np.ones(6, bool), self.saved)
        context = audit.legacy_context(self.saved, masks, values)
        self.assertEqual(context["candidate_relative_departure_from_saved_legacy_states_NOT_source_accuracy"], {})
        self.assertIn("Cv", context["unavailable_saved_legacy_state_fields"])
        self.assertEqual(context["saved_diagnostics_on_identical_accepted_queries"]["all_valid"][
            "Cv_mismatch_LEGACY_function_vs_LEGACY_column_ONLY"]["signed"]["median"], -2)
        with_reference = dict(self.saved, legacy_Cv=np.full((2, 3), 2.0))
        changed = audit.legacy_context(with_reference, masks, values)
        departure = changed["candidate_relative_departure_from_saved_legacy_states_NOT_source_accuracy"]["Cv"]
        self.assertEqual(departure["subsets"]["all_valid"]["signed"]["median"], .5)
        self.assertEqual(departure["saved_reference_field"], "legacy_Cv")
        np.testing.assert_array_equal(audit.relative_departure(np.array([3., 2.]), np.array([2., 0.])), [.5, np.nan])

    def test_positivity_accounting_and_all_failure_coordinates(self):
        values = {field: np.ones(6) for field in audit.COMPARISON_FIELDS}
        values["Cv"] = np.array([-2., 0., np.nan, 1., 2., 3.])
        masks = audit.subset_masks(np.ones(6, bool), self.saved)
        result = audit.positivity_summary(values, masks, self.saved)["all_valid"]
        counts = result["fields"]["Cv"]
        self.assertEqual((counts["finite"], counts["positive_finite"], counts["nonpositive_finite"], counts["nonfinite"]), (5, 3, 2, 1))
        self.assertEqual((result["any_field_failed"], result["all_five_finite_positive"]), (3, 3))
        self.assertEqual(counts["minimum_finite"]["flat_index"], 0)
        self.assertEqual(counts["failure_examples"][2]["nonfinite_kind"], "NaN")
        audit.save_failures(self.directory, {"test": (values, {})}, masks, self.saved)
        with np.load(self.directory / "positivity_failures.npz") as failures:
            np.testing.assert_array_equal(failures["test__Cv__flat_index"], [0, 1, 2])
            np.testing.assert_array_equal(failures["test__Cv__logQ"], self.saved["logQ"].ravel()[:3])

    def test_cell_maps_use_selected_query_medians_and_failure_maxima(self):
        shape = (2, 2)
        values = np.array([-2., 4., -6., 8.])
        mask = np.array([True, False, True, True])
        np.testing.assert_array_equal(audit.cell_median(values, mask, shape, absolute=True), [4., 8.])
        maps = audit.cell_failure_maps(values, mask, shape)
        np.testing.assert_array_equal(maps["failed_count"], [2, 0])
        np.testing.assert_array_equal(maps["negative_max_relative_to_cell_median_absolute"], [1.5, np.nan])

    def test_synthetic_end_to_end_and_selected_twist_queries(self):
        table = self.directory / "ideal.data"
        header = [51, .7, .02, self.logt.size, self.logt[0], self.logt[-1], self.logt[1] - self.logt[0],
                  self.logq.size, self.logq[0], self.logq[-1], self.logq[1] - self.logq[0]]
        with table.open("w") as stream:
            stream.write("synthetic ideal gas plus radiation\n")
            stream.write(" ".join(map(str, header)) + "\n")
            np.savetxt(stream, self.rows.reshape(-1, 22), fmt="%.16e")
        output = self.directory / "output"
        result = audit.run_audit(table, self.archive, output, chunk_size=2, suspicious_cells=1,
                                 figures=False, query_archive=True)
        self.assertEqual(result["subset_counts"]["all_valid"], 6)
        self.assertEqual(result["twist_sensitivity"]["selected_saved_queries"], 2)
        self.assertEqual(result["legacy_context"]["candidate_relative_departure_from_saved_legacy_states_NOT_source_accuracy"], {})
        for variant in result["variants"].values():
            self.assertEqual(variant["positivity"]["all_valid"]["any_field_failed"], 0)
        saved_result = json.loads((output / "results.json").read_text())
        self.assertEqual(saved_result["domain"]["strictly_off_node_queries"], 6)
        sensitivity = json.loads((output / "suspicious_cell_twist_sensitivity.json").read_text())
        self.assertEqual(len(sensitivity["schemes"]), 8)
        with np.load(output / "query_diagnostics.npz") as diagnostics:
            self.assertEqual(diagnostics["radius2_mean_total__Cv"].shape, (2, 3))
        with np.load(output / "cell_maps.npz") as maps:
            self.assertEqual(maps["radius2_mean_total__all_valid__Cv__selected_count"].sum(), 6)

    def test_cli_rejects_invalid_chunk_size_without_outputs(self):
        result = subprocess.run([sys.executable, "-B", str(Path(audit.__file__)), "--chunk-size", "0",
                                 "--output", str(self.directory / "bad_output")], capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("chunk size must be positive", result.stderr)
        self.assertFalse((self.directory / "bad_output").exists())


if __name__ == "__main__":
    unittest.main()
