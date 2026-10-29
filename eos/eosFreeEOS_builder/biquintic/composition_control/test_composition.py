"""Polynomial, radiation and gauge tests for the offline composition sidecar."""

import math
from io import StringIO
import unittest
import warnings
from pathlib import Path
from unittest.mock import mock_open, patch

import numpy as np

from composition import SUPPORTS, combine, composition_partials, composition_weights, gauge_jet
from endpoint_audit import HEADER_FIELDS, audit_printed_grid, comparison_header, read_endpoint_table
from freeeos_potential import (
    CRAD, EVAL_ORDERS, LN10, construct_twists, evaluate, radiation_potential,
    read_table as core_read_table, source_constraints as core_source_constraints, thermodynamics,
)
from run_control import SourceGridError, column_control, header_endpoint_failures, load_sources, printed_state


class EndpointMetadataTests(unittest.TestCase):
    def setUp(self):
        self.header = np.array([51, .6, .02, 3, 3.2, 3.24, .02, 3, -1, -.94, .03])
        self.logq = self.header[8] + np.arange(int(self.header[7])) * self.header[10]
        self.logt = self.header[4] + np.arange(int(self.header[3])) * self.header[6]

    def printed_stream(self, labels=None, counts=None, density_delta=None, temperature_delta=None):
        labels = self.logq if labels is None else labels
        counts = [len(self.logt)] * len(labels) if counts is None else counts
        lines = ["format-51 header", " ".join(format(value, ".17g") for value in self.header)]
        for density_index, label in enumerate(labels):
            lines.extend(["logQ = logRho - 2*logT + 12", format(label, ".17g")])
            for temperature_index in range(counts[density_index]):
                temperature = self.logt[min(temperature_index, len(self.logt) - 1)]
                row = np.zeros(22)
                row[0] = temperature + (temperature_delta or {}).get((density_index, temperature_index), 0)
                row[18] = self.logq[density_index] + 2 * temperature - 12
                row[18] += (density_delta or {}).get(density_index, 0)
                row[4], row[7] = 1, 1
                lines.append(" ".join(format(value, ".17g") for value in row))
        return StringIO("\n".join(lines) + "\n")

    def test_redundant_display_maxima_require_independent_grid_evidence(self):
        strict_header = comparison_header(self.header)
        authoritative_header = comparison_header(self.header, True)
        self.header[9], self.header[5] = 3.69, 8.2
        self.assertFalse(np.array_equal(strict_header, comparison_header(self.header)))
        np.testing.assert_array_equal(authoritative_header, comparison_header(self.header, True))
        audit = audit_printed_grid(self.printed_stream(), self.header, self.logq, self.logt)
        self.assertEqual(audit["status"], "PASSED")
        self.assertEqual(audit["printed_record_count"], 9)
        self.assertEqual(audit["Q_block_count"], 3)
        self.assertAlmostEqual(audit["derived_logQ_endpoint"], -.94)

    def test_header_comparison_ignores_only_named_display_maxima(self):
        original = comparison_header(self.header, True)
        for index, name in enumerate(HEADER_FIELDS):
            altered = self.header.copy()
            altered[index] += 1
            matches = np.array_equal(original, comparison_header(altered, True))
            self.assertEqual(matches, name in ("X", "Z", "logQ_max", "logT_max"))

    def test_reject_actual_grid_count_min_or_step_change(self):
        for index in (3, 4, 6, 7, 8, 10):
            altered = self.header.copy()
            altered[index] += 1
            with self.assertRaisesRegex(ValueError, "Actual grid/count/min/step differs"):
                audit_printed_grid(self.printed_stream(), altered, self.logq, self.logt)

    def test_reject_printed_count_discrepancy(self):
        for counts in ([3, 3, 2], [3, 3, 4]):
            with self.assertRaisesRegex(ValueError, "record count"):
                audit_printed_grid(self.printed_stream(counts=counts), self.header, self.logq, self.logt)
        with self.assertRaisesRegex(ValueError, "label/record count"):
            audit_printed_grid(self.printed_stream(labels=self.logq[:2]), self.header, self.logq, self.logt)

    def test_reject_actual_Q_labels_or_nonmonotonic_blocks(self):
        for labels in (self.logq + .001, self.logq[::-1]):
            with self.assertRaisesRegex(ValueError, "Q-block label differs"):
                audit_printed_grid(self.printed_stream(labels=labels), self.header, self.logq, self.logt)

    def test_reject_printed_T_or_rho_mismatch(self):
        with self.assertRaisesRegex(ValueError, "Printed T differs"):
            audit_printed_grid(self.printed_stream(temperature_delta={(1, 1): .001}), self.header, self.logq, self.logt)
        with self.assertRaisesRegex(ValueError, "Printed rho differs"):
            audit_printed_grid(self.printed_stream(density_delta={1: -.0001}), self.header, self.logq, self.logt)

    def test_known_coordinate_exclusion_is_narrow_and_not_a_grid_repair(self):
        self.header[1] = .7
        self.header[7:11] = [527, -10.09, 5.69, .03]
        self.logq = self.header[8] + np.arange(527) * self.header[10]
        audit = audit_printed_grid(self.printed_stream(density_delta={520: -.0001}), self.header, self.logq, self.logt)
        exclusion = audit["known_coordinate_exclusions"][0]
        self.assertEqual(exclusion["iq"], 520)
        self.assertEqual(exclusion["record_count"], len(self.logt))
        for density_index, shift in ((519, -.0001), (520, -.001)):
            with self.assertRaisesRegex(ValueError, "Printed rho differs"):
                audit_printed_grid(self.printed_stream(density_delta={density_index: shift}), self.header, self.logq, self.logt)
        self.header[1] = .6
        with self.assertRaisesRegex(ValueError, "Printed rho differs"):
            audit_printed_grid(self.printed_stream(density_delta={520: -.0001}), self.header, self.logq, self.logt)

    def test_cross_source_actual_grids_and_counts_still_fail_with_flag(self):
        for changed_field in (7, 8, 10):
            def read_source(source, validate_endpoints=True):
                coordinate = int(source.name.split("z")[1].split("x")[0]) / 100
                header = self.header.copy()
                header[1] = coordinate
                if coordinate == .6:
                    header[changed_field] += 1 if changed_field == 7 else .01
                logq = header[8] + np.arange(int(header[7])) * header[10]
                logt = header[4] + np.arange(int(header[3])) * header[6]
                header[9] = logq[-1]
                rows = np.zeros((len(logq), len(logt), 22))
                rows[:, :, 0] = logt[None, :]
                rows[:, :, 18] = logq[:, None] + 2 * logt[None, :] - 12
                rows[:, :, 4] = 1
                rows[:, :, 7] = 1
                return header, logq, logt, rows

            with patch("run_control.read_table", side_effect=read_source), patch("run_control.checksum", return_value="unchanged"), \
                    patch("pathlib.Path.open", mock_open()), patch("run_control.audit_printed_grid", return_value={"known_coordinate_exclusions": []}):
                with self.assertRaises(SourceGridError) as caught:
                    load_sources(Path("unused"), .02, True, allow_audited_endpoint_metadata=True)
            failures = caught.exception.failures
            self.assertTrue(any(entry["kind"] == "header_or_coordinate_grid_differs_from_X_050" for entry in failures))

    def test_flag_does_not_enable_endpoint_override_at_other_Z(self):
        def read_source(source, validate_endpoints=True):
            self.assertTrue(validate_endpoints)
            coordinate = int(source.name.split("z")[1].split("x")[0]) / 100
            header = self.header.copy()
            header[1:3] = [coordinate, .04]
            if coordinate == .6:
                header[9] = 3.69
            rows = np.zeros((len(self.logq), len(self.logt), 22))
            rows[:, :, 0] = self.logt[None, :]
            rows[:, :, 18] = self.logq[:, None] + 2 * self.logt[None, :] - 12
            rows[:, :, 4] = 1
            rows[:, :, 7] = 1
            return header, self.logq, self.logt, rows

        with patch("run_control.read_table", side_effect=read_source), patch("run_control.checksum", return_value="unchanged"), \
                patch("run_control.audit_printed_grid") as audit:
            with self.assertRaises(SourceGridError):
                load_sources(Path("unused"), .04, True, allow_audited_endpoint_metadata=True)
            audit.assert_not_called()

    def test_independent_loader_is_strict_by_default_and_explicit_when_audited(self):
        source = Path("in_memory")
        self.header[9] = 3.69
        text = self.printed_stream().getvalue()
        with patch("pathlib.Path.open", mock_open(read_data=text)):
            with self.assertRaisesRegex(ValueError, "Strict independent loader rejects inconsistent logQ_max"):
                read_endpoint_table(source)
            header, logq, logt, rows, audit = read_endpoint_table(source, allow_audited_endpoint_metadata=True)
        self.assertEqual(header[9], 3.69)
        np.testing.assert_array_equal(logq, self.logq)
        np.testing.assert_array_equal(logt, self.logt)
        self.assertEqual(rows.shape, (3, 3, 22))
        self.assertTrue(audit["allow_audited_endpoint_metadata"])

    def test_independent_loader_never_waives_counts_or_other_Z(self):
        source = Path("in_memory")
        text = self.printed_stream(counts=[3, 3, 2]).getvalue()
        with patch("pathlib.Path.open", mock_open(read_data=text)):
            with self.assertRaisesRegex(ValueError, "record count differs"):
                read_endpoint_table(source, allow_audited_endpoint_metadata=True)
        self.header[2] = .04
        text = self.printed_stream().getvalue()
        with patch("pathlib.Path.open", mock_open(read_data=text)):
            with self.assertRaisesRegex(ValueError, "authorized only at Z=.02"):
                read_endpoint_table(source, allow_audited_endpoint_metadata=True)

    def test_audited_core_keyword_full_flow_requires_audit_before_constraints(self):
        events = []

        def open_source(source, *arguments, **keywords):
            coordinate = int(source.name.split("z")[1].split("x")[0]) / 100
            self.header[1] = coordinate
            self.header[9] = 3.69 if coordinate == .6 else self.logq[-1]
            return self.printed_stream()

        def audited_grid(stream, header, logq, logt):
            result = audit_printed_grid(stream, header, logq, logt)
            events.append(("audit_passed", header[1]))
            return result

        def constraints(rows, logq, logt, split_radiation=False):
            self.assertEqual(events[-1], ("audit_passed", self.header[1]))
            result = core_source_constraints(rows, logq, logt, split_radiation=split_radiation)
            events.append(("constraints_used_rows", self.header[1]))
            return result

        with patch("pathlib.Path.open", open_source), patch("run_control.checksum", return_value="unchanged"), \
                patch("run_control.read_table", wraps=core_read_table) as core_reader, \
                patch("run_control.audit_printed_grid", side_effect=audited_grid), \
                patch("run_control.leadsource_constraints", side_effect=constraints):
            _, _, _, _, _, _, _, metadata = load_sources(Path("unused"), .02, True, allow_audited_endpoint_metadata=True)
            self.assertEqual(core_reader.call_count, 5)
            self.assertTrue(all(call.kwargs == {"validate_endpoints": False} for call in core_reader.call_args_list))
        self.assertEqual(events, [(stage, coordinate) for coordinate in (.5, .6, .7, .8, .9)
                                  for stage in ("audit_passed", "constraints_used_rows")])
        self.assertEqual(metadata["0.6"]["header"][9], 3.69)
        for entry in metadata.values():
            self.assertEqual(entry["source_reader"], "lead_reader_explicit_audited_endpoint_option")
            self.assertFalse(entry["core_validate_endpoints"])
            self.assertEqual(entry["endpoint_metadata_audit"]["status"], "PASSED")

    def test_failed_printed_grid_audit_cannot_use_source_rows(self):
        def open_source(source, *arguments, **keywords):
            coordinate = int(source.name.split("z")[1].split("x")[0]) / 100
            self.header[1] = coordinate
            self.header[9] = 3.69 if coordinate == .6 else self.logq[-1]
            return self.printed_stream(labels=self.logq + .001 if coordinate == .6 else self.logq)

        with patch("pathlib.Path.open", open_source), patch("run_control.checksum", return_value="unchanged"), \
                patch("run_control.read_table", wraps=core_read_table), \
                patch("run_control.leadsource_constraints", wraps=core_source_constraints) as constraints:
            with self.assertRaises(SourceGridError) as caught:
                load_sources(Path("unused"), .02, True, allow_audited_endpoint_metadata=True)
            self.assertEqual(constraints.call_count, 4)
        self.assertTrue(any(entry["X"] == .6 and entry["kind"] == "audited_endpoint_reader_rejected_source"
                            for entry in caught.exception.failures))


def polynomial_nodes(logq, logt):
    native = np.zeros((3, 3, len(logq), len(logt)))
    density_coordinate = logq[:, None] * LN10
    temperature_coordinate = logt[None, :] * LN10
    for density_order in range(3):
        for temperature_order in range(3):
            native[density_order, temperature_order] = (
                1e7 * math.factorial(4) / math.factorial(4 - density_order) * density_coordinate**(4 - density_order)
                * math.factorial(5) / math.factorial(5 - temperature_order) * temperature_coordinate**(5 - temperature_order)
            )
    return native


class CompositionTests(unittest.TestCase):
    def setUp(self):
        self.logq, self.logt = np.array([.2, .5]), np.array([3.3, 3.4])
        self.complete = np.ones((2, 2), bool)
        self.nodes = polynomial_nodes(self.logq, self.logt)
        self.sample_q, self.sample_t = np.array([.27, .38, .45]), np.array([3.32, 3.35, 3.38])

    def test_expected_weights_and_negative_warning(self):
        linear, _ = composition_weights("linear")
        np.testing.assert_allclose(list(linear.values()), [.5, .5], atol=1e-15)
        with self.assertWarnsRegex(RuntimeWarning, "NEGATIVE"):
            cubic, _ = composition_weights("cubic_lagrange")
        np.testing.assert_allclose(list(cubic.values()), [-1/6, 2/3, 2/3, -1/6], atol=1e-15)
        self.assertNotIn(.7, cubic)
        self.assertNotIn(.7, linear)

    def test_composition_polynomial_values_and_partials(self):
        for scheme, support in SUPPORTS.items():
            for target in (.61, .7, .79):
                weights, partials = composition_weights(scheme, target, warn=False)
                self.assertAlmostEqual(sum(weights.values()), 1)
                self.assertAlmostEqual(sum(partials.values()), 0)
                for power in range(len(support)):
                    values = {coordinate: coordinate**power for coordinate in support}
                    self.assertAlmostEqual(combine(values, weights), target**power)
                    expected = power * target**(power - 1) if power else 0
                    self.assertAlmostEqual(combine(values, partials), expected)

    def test_bounded_no_extrapolation(self):
        for target in (.59, .81, np.nan):
            with self.assertRaises(ValueError):
                composition_weights("linear", target)

    def test_centered_two_stage_temperature_coverage_starts_at_328(self):
        logt = 3.1 + np.arange(29) * .02
        valid = np.broadcast_to(logt[None, :] >= 3.2 - 1e-12, (29, 29)).copy()
        low = np.zeros((3, 3, 29, 29))
        _, complete, _ = construct_twists(low, valid, .03 * LN10, .02 * LN10, radius=2, route="mean")
        self.assertAlmostEqual(logt[np.any(complete, axis=0)].min(), 3.28)
        self.assertFalse(np.any(complete[:, logt < 3.28 - 1e-12]))

    def test_reject_inconsistent_header_endpoint(self):
        header = np.array([51, .6, .02, 306, 2.1, 8.2, .02, 527, -10.09, 5.69, .03])
        logq, logt = header[8] + np.arange(527) * .03, header[4] + np.arange(306) * .02
        self.assertEqual(header_endpoint_failures(header, logq, logt), [])
        header[9] = 3.69
        failures = header_endpoint_failures(header, logq, logt)
        self.assertEqual(failures[0]["field"], "logQ_max")
        self.assertAlmostEqual(failures[0]["declared"], 3.69)
        self.assertAlmostEqual(failures[0]["count_start_step_endpoint"], 5.69)

    def test_tensor_polynomial_nodal_mix_and_jets(self):
        base = evaluate(self.nodes, self.complete, self.logq, self.logt, self.sample_q, self.sample_t)
        for scheme, support in SUPPORTS.items():
            weights, _ = composition_weights(scheme, warn=False)
            for power in range(len(support)):
                source_nodes = {coordinate: self.nodes * coordinate**power for coordinate in support}
                mixed = evaluate(combine(source_nodes, weights), self.complete, self.logq, self.logt,
                                 self.sample_q, self.sample_t)
                source_jets = {coordinate: evaluate(nodes, self.complete, self.logq, self.logt,
                                                    self.sample_q, self.sample_t)
                               for coordinate, nodes in source_nodes.items()}
                composed = combine(source_jets, weights)
                for order in EVAL_ORDERS:
                    np.testing.assert_allclose(mixed[order], base[order] * .7**power, rtol=5e-9)
                    np.testing.assert_allclose(mixed[order], composed[order], rtol=5e-9)

    def test_radiation_partition_unity_and_zero_X_partial(self):
        radiation = radiation_potential(self.sample_q, self.sample_t)
        for scheme, support in SUPPORTS.items():
            weights, partials = composition_weights(scheme, warn=False)
            mixed = combine(dict.fromkeys(support, radiation), weights)
            derivative = combine(dict.fromkeys(support, radiation), partials)
            for order in EVAL_ORDERS:
                np.testing.assert_allclose(mixed[order], radiation[order], rtol=5e-15)
                np.testing.assert_allclose(derivative[order] / radiation[order], 0, atol=1e-14)

    def test_exact_composition_gauges_fixed_X_invariants(self):
        base = evaluate(self.nodes, self.complete, self.logq, self.logt, self.sample_q, self.sample_t)
        baseline = thermodynamics(base, self.sample_q, self.sample_t)
        for scheme, support in SUPPORTS.items():
            weights, _ = composition_weights(scheme, warn=False)
            constants = {coordinate: 1e10 * coordinate**2 for coordinate in support}
            slopes = {coordinate: 1e6 * (1 + coordinate**3) for coordinate in support}
            for use_constant, use_slope in ((True, False), (False, True), (True, True)):
                shifts = {coordinate: gauge_jet(self.sample_q, self.sample_t,
                                               constants[coordinate] if use_constant else 0,
                                               slopes[coordinate] if use_slope else 0) for coordinate in support}
                shifted = thermodynamics(base + combine(shifts, weights), self.sample_q, self.sample_t)
                for field in ("P", "Cv"):
                    np.testing.assert_allclose(shifted[field], baseline[field], rtol=5e-13)
                expected_e = combine(constants, weights) if use_constant else 0
                expected_s = -combine(slopes, weights) if use_slope else 0
                np.testing.assert_allclose(shifted["E"], baseline["E"] + expected_e, rtol=5e-13, atol=1e-4)
                np.testing.assert_allclose(shifted["S"], baseline["S"] + expected_s, rtol=5e-13, atol=1e-7)

    def test_nodal_linear_T_gauge_is_exact_at_knots_not_offnode(self):
        temperature_grid = np.broadcast_to(self.logt[None, :], (2, 2))
        shifted_nodes = self.nodes + gauge_jet(None, temperature_grid, 1e10, 1e8)[:3, :3]
        node_q, node_t = np.meshgrid(self.logq, self.logt, indexing="ij")
        node_q, node_t = node_q.ravel(), node_t.ravel()
        baseline = thermodynamics(evaluate(self.nodes, self.complete, self.logq, self.logt, node_q, node_t), node_q, node_t)
        shifted = thermodynamics(evaluate(shifted_nodes, self.complete, self.logq, self.logt, node_q, node_t), node_q, node_t)
        for field in ("P", "Cv"):
            np.testing.assert_allclose(shifted[field], baseline[field], rtol=2e-12)
        baseline = thermodynamics(evaluate(self.nodes, self.complete, self.logq, self.logt, self.sample_q, self.sample_t),
                                  self.sample_q, self.sample_t)
        shifted = thermodynamics(evaluate(shifted_nodes, self.complete, self.logq, self.logt, self.sample_q, self.sample_t),
                                 self.sample_q, self.sample_t)
        self.assertGreater(np.max(abs(shifted["Cv"] - baseline["Cv"])), 1)
        np.testing.assert_allclose(shifted["P"], baseline["P"], rtol=1e-10)

    def test_same_physical_columns_not_log_values(self):
        rows = np.zeros((2, 22))
        rows[:, 1:4] = np.array([[2, 3, 4], [4, 5, 6]])
        temperature = np.array([3.3, 3.3])
        states = {.6: printed_state(rows[:1], temperature[:1]), .8: printed_state(rows[1:], temperature[1:])}
        weights, _ = composition_weights("linear")
        mixed = column_control(states, weights)
        self.assertAlmostEqual(mixed["Pgas"][0], (1e2 + 1e4) / 2)
        self.assertAlmostEqual(mixed["P"][0] - mixed["Pgas"][0], CRAD * 10**(4 * 3.3) / 3)
        self.assertNotEqual(mixed["Pgas"][0], 1e3)

    def test_composition_response_partials_of_same_potential(self):
        for scheme, support in SUPPORTS.items():
            models = {coordinate: (1 + coordinate**2) * self.nodes for coordinate in support}
            weights, partial_weights = composition_weights(scheme, warn=False)
            native = evaluate(combine(models, weights), self.complete, self.logq, self.logt,
                              self.sample_q, self.sample_t, add_radiation=True)
            partial_native = evaluate(combine(models, partial_weights), self.complete, self.logq, self.logt,
                                      self.sample_q, self.sample_t)
            state = thermodynamics(native, self.sample_q, self.sample_t)
            analytic = composition_partials(native, partial_native, state)
            step = 1e-4
            plus_weights, _ = composition_weights(scheme, .7 + step, warn=False)
            minus_weights, _ = composition_weights(scheme, .7 - step, warn=False)
            plus = thermodynamics(evaluate(combine(models, plus_weights), self.complete, self.logq, self.logt,
                                           self.sample_q, self.sample_t, add_radiation=True), self.sample_q, self.sample_t)
            minus = thermodynamics(evaluate(combine(models, minus_weights), self.complete, self.logq, self.logt,
                                            self.sample_q, self.sample_t, add_radiation=True), self.sample_q, self.sample_t)
            for field, partial in analytic.items():
                scale = np.maximum.reduce([np.ones(3), abs(state[field]), abs(partial)])
                np.testing.assert_allclose(((plus[field] - minus[field]) / (2 * step) - partial) / scale, 0, atol=2e-7)


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        unittest.main()
