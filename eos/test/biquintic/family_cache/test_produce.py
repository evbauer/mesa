"""Audited reader regression cases for source-grid metadata and coordinates."""

import tempfile
import unittest
from pathlib import Path

import numpy as np

from produce import read_audited


class AuditedReaderTests(unittest.TestCase):
    def make_source(self, path, shifted_density=False, short_block=False):
        density_axis = np.linspace(-1.0, -0.2, 9)
        temperature_axis = np.linspace(3.2, 3.8, 7)
        lines = ["synthetic source", "51 .6 .02 7 3.2 3.8 .1 8 -1 -.3 .1"]
        for density_index, coordinate_q in enumerate(density_axis):
            lines.extend(("logQ = logRho - 2*logT + 12", f"{coordinate_q:.17e}"))
            for temperature_index, coordinate_t in enumerate(temperature_axis):
                if short_block and density_index == 4 and temperature_index == 2:
                    continue
                record = np.ones(22)
                record[0] = coordinate_t
                record[18] = coordinate_q + 2 * coordinate_t - 12
                if shifted_density and density_index == 4 and temperature_index == 2:
                    record[18] += 1e-4
                lines.append(" ".join(f"{value:.17e}" for value in record))
        path.write_text("\n".join(lines) + "\n")

    def test_reconstructs_redundant_header_only_from_complete_grid(self):
        with tempfile.TemporaryDirectory() as directory:
            source_path = Path(directory) / "mesa-FreeEOS_02z60x.data"
            self.make_source(source_path)
            observed, logq, logt, rows, audit = read_audited(source_path)
            self.assertEqual(rows.shape, (9, 7, 22))
            self.assertEqual((len(logq), len(logt)), (9, 7))
            self.assertEqual(observed[7], 9)
            self.assertAlmostEqual(observed[9], -.2)
            self.assertEqual(set(audit["header_anomalies"]), {"num_logQ", "logQ_max"})

    def test_mismatched_density_is_audited_without_modifying_record(self):
        with tempfile.TemporaryDirectory() as directory:
            source_path = Path(directory) / "mesa-FreeEOS_02z60x.data"
            self.make_source(source_path, shifted_density=True)
            _, logq, logt, rows, audit = read_audited(source_path)
            self.assertEqual(audit["density_coordinate_mismatch_count"], 1)
            self.assertAlmostEqual(rows[4, 2, 18] - (logq[4] + 2 * logt[2] - 12), 1e-4)

    def test_incomplete_printed_rectangle_is_unavailable(self):
        with tempfile.TemporaryDirectory() as directory:
            source_path = Path(directory) / "mesa-FreeEOS_02z60x.data"
            self.make_source(source_path, short_block=True)
            with self.assertRaisesRegex(ValueError, "rectangle"):
                read_audited(source_path)


if __name__ == "__main__":
    unittest.main()
