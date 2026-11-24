"""Synthetic format/provenance tests; not a physical source validation."""

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from export_cache import export_cache
from freeeos_potential import checksum, read_table


class CacheExportTests(unittest.TestCase):
    def write_source(self, path):
        logq, logt = np.linspace(-1., -.6, 9), np.linspace(3.2, 3.6, 9)
        lines = ["synthetic format fixture", "51 .7 .02 9 3.2 3.6 .05 9 -1 -.6 .05"]
        for coordinate_q in logq:
            for coordinate_t in logt:
                row = np.ones(22)
                row[0], row[1], row[2], row[3] = coordinate_t, 7., 13., 9.
                row[4], row[5], row[6], row[7] = 1., 1., 1.2e8, 1.e8
                row[8], row[9], row[10] = 0., 1.e8 / 10**coordinate_t, -1.e8
                row[17] = 0
                row[18] = coordinate_q + 2 * coordinate_t - 12
                lines.append(" ".join(f"{value:.17e}" for value in row))
        path.write_text("\n".join(lines) + "\n")

    def test_export_metadata_records_and_invalid_zero_fill(self):
        with tempfile.TemporaryDirectory() as directory:
            source, cache = Path(directory) / "source.data", Path(directory) / "prototype.cache"
            self.write_source(source)
            result = export_cache(source, cache, radius=1)
            lines = cache.read_text().splitlines()
            self.assertEqual(lines[0], "FREEEOS_BIQUINTIC_V1")
            self.assertEqual(lines[1], "9 9")
            self.assertEqual(lines[3], checksum(source))
            self.assertEqual(lines[4].split()[1:], ["1", "mean"])
            records = np.stack([np.fromstring(line, sep=" ") for line in lines[6:]])
            self.assertEqual(records.shape, (81, 12))
            np.testing.assert_array_equal(records[:, 0], np.tile(np.arange(1, 10), 9))
            np.testing.assert_array_equal(records[:, 1], np.repeat(np.arange(1, 10), 9))
            np.testing.assert_array_equal(records[records[:, 2] == 0, 3:], 0)
            self.assertEqual(int(records[:, 2].sum()), result["valid_nodes"])
            self.assertTrue(np.all(np.isfinite(records)))
            self.assertEqual(json.loads(cache.with_suffix(".cache.json").read_text())["cache_sha256"], checksum(cache))

    def test_header_rejection(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "source.data"
            self.write_source(source)
            original = source.read_text()
            for header in ("51 .7 .02 9 3.2 3.6 0 9 -1 -.6 .05", "51 .7 .02 1 3.2 3.6 .05 9 -1 -.6 .05", "50 .7 .02 9 3.2 3.6 .05 9 -1 -.6 .05"):
                lines = original.splitlines()
                lines[1] = header
                source.write_text("\n".join(lines))
                with self.assertRaises(ValueError):
                    read_table(source)


if __name__ == "__main__":
    unittest.main()
