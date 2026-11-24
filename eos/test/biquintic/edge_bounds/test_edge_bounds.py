"""Resource/mask/native-floor and polynomial/restoration guard checks."""

import unittest
import importlib.util
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np

from run_control import (DTYPE, adaptive_split, bernstein_bounds, check_limits, finish_fence, guarded_powers,
                         new_cell_mask, polynomial_value, radiation_cv_numerator, radiation_guard, support_guard)

verifier_spec = importlib.util.spec_from_file_location('edge_bounds_recorded_verifier', Path(__file__).with_name('verify_artifacts.py'))
verifier = importlib.util.module_from_spec(verifier_spec)
verifier_spec.loader.exec_module(verifier)
coordinate_population_pair = verifier.coordinate_population_pair


class EdgeBoundTests(unittest.TestCase):
    def test_original_coordinate_pairing_reorders_and_rejects_duplicates(self):
        first = np.array([[4, 2], [1, 3], [9, 1]])
        second = np.array([[7, 0], [1, 3], [4, 2]])
        common_first, common_second, only_first, only_second = coordinate_population_pair(first, second)
        np.testing.assert_array_equal(first[common_first], second[common_second])
        np.testing.assert_array_equal(first[only_first], [[9, 1]])
        np.testing.assert_array_equal(second[only_second], [[7, 0]])
        for invalid in (np.array([[1, 2], [1, 2]]), np.array([[-1, 2]]), np.array([[1., 2.]])):
            with self.assertRaises(ValueError):
                coordinate_population_pair(invalid, second)

    def test_reproduction_fence_moves_daily_but_preserves_today_deadline(self):
        today = datetime(2026, 10, 28, 19, 0, tzinfo=ZoneInfo('US/Pacific'))
        tomorrow = datetime(2026, 10, 29, 19, 0, tzinfo=ZoneInfo('US/Pacific'))
        self.assertEqual(finish_fence(today).hour, 20)
        self.assertEqual(finish_fence(tomorrow).day, 29)
        self.assertEqual(finish_fence(today, '2026-10-28T21:00:00-07:00').hour, 20)
        self.assertEqual(finish_fence(tomorrow, '2026-10-29T21:00:00-07:00').hour, 21)
        with self.assertRaises(ValueError):
            finish_fence(tomorrow, '2026-10-29T21:00:00')

    def test_resource_caps_reject_unbounded_parameters(self):
        check_limits(20000, 3, 85, 180)
        for controls in [(20001, 3, 85, 180), (20000, 4, 85, 180), (20000, 3, 86, 180), (20000, 3, 85, 181)]:
            with self.assertRaises(ValueError):
                check_limits(*controls)

    def test_new_cells_exclude_baseline_and_require_every_source_corner(self):
        complete = np.ones((4, 4), dtype=bool)
        source = complete.copy()
        baseline = np.zeros_like(complete)
        baseline[1:3, 1:3] = True
        selected = new_cell_mask(complete, source, baseline, 8)
        self.assertEqual(int(selected.sum()), 8)
        self.assertFalse(selected[1, 1])
        with self.assertRaises(ValueError):
            new_cell_mask(complete, source, baseline, 7)
        source[0, 0] = False
        with self.assertRaises(ValueError):
            new_cell_mask(complete, source, baseline, 9)

    def test_recursive_donors_cannot_cross_declared_or_known_bug_floor(self):
        logt = np.array([3., 3.1, 3.2, 3.3, 3.4])
        support = np.array([[[0, 4, 1, 4]]])
        self.assertEqual(support_guard(support, logt, 3.1), 3.1)
        with self.assertRaises(ValueError):
            support_guard(support, logt, 3.2)
        with self.assertRaises(ValueError):
            support_guard(support, logt, 3.0)

    def test_polynomial_closed_cell_identity_and_nonfinite_guard(self):
        jets = np.zeros((1, 2, 2, 3, 3))
        for density_end in range(2):
            for temperature_end in range(2):
                jets[0, density_end, temperature_end, 0, 0] = 3 * density_end + density_end**2 + temperature_end**2
                jets[0, density_end, temperature_end, 1, 0] = 3 + 2 * density_end
                jets[0, density_end, temperature_end, 2, 0] = 2
                jets[0, density_end, temperature_end, 0, 1] = 2 * temperature_end
                jets[0, density_end, temperature_end, 0, 2] = 2
        powers, magnitudes = guarded_powers(jets, 1., 1., np.array([.99]), np.array([1.01]))
        np.testing.assert_allclose(polynomial_value(powers['P_rho'], [.23], [.37]), [5 + 2 * .99 * .23], atol=1e-12)
        coefficients, margin, *_ = bernstein_bounds(powers['P_rho'], magnitudes['P_rho'])
        self.assertGreater(float(coefficients.min() - margin[0]), 0)
        jets[0, 0, 0, 2, 2] = np.nan
        with self.assertRaises(ValueError):
            guarded_powers(jets, 1., 1., np.ones(1), np.ones(1))

    def test_prho_restoration_guard_and_ambiguous_polynomial_not_witness(self):
        box = tuple(map(DTYPE, (-10., -9.9, 7., 7.1)))
        maximum = radiation_cv_numerator(box[0], box[3])
        self.assertEqual(radiation_guard('P_rho', maximum, False), 0)
        guard = radiation_guard('P_rho', maximum, True)
        self.assertGreater(guard, 0)
        coefficients = np.full((6, 6), guard / 2, dtype=DTYPE)
        representation = adaptive_split(coefficients, DTYPE('1e-12'), box, 'P_rho', restored_core=False, max_depth=0, max_nodes=1)
        restored = adaptive_split(coefficients, DTYPE('1e-12'), box, 'P_rho', restored_core=True, max_depth=0, max_nodes=1)
        self.assertEqual(representation['status'], 'positive')
        self.assertEqual(restored['status'], 'unresolved')
        self.assertIsNone(restored['witness'])
        self.assertEqual(restored['visited_nodes'], 1)


if __name__ == '__main__':
    unittest.main()
