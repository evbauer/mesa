"""Leakage, shared-sample, local NAK, and smooth full-pipeline controls."""

import unittest

import numpy as np
from scipy.interpolate import CubicSpline

from refinement import (
    CRAD, STRIDES, build_grids, common_queries, evaluate_candidate,
    lagrange_four, own_use, paired_ratios, table_reference,
)


def analytic_table():
    logq, logt = np.linspace(-4, 1, 257), np.linspace(3.5, 5, 193)
    temperature = np.broadcast_to(10.0**logt, (len(logq), len(logt)))
    rho = 10.0**(logq[:, None] + 2 * logt - 12)
    gas_constant = 1e8
    pgas = rho * gas_constant * temperature
    prad = CRAD * temperature**4 / 3
    total = pgas + prad
    cv = 1.5 * gas_constant + 4 * CRAD * temperature**3 / rho
    chit, chirho = (pgas + 4 * prad) / total, pgas / total
    rows = np.zeros((*rho.shape, 22))
    rows[..., 0] = logt
    rows[..., 1] = np.log10(pgas)
    rows[..., 2] = np.log10(1.5 * gas_constant * temperature + CRAD * temperature**4 / rho)
    rows[..., 3] = np.log10(gas_constant * (1.5 * np.log(temperature) - np.log(rho))
                            + 4 * CRAD * temperature**3 / (3 * rho))
    rows[..., 4], rows[..., 5] = chirho, chit
    rows[..., 6] = cv + total * chit**2 / (rho * temperature * chirho)
    rows[..., 7] = cv
    rows[..., 8] = -CRAD * temperature**4 / rho**2
    rows[..., 9] = cv / temperature
    rows[..., 10] = -gas_constant / rho - 4 * CRAD * temperature**3 / (3 * rho**2)
    rows[..., 18] = np.log10(rho)
    return rows, logq, logt


class RefinementTests(unittest.TestCase):
    def test_four_node_notaknot_equivalence(self):
        generator = np.random.default_rng(14)
        values = generator.normal(size=4)
        positions = np.linspace(0, 1, 31)
        expected = CubicSpline([-1, 0, 1, 2], values, bc_type='not-a-knot')(positions)
        np.testing.assert_allclose(lagrange_four(positions) @ values, expected, atol=1e-14)

    def test_common_holdouts_and_smooth_pipeline(self):
        rows, logq, logt = analytic_table()
        grids = build_grids(rows, logq, logt)
        density, temperature, common, _ = common_queries(rows, logq, logt, grids)
        selected = np.flatnonzero(common)[::max(1, int(common.sum()) // 120)]
        self.assertGreater(len(selected), 10)
        self.assertTrue(np.all(density[selected] % 2 == 1))
        self.assertTrue(np.all(temperature[selected] % 2 == 1))
        reference = table_reference(rows[density[selected], temperature[selected]])
        errors = []
        for stride in STRIDES:
            lower_density = density[selected] // stride
            lower_temperature = temperature[selected] // stride
            for density_offset in [0, 1]:
                for temperature_offset in [0, 1]:
                    self.assertTrue(np.all(own_use(
                        grids[stride]['q'][lower_density + density_offset],
                        grids[stride]['t'][lower_temperature + temperature_offset])))
            actual = evaluate_candidate('BQ_split', grids[stride], logq[density[selected]], logt[temperature[selected]])
            errors.append(max(abs((actual['Cv'] - reference['Cv']) / reference['Cv'])))
        self.assertLess(errors[0], errors[1])
        self.assertLess(errors[1], errors[2])

    def test_no_withheld_constraint_leak(self):
        rows, logq, logt = analytic_table()
        original = build_grids(rows, logq, logt)
        changed = rows.copy()
        changed[1::2, 1::2, 2] += .5
        changed[1::2, 1::2, 7] *= 100
        perturbed = build_grids(changed, logq, logt)
        for stride in STRIDES:
            for name in ['BQ_total', 'BQ_split']:
                np.testing.assert_array_equal(original[stride]['candidates'][name]['nodes'],
                                              perturbed[stride]['candidates'][name]['nodes'])
            np.testing.assert_array_equal(original[stride]['rows'], perturbed[stride]['rows'])

    def test_ratio_zero_cases(self):
        ratios, summary = paired_ratios(np.array([1., 0., 2., 1.]), np.array([0., 0., 1., 2.]))
        self.assertTrue(np.isinf(ratios[0]))
        self.assertTrue(np.isnan(ratios[1]))
        self.assertEqual(summary['fine_zero_coarse_positive'], 1)
        self.assertEqual(summary['both_zero_undefined'], 1)
        self.assertEqual(summary['fraction_improved_with_refinement'], .5)


if __name__ == '__main__':
    unittest.main()
