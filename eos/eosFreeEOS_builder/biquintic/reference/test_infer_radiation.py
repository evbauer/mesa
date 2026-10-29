"""Small analytic tests of table-audit algebra, not source-EOS validation."""

import unittest

import numpy as np

from infer_radiation import CRAD, clean_mask, infer_fields, print_half_unit


def analytic_rows(source_a, strip_a):
    rho = np.array([1e-10, 1e-9, 1e-8])
    temperature = np.array([1e6, 2e6, 3e6])
    gas_constant = 1e8
    material_pressure = rho * gas_constant * temperature
    prad = source_a * temperature**4 / 3
    pressure = material_pressure + prad
    pgas_stored = pressure - strip_a * temperature**4 / 3
    chit = (material_pressure + 4 * prad) / pressure
    chirho = material_pressure / pressure
    cv = 1.5 * gas_constant + 4 * source_a * temperature**3 / rho
    rows = np.ones((3, 22))
    rows[:, 0] = np.log10(temperature)
    rows[:, 1] = np.log10(pgas_stored)
    rows[:, 2] = np.log10(1.5 * gas_constant * temperature + source_a * temperature**4 / rho)
    rows[:, 3] = np.log10(gas_constant * (1.5 * np.log(temperature) - np.log(rho))
                         + 4 * source_a * temperature**3 / (3 * rho))
    rows[:, 4], rows[:, 5] = chirho, chit
    rows[:, 6] = cv + pressure * chit**2 / (rho * temperature * chirho)
    rows[:, 7] = cv
    rows[:, 8] = -source_a * temperature**4 / rho**2
    rows[:, 9] = cv / temperature
    rows[:, 10] = -gas_constant / rho - 4 * source_a * temperature**3 / (3 * rho**2)
    rows[:, 17] = 0
    rows[:, 18] = np.log10(rho)
    return rows


class RadiationAuditTests(unittest.TestCase):
    def test_print_units(self):
        self.assertEqual(print_half_unit('-9.00999999998238E-02'), 5e-17)
        self.assertEqual(print_half_unit('2.10000000000009E+00'), 5e-15)

    def test_both_routes_and_cp(self):
        rows = analytic_rows(CRAD, CRAD)
        fields = infer_fields(rows, np.zeros_like(rows))
        np.testing.assert_allclose(fields['a_maxwell'], CRAD, rtol=2e-14)
        np.testing.assert_allclose(fields['a_energy'], CRAD, rtol=2e-14)
        np.testing.assert_allclose(fields['cp_identity_current'], 0, atol=2e-14)

    def test_strip_is_not_internal_source_constant(self):
        source_a = CRAD * 1.000001
        rows = analytic_rows(source_a, CRAD)
        fields = infer_fields(rows, np.zeros_like(rows))
        np.testing.assert_allclose(fields['a_maxwell'], CRAD, rtol=2e-14)
        np.testing.assert_allclose(fields['a_energy'], CRAD, rtol=2e-14)
        np.testing.assert_allclose(fields['a_erho_proxy_over_crad_minus1'], 1e-6, atol=2e-14)
        np.testing.assert_allclose(fields['a_linear_density_pressure_proxy_over_crad_minus1'], 1e-6, atol=2e-14)
        np.testing.assert_allclose(fields['a_erho_minus_linear_density_pressure_proxy_over_crad'], 0, atol=2e-14)

    def test_print_rounding_encloses_known_constant(self):
        rows = analytic_rows(CRAD, CRAD)
        tokens = [[f'{value:.14E}' for value in record] for record in rows]
        rounded = np.asarray(tokens, dtype=float)
        bounds = np.array([[print_half_unit(token) for token in record] for record in tokens])
        fields = infer_fields(rounded, bounds)
        self.assertTrue(np.all(abs(fields['a_maxwell'] - CRAD)
                               <= CRAD * fields['a_maxwell_rounding_bound_over_crad']))
        self.assertTrue(np.all(abs(fields['a_energy'] - CRAD)
                               <= CRAD * fields['a_energy_rounding_bound_over_crad']))

    def test_exclusion_and_clean_stencil(self):
        rows = np.ones((527, 7, 22))
        header = np.array([51, .7, .02, 7, 3.2, 3.32, .02, 527, -10.09, 5.69, .03])
        rows[..., 0] = 3.2 + np.arange(7)[None, :] * .02
        rows[..., 17] = 0
        rows[..., 18] = -10.09 + np.arange(527)[:, None] * .03 + 2 * rows[..., 0] - 12
        valid, clean, _ = clean_mask(rows, header, 2)
        self.assertFalse(valid[520].any())
        self.assertFalse(clean[518:523].any())
        self.assertTrue(clean[517, 3])


if __name__ == '__main__':
    unittest.main()
