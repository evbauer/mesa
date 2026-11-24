"""Analytic identity and profile-parser controls for the stellar replay."""

from pathlib import Path
import tempfile
import unittest

import numpy as np

from audit import CRAD, FIELDS, LN10, diagnostics, profile


class IdentityTests(unittest.TestCase):
    def ideal_radiation(self, point):
        rho, temperature = np.exp(np.asarray(point) * LN10)
        gas_constant = 8.314462618e7
        gas = rho * gas_constant * temperature
        gas_energy = 1.5 * gas_constant * temperature
        radiation_energy = CRAD * temperature**4 / rho
        gas_entropy = gas_constant * (1.5 * np.log(temperature) - np.log(rho) + 20)
        radiation_entropy = 4 * radiation_energy / (3 * temperature)
        energy, entropy = gas_energy + radiation_energy, gas_entropy + radiation_entropy
        result = np.zeros(125)
        result[4:7] = np.log([gas, energy, entropy])
        result[31:34] = [1, -radiation_energy / energy, (-gas_constant - radiation_entropy) / entropy]
        result[58:61] = [1, (gas_energy + 4 * radiation_energy) / energy,
                         (1.5 * gas_constant + 3 * radiation_entropy) / entropy]
        result[4 + FIELDS["Cv"]] = (gas_energy + 4 * radiation_energy) / temperature
        return result

    def test_analytic_gas_plus_radiation(self):
        for point in ((0., 4.), (-3., 6.), (-7., 7.)):
            result = diagnostics(self.ideal_radiation(point), point)
            self.assertLess(max(abs(value) for values in result.values() for value in values.values()), 3e-14)

    def test_independent_energy_gradient_breaks_closure(self):
        point = (0., 4.)
        result = self.ideal_radiation(point)
        result[59] *= 1.01
        measured = diagnostics(result, point)["ratios"]
        self.assertGreater(abs(measured["Dse"]), .009)
        self.assertGreater(abs(measured["energy_Cv"]), .009)

    def test_profile_parser(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "profile.data"
            path.write_text("1 2\nmodel_number star_age\n8 1e6\n\n1 2 3\nzone logRho logT\n1 -3 6\n")
            metadata, data = profile(path)
            self.assertEqual(metadata["model_number"], "8")
            np.testing.assert_array_equal(data["logT"], [6.])


if __name__ == "__main__":
    unittest.main()
