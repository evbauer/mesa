"""Independent extended-precision check of cancellation-sensitive dE/dRho."""

import json
import sys

import numpy as np

from api_client import CACHE, EOSClient, LN10, CRAD
from runtime_environment import CODE, WORKSPACE, EXPERIMENT

sys.path.insert(0, str(CODE.parent))
from positivity_bounds.bounds import cell_power, derivative, polynomial_value


REAL = np.longdouble


class ExtendedCache:
    def __init__(self):
        with CACHE.open() as stream:
            header = [stream.readline() for _ in range(6)]
        self.nq, self.nt = map(int, header[1].split())
        minimum_u, spacing_u, minimum_t, spacing_t = map(float, header[5].split())
        self.u_axis = minimum_u + np.arange(self.nq) * spacing_u
        self.t_axis = minimum_t + np.arange(self.nt) * spacing_t
        records = np.loadtxt(CACHE, skiprows=6).reshape(self.nt, self.nq, 12).swapaxes(0, 1)
        self.valid = records[:, :, 2] == 1
        self.jets = records[:, :, 3:].reshape(self.nq, self.nt, 3, 3).swapaxes(2, 3)
        self.powers = {}

    def cell(self, point):
        density_log, temperature_log = point
        native_u = density_log - 2 * temperature_log + REAL(12 * LN10)
        density_index = int(np.searchsorted(self.u_axis, native_u, side="right") - 1)
        temperature_index = int(np.searchsorted(self.t_axis, temperature_log, side="right") - 1)
        if not 0 <= density_index < self.nq - 1 or not 0 <= temperature_index < self.nt - 1:
            raise ValueError("Outside native cache")
        if not np.all(self.valid[density_index:density_index + 2, temperature_index:temperature_index + 2]):
            raise ValueError("Incomplete cache cell")
        return density_index, temperature_index

    def evaluate(self, point):
        point = np.asarray(point, dtype=REAL)
        density_log, temperature_log = point
        native_u = density_log - 2 * temperature_log + REAL(12 * LN10)
        density_index, temperature_index = self.cell(point)
        spacing_u = REAL(self.u_axis[density_index + 1] - self.u_axis[density_index])
        spacing_t = REAL(self.t_axis[temperature_index + 1] - self.t_axis[temperature_index])
        key = (density_index, temperature_index)
        if key not in self.powers:
            jets = self.jets[density_index:density_index + 2, temperature_index:temperature_index + 2]
            self.powers[key] = cell_power(jets, spacing_u, spacing_t)[0]
        phase_u = (native_u - REAL(self.u_axis[density_index])) / spacing_u
        phase_t = (temperature_log - REAL(self.t_axis[temperature_index])) / spacing_t
        native = {(order_u, order_t): polynomial_value(
            derivative(self.powers[key], order_u, order_t, spacing_u, spacing_t), phase_u, phase_t)
            for order_u, order_t in ((1, 0), (2, 0), (1, 1), (3, 0), (2, 1), (1, 2))}
        radiation_f = -REAL(CRAD) * np.exp(4 * temperature_log - density_log) / 3
        first = native[1, 0] - radiation_f
        second = native[2, 0] + radiation_f
        mixed = native[1, 1] - 2 * native[2, 0] - 4 * radiation_f
        rho_mixed = native[2, 1] - 2 * native[3, 0] + 4 * radiation_f
        temp_mixed = native[1, 2] - 4 * native[2, 1] + 4 * native[3, 0] - 16 * radiation_f
        rho = np.exp(density_log)
        return np.array([(first - mixed) / rho, (second - rho_mixed - first + mixed) / rho,
                         (mixed - temp_mixed) / rho], dtype=REAL)


def run():
    cache = ExtendedCache()
    manifest = WORKSPACE / "freeeos_biquintic_progress/experiments/038_solver_manifest/queries_v2.json"
    cases = json.loads(manifest.read_text())["cases"]
    records = []
    with EOSClient(EXPERIMENT / "precision_check") as client:
        for case in cases:
            for side in (-1, 1):
                point = np.array([case["logRho"], case["logT"]]) * LN10
                point[0 if case["direction"] == "u" else 1] += side * .001
                extended = cache.evaluate(point)
                result = client.query("candidate", "component", *(point / LN10))
                if result[1] != 0:
                    raise ValueError("Actual component failure")
                returned = result[[4 + 11, 31 + 11, 58 + 11]]
                scale = np.maximum.reduce([abs(extended), np.full(3, abs(extended[0])), np.ones(3)])
                difference = abs(extended - returned) / scale
                fd_checks = []
                for axis in range(2):
                    errors = []
                    for step_text in ("0.0002", "0.0001", "0.00003", "0.00001"):
                        step = REAL(step_text)
                        offset = np.zeros(2, dtype=REAL)
                        offset[axis] = step
                        if any(cache.cell(point + factor * offset) != cache.cell(point) for factor in (-2, -1, 1, 2)):
                            raise ValueError("Precision finite difference leaves native cell")
                        sampled = {factor: cache.evaluate(point + factor * offset)[0] for factor in (-2, -1, 1, 2)}
                        finite = (8 * (sampled[1] - sampled[-1]) - (sampled[2] - sampled[-2])) / (12 * step)
                        errors.append(float(abs(finite - extended[axis + 1]) / scale[axis + 1]))
                    fd_checks.append(errors)
                records.append({"case_id": case["case_id"], "side": side,
                    "returned": returned.tolist(), "extended": [float(value) for value in extended],
                    "scaled_difference": [float(value) for value in difference],
                    "extended_finite_difference_error": fd_checks})
    maximum = max(max(record["scaled_difference"]) for record in records)
    finite_maximum = max(max(min(errors) for errors in record["extended_finite_difference_error"]) for record in records)
    summary = {"scope": "same archived double-precision nodal data; independent long-double power assembly and physical shear",
        "longdouble_mantissa_bits": np.finfo(REAL).nmant, "queries": len(records),
        "maximum_scaled_difference": maximum, "maximum_extended_fd_error": finite_maximum, "records": records}
    (EXPERIMENT / "precision_check/summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key != "records"}, indent=2))
    if maximum > 3e-6 or finite_maximum > 1e-7:
        raise AssertionError("Independent precision check exceeds budget")


if __name__ == "__main__":
    run()
