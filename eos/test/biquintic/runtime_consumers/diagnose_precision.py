"""Separate frozen-failure precision/tolerance controls; does not replace full results."""

import json
import math
from pathlib import Path
import resource
import time

import numpy as np

from run_consumers import CODE, EXPERIMENTS, OUTPUT, Geometry, LoggedClient, LN10, scalar_solver, write_json
from api_client import CACHE, CRAD


class ExtendedCache:
    def __init__(self, precision=np.longdouble):
        self.real = precision
        with CACHE.open() as stream:
            header = [next(stream).split() for _ in range(6)]
        nq, nt = map(int, header[1])
        umin, du, tmin, dt = map(float, header[5])
        self.u = umin + du * np.arange(nq)
        self.t = tmin + dt * np.arange(nt)
        self.rows = np.loadtxt(CACHE, skiprows=6).reshape(nt, nq, 12)
        self.coefficients = np.array([[1, 0, 0, -10, 15, -6], [0, 1, 0, -6, 8, -3],
                                      [0, 0, .5, -1.5, 1.5, -.5], [0, 0, 0, 10, -15, 6],
                                      [0, 0, 0, -4, 7, -3], [0, 0, 0, .5, -1, .5]], dtype=self.real)

    def basis(self, coordinate, lower, upper, order):
        spacing = self.real(upper) - self.real(lower)
        position = (self.real(coordinate) - self.real(lower)) / spacing
        weights = np.empty((2, 3), dtype=self.real)
        for endpoint in range(2):
            for derivative in range(3):
                coefficients = self.coefficients[3 * endpoint + derivative]
                polynomial = self.real(0)
                for degree in range(5, order - 1, -1):
                    multiplier = math.factorial(degree) // math.factorial(degree - order)
                    polynomial = polynomial * position + coefficients[degree] * multiplier
                weights[endpoint, derivative] = polynomial * spacing**(derivative - order)
        return weights

    def evaluate(self, density_log, temperature_log):
        density = 10.0**float(density_log / LN10)
        temperature = 10.0**float(temperature_log / LN10)
        physical_x, physical_t = math.log(density), math.log(temperature)
        native_u = physical_x - 2 * physical_t + 12 * LN10
        column = int(np.searchsorted(self.u, native_u, side="right") - 1)
        row = int(np.searchsorted(self.t, physical_t, side="right") - 1)
        corners = self.rows[row:row + 2, column:column + 2]
        if not np.all(corners[:, :, 2] == 1):
            raise ValueError("Incomplete extended-control cell")
        density_weights = {order: self.basis(native_u, *self.u[column:column + 2], order) for order in range(4)}
        thermal_weights = {order: self.basis(physical_t, *self.t[row:row + 2], order) for order in range(4)}
        native = {}
        for density_order in range(4):
            for thermal_order in range(4 - density_order):
                value = self.real(0)
                for thermal_end in range(2):
                    for density_end in range(2):
                        for thermal_node in range(3):
                            for density_node in range(3):
                                value += (self.real(corners[thermal_end, density_end, 3 + density_node + 3 * thermal_node]) *
                                          density_weights[density_order][density_end, density_node] *
                                          thermal_weights[thermal_order][thermal_end, thermal_node])
                native[density_order, thermal_order] = value
        first = native[0, 1] - 2 * native[1, 0]
        second = native[0, 2] - 4 * native[1, 1] + 4 * native[2, 0]
        third = native[0, 3] - 6 * native[1, 2] + 12 * native[2, 1] - 8 * native[3, 0]
        radiation_cv = 4 * self.real(CRAD) * np.exp(3 * self.real(physical_t) - self.real(physical_x))
        cv = (first - second) / np.exp(self.real(physical_t)) + radiation_cv
        cv_t = (2 * second - third - first) / np.exp(self.real(physical_t)) + 3 * radiation_cv
        return cv, cv_t


def main():
    output = OUTPUT / "precision_control"
    output.mkdir(parents=True, exist_ok=True)
    complete = [json.loads(line) for line in (OUTPUT / "full/trials.jsonl").read_text().splitlines()]
    failures = [record for record in complete if record["status"] != "converged"]
    manifest = json.loads((EXPERIMENTS / "038_solver_manifest/queries_v2.json").read_text())
    cases = {case["case_id"]: case for case in manifest["cases"]}
    write_json(output / "frozen_failures.json", failures)
    extended = ExtendedCache()
    geometry = Geometry()
    summaries, samples, controls, derivatives = [], [], [], []
    started = time.time()
    with LoggedClient(output / "api") as client:
        for failure in failures:
            case = cases[failure["case_id"]]
            density_log = case["logRho"] * LN10
            center = failure["newton"]["history"][-1]["y"]
            initial, heating = failure["initial_y"], failure["H_erg_per_g_s"] * failure["dt_seconds"]
            points = [center]
            lower = upper = center
            for position in range(256):
                lower = float(np.nextafter(lower, -np.inf))
                upper = float(np.nextafter(upper, np.inf))
                points.extend([lower, upper])
            points = sorted(points)
            client.begin(dict(consumer="precision_control", route=failure["route"], model="candidate",
                              case_id=case["case_id"], trial_id=failure["trial_id"]))
            client.stage = "adjacent_representable_state_scan"
            local = []
            for point in points:
                actual = client.query("candidate", failure["route"], density_log / LN10, point / LN10)
                direct = client.query("direct", "component", density_log / LN10, point / LN10)
                if actual[1] or direct[1]:
                    raise ValueError("Unexpected API error in preserved-failure neighborhood")
                cv, cv_t = actual[14], actual[68]
                residual = point - initial - heating / (np.exp(point) * cv)
                wide_residual = np.longdouble(point) - np.longdouble(initial) - np.longdouble(heating) / (
                    np.exp(np.longdouble(point)) * np.longdouble(cv))
                extended_cv, extended_cv_t = extended.evaluate(density_log, point)
                extended_residual = np.longdouble(point) - np.longdouble(initial) - np.longdouble(heating) / (
                    np.exp(np.longdouble(point)) * extended_cv)
                record = dict(route=failure["route"], trial_id=failure["trial_id"], y=point,
                              decimal_logT=float(point / LN10), effective_y=math.log(10.0**float(point / LN10)),
                              Cv=float(cv), Cv_t=float(cv_t), residual=float(residual),
                              residual_extended_consumer_only=float(wide_residual),
                              Cv_extended_cell=float(extended_cv), Cv_t_extended_cell=float(extended_cv_t),
                              residual_extended_cell=float(extended_residual),
                              component_direct_bitwise_equal=np.array_equal(actual[4:85], direct[4:85]),
                              fractions=actual[23:31].tolist(), native_cell=geometry.cell(density_log, point))
                local.append(record)
                samples.append(record)
            residuals = np.array([record["residual"] for record in local])
            wide = np.array([record["residual_extended_consumer_only"] for record in local])
            corrected = np.array([record["residual_extended_cell"] for record in local])
            effective_steps = np.diff([record["effective_y"] for record in local])
            cv_changes = np.diff([record["Cv"] for record in local])
            cv_slopes = np.array([record["Cv_t"] for record in local])
            jitter = cv_changes - .5 * (cv_slopes[:-1] + cv_slopes[1:]) * effective_steps
            normal_coordinates = [density_log - 2 * center + 12 * LN10, center]
            distance = min(float(np.min(abs(knots - coordinate))) for knots, coordinate in zip((geometry.u, geometry.t), normal_coordinates))
            summaries.append(dict(route=failure["route"], trial_id=failure["trial_id"], case_id=case["case_id"],
                                  tolerance=failure["tolerance_absolute_R"], samples=len(local),
                                  cell_count=len({tuple(record["native_cell"]) for record in local}),
                                  nearest_native_face_distance=distance,
                                  min_abs_residual=float(min(abs(residuals))),
                                  samples_meeting_frozen_tolerance=int(np.count_nonzero(abs(residuals) <= failure["tolerance_absolute_R"])),
                                  min_abs_extended_consumer_residual=float(min(abs(wide))),
                                  consumer_arithmetic_max_change=float(max(abs(wide - residuals))),
                                  min_abs_extended_cell_residual=float(min(abs(corrected))),
                                  max_adjacent_residual_jump=float(max(abs(np.diff(residuals)))),
                                  decreasing_residual_steps=int(np.count_nonzero(np.diff(residuals) < 0)),
                                  max_Cv_adjacent_change_minus_analytic_slope=float(max(abs(jitter))),
                                  pure_FreeEOS=all(record["fractions"] == [0, 0, 1, 0, 0, 0, 0, 0] for record in local),
                                  direct_bitwise_parity=all(record["component_direct_bitwise_equal"] for record in local)))

            def api_callback(point):
                state = client.state("candidate", failure["route"], density_log, point)
                return state["Cv"], state["Cv_t"]

            client.stage = "residual_jacobian_step_sweep"
            analytic = scalar_solver.residual(api_callback, center, initial, heating)[1]
            for step in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8):
                if geometry.cell(density_log, center - step) != geometry.cell(density_log, center + step):
                    raise ValueError("Failure-control derivative crosses cell")
                finite = (scalar_solver.residual(api_callback, center + step, initial, heating)[0] -
                          scalar_solver.residual(api_callback, center - step, initial, heating)[0]) / (2 * step)
                derivatives.append(dict(route=failure["route"], step=step, analytic=float(analytic),
                                        finite_difference=float(finite), scaled_error=float((finite - analytic) / max(1, abs(analytic)))))
            client.stage = "separately_named_tolerance_control"
            for tolerance in (1e-13, 3e-13, 1e-12, 1e-11):
                outcome = scalar_solver.solve(api_callback, initial, heating, *failure["bracket"],
                                              failure["solver_start_y"], tolerance=tolerance, max_iterations=40)
                controls.append(dict(route=failure["route"], method="actual_API", tolerance=tolerance,
                                     status=outcome.status, iterations=outcome.iterations,
                                     final_residual=outcome.history[-1]["residual"], history=outcome.history,
                                     rejected_history=outcome.rejected_history))
            def extended_callback(point):
                return tuple(float(value) for value in extended.evaluate(density_log, point))
            outcome = scalar_solver.solve(extended_callback, initial, heating, *failure["bracket"],
                                          failure["solver_start_y"], tolerance=1e-13, max_iterations=40)
            controls.append(dict(route=failure["route"], method="offline_extended_cell", tolerance=1e-13,
                                 status=outcome.status, iterations=outcome.iterations,
                                 final_residual=outcome.history[-1]["residual"], history=outcome.history,
                                 rejected_history=outcome.rejected_history,
                                 scope="Same binary64 cache nodes and effective API coordinates; longdouble polynomial arithmetic; diagnostic only"))
    with (output / "representable_states.jsonl").open("w") as stream:
        for sample in samples:
            stream.write(json.dumps(sample) + "\n")
    write_json(output / "tolerance_controls.json", controls)
    write_json(output / "derivative_checks.json", derivatives)
    result = dict(failures=summaries, float_significand_bits=np.finfo(float).nmant + 1,
                  longdouble_significand_bits=np.finfo(np.longdouble).nmant + 1,
                  elapsed_seconds=time.time() - started, maxrss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  controls=[{key: record[key] for key in ("route", "method", "tolerance", "status", "iterations", "final_residual")}
                            for record in controls], source_or_core_changes=False,
                  scope="Separate diagnosis; frozen full-run failures remain failures")
    write_json(output / "results.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
