"""Independent same-polynomial, matched API, and frozen-failure controls."""

from decimal import Decimal, localcontext
import json
import math
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

DIRECTORY = Path(__file__).resolve().parent
CODE = DIRECTORY.parent
WORKSPACE = DIRECTORY.parents[4]
EXPERIMENTS = WORKSPACE / "freeeos_biquintic_progress/experiments"
OUTPUT = EXPERIMENTS / "053_arithmetic_control"
sys.path.insert(0, str(CODE))
sys.path.insert(0, str(CODE / "runtime_bridge"))
from runtime_consumers import run_consumers as runner
sys.modules["run_consumers"] = runner
from runtime_consumers.diagnose_precision import ExtendedCache
import api_client
from runtime_environment import private_environment, installed_inventory, verify_private_cache, sha256


ORDERS = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2), (3, 0), (2, 1), (1, 2), (0, 3)]
COEFFICIENTS = [[1, 0, 0, -10, 15, -6], [0, 1, 0, -6, 8, -3],
                [0, 0, .5, -1.5, 1.5, -.5], [0, 0, 0, 10, -15, 6],
                [0, 0, 0, -4, 7, -3], [0, 0, 0, .5, -1, .5]]


def decimal_weights(coordinate, lower, upper, order):
    spacing = Decimal.from_float(upper) - Decimal.from_float(lower)
    fraction = (Decimal.from_float(coordinate) - Decimal.from_float(lower)) / spacing
    weights = []
    for basis, coefficients in enumerate(COEFFICIENTS):
        value = Decimal(0)
        for degree in range(order, 6):
            multiplier = math.factorial(degree) // math.factorial(degree - order)
            power = Decimal(1) if degree == order else fraction**(degree - order)
            value += Decimal(coefficients[degree]) * multiplier * power
        weights.append(value * spacing**(basis % 3 - order))
    return weights


def reference_native(bounds, coordinate, corners):
    with localcontext() as context:
        context.prec = 80
        density = {order: decimal_weights(coordinate[0], *bounds[:2], order) for order in range(4)}
        thermal = {order: decimal_weights(coordinate[1], *bounds[2:], order) for order in range(4)}
        values = []
        for density_order, thermal_order in ORDERS:
            value = Decimal(0)
            for endpoint_t in range(2):
                for endpoint_u in range(2):
                    for derivative_t in range(3):
                        for derivative_u in range(3):
                            datum = Decimal.from_float(float(corners[endpoint_t, endpoint_u, derivative_u + 3 * derivative_t]))
                            value += datum * density[density_order][3 * endpoint_u + derivative_u] * thermal[thermal_order][3 * endpoint_t + derivative_t]
            values.append(float(value))
    return values


def native_checks(output):
    cache = ExtendedCache()
    cases = json.loads((EXPERIMENTS / "038_solver_manifest/queries_v2.json").read_text())["cases"]
    samples = []
    for case in cases:
        for offset in (-1e-8, -.001, .001, 1e-8):
            native = np.array([case["logQ"], case["logT"]]) * runner.LN10
            native[0 if case["direction"] == "u" else 1] += offset
            samples.append((case["case_id"], native))
    valid = (cache.rows[:-1, :-1, 2] == 1) & (cache.rows[1:, :-1, 2] == 1) & (cache.rows[:-1, 1:, 2] == 1) & (cache.rows[1:, 1:, 2] == 1)
    eligible = np.argwhere(valid & (cache.t[:-1, None] >= 3 * runner.LN10))
    generator = np.random.default_rng(53029)
    for chosen in generator.choice(len(eligible), 64, replace=False):
        row, column = eligible[chosen]
        fraction = generator.uniform(.05, .95, 2)
        samples.append(("random_clean_cell", np.array([cache.u[column], cache.t[row]]) + fraction *
                        np.array([cache.u[column + 1] - cache.u[column], cache.t[row + 1] - cache.t[row]])))
    records, lines = [], []
    for label, coordinate in samples:
        column = int(np.searchsorted(cache.u, coordinate[0], side="right") - 1)
        row = int(np.searchsorted(cache.t, coordinate[1], side="right") - 1)
        corners = cache.rows[row:row + 2, column:column + 2, 3:]
        if not np.all(cache.rows[row:row + 2, column:column + 2, 2] == 1):
            raise AssertionError("Native check requires a complete clean cell")
        bounds = [float(cache.u[column]), float(cache.u[column + 1]), float(cache.t[row]), float(cache.t[row + 1])]
        lines.extend([" ".join(format(value, ".17g") for value in bounds), " ".join(format(value, ".17g") for value in coordinate)])
        for corner in corners.reshape(4, 9):
            lines.append(" ".join(format(value, ".17g") for value in corner))
        records.append(dict(label=label, coordinate=coordinate.tolist(), cell=[column, row], bounds=bounds,
                            reference=reference_native(bounds, coordinate, corners)))
    inputs = "\n".join(lines) + "\n"
    (output / "native_queries.txt").write_text(inputs)
    summary = {}
    for name, executable in (("baseline", "baseline_cell_driver"), ("variant", "cell_driver")):
        binary = OUTPUT / "build/driver" / executable
        response = subprocess.run([str(binary)], input=inputs, text=True, capture_output=True, check=True,
                                  env=private_environment(output / (name + "_native")))
        (output / (name + "_native.stdout")).write_text(response.stdout)
        (output / (name + "_native.stderr")).write_text(response.stderr)
        data = response.stdout.splitlines()
        if len(data) != 3 * len(records):
            raise AssertionError("Native protocol length differs")
        errors, absolute = [], []
        for index, record in enumerate(records):
            value = np.fromstring(data[3 * index + 1], sep=" ")
            reference = np.array(record["reference"])
            difference = value - reference
            record[name] = value.tolist()
            record[name + "_scaled_error"] = (difference / np.maximum(1, np.abs(reference))).tolist()
            errors.append(np.abs(difference) / np.maximum(1, np.abs(reference)))
            absolute.append(np.abs(difference))
        summary[name] = dict(maximum_scaled_by_partial=np.max(errors, axis=0).tolist(),
                             maximum_absolute_by_partial=np.max(absolute, axis=0).tolist(),
                             maximum_scaled=float(np.max(errors)), binary_sha256=sha256(binary))
    reference = np.array([record["reference"] for record in records])
    variant = np.array([record["variant"] for record in records])
    baseline = np.array([record["baseline"] for record in records])
    for name, values in (("baseline", baseline), ("variant", variant)):
        summary[name]["rms_scaled_by_partial"] = np.sqrt(np.mean(((values - reference) / np.maximum(1, np.abs(reference)))**2, axis=0)).tolist()
    summary.update(samples=len(records), decimal_digits=80, derivative_orders=ORDERS,
                   nodes="exact binary64 cached corner data; no source re-estimation", forced_cells="native isolated driver only, not EOS API")
    runner.write_json(output / "native_records.json", records)
    runner.write_json(output / "native_results.json", summary)
    return summary


def api_checks(output):
    baseline_experiment = api_client.EXPERIMENT
    frozen_samples = [json.loads(line) for line in (EXPERIMENTS / "051_runtime_consumers/precision_control/representable_states.jsonl").read_text().splitlines()]
    failures = json.loads((EXPERIMENTS / "051_runtime_consumers/precision_control/frozen_failures.json").read_text())
    failure = next(record for record in failures if record["route"] == "component")
    cases = json.loads((EXPERIMENTS / "038_solver_manifest/queries_v2.json").read_text())["cases"]
    case = next(record for record in cases if record["case_id"] == failure["case_id"])
    density = case["logRho"] * runner.LN10
    api_client.EXPERIMENT = OUTPUT
    records = []
    with api_client.EOSClient(output / "variant_api") as variant:
        for sample in frozen_samples:
            point = sample["y"]
            values = variant.query("candidate", sample["route"], density / runner.LN10, point / runner.LN10)
            if values[1] != 0:
                raise AssertionError("Failure-neighborhood query rejected")
            residual = point - failure["initial_y"] - failure["H_erg_per_g_s"] * failure["dt_seconds"] / (np.exp(point) * values[14])
            records.append(dict(route=sample["route"], y=point, baseline_residual=sample["residual"],
                                variant_residual=float(residual), Cv=float(values[14]), Cv_t=float(values[68]),
                                callback_count=int(values[2]), fractions=values[23:31].tolist()))
        api_client.EXPERIMENT = baseline_experiment
        matched = []
        with api_client.EOSClient(output / "baseline_api") as baseline:
            states = [(case["case_id"], case["logRho"], case["logT"] + offset) for case in cases for offset in (-.0001, .0001)]
            states.extend([("inactive", 4., 7.), ("blend", -2., 5.)])
            for label, log_rho, log_temp in states:
                for route in ("component", "final"):
                    for scenario in ("candidate", "off", "composition_fault", "input_fault"):
                        first = baseline.query(scenario, route, log_rho, log_temp)
                        second = variant.query(scenario, route, log_rho, log_temp)
                        if not np.array_equal(first[1:4], second[1:4]):
                            raise AssertionError("Arithmetic variant changed status/callback/retention")
                        if first[1] == 0:
                            difference = (second[4:] - first[4:]) / np.maximum(1, np.abs(first[4:]))
                            if not np.array_equal(first[85:], second[85:]):
                                raise AssertionError("Hybrid abundance changed")
                            if scenario == "off" and not np.array_equal(first[4:], second[4:]):
                                raise AssertionError("Opt-out numerical response changed")
                        else:
                            difference = np.zeros(121)
                        matched.append(dict(label=label, scenario=scenario, route=route, logRho=log_rho, logT=log_temp,
                                            status=int(first[1]), callbacks=int(first[2]),
                                            scaled_changes=difference.tolist(),
                                            values_baseline=first[4:].tolist() if first[1] == 0 else None,
                                            values_variant=second[4:].tolist() if second[1] == 0 else None))
    api_client.EXPERIMENT = baseline_experiment
    summary = {}
    for route in ("component", "final"):
        points = [record for record in records if record["route"] == route]
        values = np.array([record["variant_residual"] for record in points])
        summary[route] = dict(samples=len(points), max_adjacent_jump=float(max(abs(np.diff(values)))),
                              decreasing_steps=int(np.count_nonzero(np.diff(values) < 0)),
                              minimum_absolute=float(min(abs(values))),
                              passing_frozen_tolerance=int(np.count_nonzero(abs(values) <= 1e-13)))
    runner.write_json(output / "neighborhood_records.json", records)
    runner.write_json(output / "neighborhood_results.json", summary)
    runner.write_json(output / "matched_api_records.json", matched)
    runner.write_json(output / "matched_api_results.json", dict(queries=len(matched),
                      nonzero_statuses=sum(record["status"] != 0 for record in matched),
                      maximum_scaled_by_column=np.max([np.abs(record["scaled_changes"]) for record in matched], axis=0).tolist(),
                      zero_status_callback_hybrid_and_optout_violations=True))
    for label in ("variant_api", "baseline_api"):
        verify_private_cache(output / label)
    return summary


def main():
    output = OUTPUT / "validation"
    if output.exists():
        raise RuntimeError("Preserve existing validation output")
    output.mkdir()
    started = time.time()
    before = installed_inventory()
    native = native_checks(output)
    api = api_checks(output)
    after = installed_inventory()
    if before != after:
        raise AssertionError("Installed EOS inventory changed")
    runner.write_json(output / "results.json", dict(native=native, api=api, installed_unchanged=True,
                      elapsed_seconds=time.time() - started, source_sha256=sha256(Path(__file__))))
    print(json.dumps(dict(native=native, neighborhood=api, elapsed_seconds=time.time() - started), indent=2))


if __name__ == "__main__":
    main()
