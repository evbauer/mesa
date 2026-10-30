"""Independent third-order seed composition and bounded actual-API parity."""

from collections import Counter
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
OUTPUT = EXPERIMENTS / "059_arithmetic_ad"
BASELINE = EXPERIMENTS / "053_arithmetic_control"
sys.path.insert(0, str(CODE))
from arithmetic_control import validate as prior
from runtime_adapter import run_tests
import api_client
from runtime_environment import private_environment, installed_inventory, sha256

ORDERS = prior.ORDERS


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def add(first, second):
    result = first.copy()
    for order, value in second.items():
        result[order] = result.get(order, Decimal(0)) + value
    return result


def multiply(first, second):
    result = {}
    for first_order, first_value in first.items():
        for second_order, second_value in second.items():
            order = (first_order[0] + second_order[0], first_order[1] + second_order[1])
            if sum(order) <= 3:
                result[order] = result.get(order, Decimal(0)) + first_value * second_value
    return result


def scale(polynomial, value):
    return {order: coefficient * value for order, coefficient in polynomial.items()}


def series(seed):
    return {order: Decimal.from_float(float(value)) / (math.factorial(order[0]) * math.factorial(order[1]))
            for order, value in zip(ORDERS, seed)}


def derivative_vector(polynomial):
    return [float(polynomial.get(order, Decimal(0)) * math.factorial(order[0]) * math.factorial(order[1]))
            for order in ORDERS]


def basis_series(seed, lower, upper):
    spacing = Decimal.from_float(upper) - Decimal.from_float(lower)
    position = scale(add(series(seed), {(0, 0): -Decimal.from_float(lower)}), 1 / spacing)
    weights = []
    for basis, coefficients in enumerate(prior.COEFFICIENTS):
        polynomial = {}
        for coefficient in reversed(coefficients):
            polynomial = add(multiply(polynomial, position), {(0, 0): Decimal(coefficient)})
        weights.append(scale(polynomial, spacing**(basis % 3)))
    return weights


def composed_references(bounds, corners, seed_u, seed_tau, stable_jets):
    with localcontext() as context:
        context.prec = 80
        weights_u = basis_series(seed_u, *bounds[:2])
        weights_tau = basis_series(seed_tau, *bounds[2:])
        source = {}
        for endpoint_tau in range(2):
            for endpoint_u in range(2):
                for derivative_tau in range(3):
                    for derivative_u in range(3):
                        datum = Decimal.from_float(float(corners[endpoint_tau, endpoint_u, derivative_u + 3 * derivative_tau]))
                        product = multiply(weights_u[3 * endpoint_u + derivative_u], weights_tau[3 * endpoint_tau + derivative_tau])
                        source = add(source, scale(product, datum))
        delta_u, delta_tau = series(seed_u), series(seed_tau)
        delta_u[0, 0] = Decimal(0)
        delta_tau[0, 0] = Decimal(0)
        powers_u = [{(0, 0): Decimal(1)}]
        powers_tau = [{(0, 0): Decimal(1)}]
        for order in range(3):
            powers_u.append(multiply(powers_u[-1], delta_u))
            powers_tau.append(multiply(powers_tau[-1], delta_tau))
        lifted = {}
        for (order_u, order_tau), value in zip(ORDERS, stable_jets):
            coefficient = Decimal.from_float(float(value)) / (math.factorial(order_u) * math.factorial(order_tau))
            lifted = add(lifted, scale(multiply(powers_u[order_u], powers_tau[order_tau]), coefficient))
        return derivative_vector(source), derivative_vector(lifted)


def driver_call(executable, inputs, output, label):
    binary = OUTPUT / "build/driver" / executable
    completed = subprocess.run([str(binary)], input=inputs, text=True, capture_output=True, check=True,
                               env=private_environment(output / label), timeout=120)
    (output / (label + ".stdout")).write_text(completed.stdout)
    (output / (label + ".stderr")).write_text(completed.stderr)
    return completed.stdout.splitlines()


def seed_checks(output):
    native = json.loads((BASELINE / "validation/native_records.json").read_text())
    inputs = (BASELINE / "validation/native_queries.txt").read_text()
    (output / "native_queries.txt").write_text(inputs)
    lines = driver_call("cell_driver", inputs, output, "identity")
    real = np.array([np.fromstring(lines[3 * index + 1], sep=" ") for index in range(len(native))])
    ad = np.array([np.fromstring(lines[3 * index + 2], sep=" ") for index in range(len(native))])
    previous = np.array([record["variant"] for record in native])
    reference = np.array([record["reference"] for record in native])
    if not np.array_equal(real, previous):
        raise AssertionError("Frozen053 real evaluation changed")
    real_ad_error = abs(ad - real) / np.maximum(1, abs(real))
    if np.max(real_ad_error) > 3e-15:
        raise AssertionError("Identity-seed AD differs from stable real jets")
    patterns = {
        "identity": ([1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0]),
        "linear_shear": ([.7, -.4, 0, 0, 0, 0, 0, 0, 0], [-.2, 1.1, 0, 0, 0, 0, 0, 0, 0]),
        "nonlinear_all_orders": ([.7, -.4, .2, -.3, .6, -.5, .9, -.8, .4], [-.2, 1.1, -.4, .5, -.6, .7, -.8, .9, -.3]),
        "zero_first_nonlinear": ([0, 0, .2, -.3, .6, -.5, .9, -.8, .4], [0, 0, -.4, .5, -.6, .7, -.8, .9, -.3]),
        "constant": ([0] * 9, [0] * 9),
    }
    query_lines = inputs.splitlines()
    seed_lines, records = [], []
    for index, record in enumerate(native):
        bounds = record["bounds"]
        corners = np.array([np.fromstring(query_lines[6 * index + 2 + corner], sep=" ") for corner in range(4)]).reshape(2, 2, 9)
        for name, (derivatives_u, derivatives_tau) in patterns.items():
            seed_u = [record["coordinate"][0], *derivatives_u]
            seed_tau = [record["coordinate"][1], *derivatives_tau]
            source_reference, stable_reference = composed_references(bounds, corners, seed_u, seed_tau, real[index])
            seed_lines.append(query_lines[6 * index])
            seed_lines.extend(" ".join(format(value, ".17g") for value in seed) for seed in (seed_u, seed_tau))
            seed_lines.extend(query_lines[6 * index + 2:6 * index + 6])
            records.append(dict(native_index=index, pattern=name, seed_u=seed_u, seed_tau=seed_tau,
                                source_reference=source_reference, stable_reference=stable_reference))
    seed_input = "\n".join(seed_lines) + "\n"
    (output / "seed_queries.txt").write_text(seed_input)
    returned = driver_call("jet_driver", seed_input, output, "nonlinear")
    if len(returned) != len(records):
        raise AssertionError("Seed protocol count mismatch")
    for record, line in zip(records, returned):
        raw = np.fromstring(line, sep=" ")
        if len(raw) != 11 or raw[0] != 0 or not np.all(np.isfinite(raw)):
            raise AssertionError("Seed evaluation failed")
        record["actual"] = raw[1:].tolist()
        for reference_kind in ("source", "stable"):
            expected = np.array(record[reference_kind + "_reference"])
            record[reference_kind + "_scaled_error"] = (abs(raw[1:] - expected) / np.maximum(1, abs(expected))).tolist()
        if max(record["stable_scaled_error"]) > 3e-13:
            raise AssertionError("Nonlinear composition differs from independently composed stable jets")
    groups = {name: dict(samples=sum(record["pattern"] == name for record in records),
                        maximum_stable_jet_error=max(max(record["stable_scaled_error"]) for record in records if record["pattern"] == name),
                        maximum_source_polynomial_error=max(max(record["source_scaled_error"]) for record in records if record["pattern"] == name))
              for name in patterns}
    summary = dict(native_points=len(native), identity_real_unchanged_bitwise=True,
                   identity_real_AD_max_scaled=float(np.max(real_ad_error)),
                   identity_AD_source_max_by_partial=np.max(abs(ad - reference) / np.maximum(1, abs(reference)), axis=0).tolist(),
                   composed_records=len(records), composed_partial_entries=10 * len(records), groups=groups,
                   reference_method="80-digit direct Hermite polynomial composition in the two-variable degree3 series ring; separate stable-jet composition",
                   scope="Third-order composition only; same source polynomial and frozen053real jets")
    write_json(output / "seed_records.json", records)
    write_json(output / "seed_results.json", summary)
    return summary


def api_checks(output):
    cases = json.loads((BASELINE / "validation/matched_api_records.json").read_text())
    original = api_client.EXPERIMENT
    api_client.EXPERIMENT = BASELINE
    records = []
    with api_client.EOSClient(output / "api_053") as baseline:
        api_client.EXPERIMENT = OUTPUT
        with api_client.EOSClient(output / "api_059") as variant:
            for case in cases:
                query = (case["scenario"], case["route"], case["logRho"], case["logT"])
                first, second = baseline.query(*query), variant.query(*query)
                if not np.array_equal(first[1:], second[1:], equal_nan=True):
                    raise AssertionError("AD-only change altered actual real-path API result")
                records.append(dict(scenario=case["scenario"], route=case["route"], label=case["label"],
                                    status=int(first[1]), callback_count=int(first[2]), all124_payload_fields_equal=True))
    api_client.EXPERIMENT = original
    summary = dict(queries=len(records), statuses=dict(Counter(str(record["status"]) for record in records)),
                   full_payload_bitwise_equal=True, fallback="none", scope="Bounded actual real-path API repeat; not a stellar or full consumer rerun")
    write_json(output / "api_records.json", records)
    write_json(output / "api_results.json", summary)
    return summary


def main():
    output = OUTPUT / "validation"
    if output.exists():
        raise RuntimeError("Preserve existing059 validation")
    output.mkdir()
    started = time.time()
    before = installed_inventory()
    original = (CODE / "arithmetic_control/freeeos_biquintic.f90").read_text()
    changed = (DIRECTORY / "freeeos_biquintic.f90").read_text()
    real_body = lambda source: source[source.index("   subroutine evaluate_cell_real"):source.index("   end subroutine evaluate_cell_real")]
    if real_body(original) != real_body(changed):
        raise AssertionError("Real evaluator source differs from frozen053")
    seed_summary = seed_checks(output)
    api_summary = api_checks(output)
    if before != installed_inventory():
        raise AssertionError("Installed EOS inventory changed")
    write_json(output / "results.json", dict(seed=seed_summary, api=api_summary, installed_unchanged=True,
               real_source_unchanged=True, elapsed_seconds=time.time() - started,
               source_sha256=sha256(Path(__file__)), driver_sha256=sha256(OUTPUT / "build/driver/candidate_driver")))
    print(json.dumps(dict(seed=seed_summary, api=api_summary, elapsed_seconds=time.time() - started), indent=2))


if __name__ == "__main__":
    main()
