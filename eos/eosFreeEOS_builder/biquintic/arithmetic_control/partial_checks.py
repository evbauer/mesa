"""All27 state-partial finite differences, unchanged AD, and exact nodal basis checks."""

from decimal import Decimal
import json
from pathlib import Path
import sys

import numpy as np

from validate import CODE, OUTPUT, EXPERIMENTS, decimal_weights, runner
import api_client


def main():
    interior_extra = "--interior-extra" in sys.argv
    output = OUTPUT / ("partial_validation_interior" if interior_extra else "partial_validation")
    if output.exists():
        raise RuntimeError("Preserve existing partial checks")
    output.mkdir()
    endpoint_checks = []
    for endpoint in range(2):
        for order in range(3):
            weights = decimal_weights(float(endpoint), 0., 1., order)
            for basis, value in enumerate(weights):
                expected = Decimal(int(basis == 3 * endpoint + order))
                endpoint_checks.append(value == expected)
    if not all(endpoint_checks):
        raise AssertionError("Hermite nodal jet constraint changed")
    native = json.loads((OUTPUT / "validation/native_records.json").read_text())
    baseline_lines = (OUTPUT / "validation/baseline_native.stdout").read_text().splitlines()
    variant_lines = (OUTPUT / "validation/variant_native.stdout").read_text().splitlines()
    baseline_ad = np.array([np.fromstring(baseline_lines[3 * index + 2], sep=" ") for index in range(len(native))])
    variant_ad = np.array([np.fromstring(variant_lines[3 * index + 2], sep=" ") for index in range(len(native))])
    if not np.array_equal(baseline_ad, variant_ad):
        raise AssertionError("Unchanged AD path changed numerically")
    reference = np.array([record["reference"] for record in native])
    ad_summary = dict(unchanged_AD_bitwise=True,
        ad_max_scaled_by_partial=np.max(abs(baseline_ad - reference) / np.maximum(1, abs(reference)), axis=0).tolist(),
        variant_real_AD_difference_by_partial=np.max(abs(variant_ad - np.array([record["variant"] for record in native])) /
                                                     np.maximum(1, abs(reference)), axis=0).tolist())
    cases = json.loads((EXPERIMENTS / "038_solver_manifest/queries_v2.json").read_text())["cases"]
    states = []
    for case in cases:
        center = np.array([case["logRho"], case["logT"]]) * runner.LN10
        for side in (-1, 1):
            point = center.copy()
            point[0 if case["direction"] == "u" else 1] += side * .001
            states.append((case["case_id"], point))
    states.extend([("inactive", np.array([4.0031, 7.0017] if interior_extra else [4., 7.]) * runner.LN10),
                   ("blend", np.array([-2.0031, 5.0017] if interior_extra else [-2., 5.]) * runner.LN10)])
    geometry = runner.Geometry()
    records = []
    original = api_client.EXPERIMENT
    with (output / "fd_records.jsonl").open("w") as stream:
        for name, experiment in (("baseline", original), ("variant", OUTPUT)):
            api_client.EXPERIMENT = experiment
            with api_client.EOSClient(output / name) as client:
                for label, center in states:
                    for route in ("component", "final"):
                        nominal = client.query("candidate", route, *(center / runner.LN10))
                        if nominal[1] != 0:
                            raise AssertionError("Unexpected FD center error")
                        for axis in (0, 1):
                            analytic = nominal[31 + 27 * axis:58 + 27 * axis]
                            errors = []
                            steps = []
                            for step in (1e-4, 3e-5, 1e-5, 3e-6):
                                offset = np.zeros(2)
                                offset[axis] = step
                                points = [center - offset, center + offset]
                                crosses_cell = any(geometry.cell(*point) != geometry.cell(*center) for point in points)
                                if (interior_extra or label not in ("inactive", "blend")) and crosses_cell:
                                    raise AssertionError("Native FD stencil crossed another cell")
                                first, second = [client.query("candidate", route, *(point / runner.LN10)) for point in points]
                                if first[1] or second[1]:
                                    raise AssertionError("Unexpected FD stencil status")
                                finite = (second[4:31] - first[4:31]) / (2 * step)
                                errors.append(abs(finite - analytic) / np.maximum(1, abs(analytic)))
                                steps.append(dict(step=step, finite=finite.tolist(), scaled_error=errors[-1].tolist(),
                                                  native_cells=[geometry.cell(*point) for point in points], crosses_native_cell=crosses_cell,
                                                  callback_counts=[int(first[2]), int(second[2])],
                                                  fractions=[first[23:31].tolist(), second[23:31].tolist()]))
                            record = dict(model=name, label=label, route=route, axis=axis, center=center.tolist(),
                                          analytic=analytic.tolist(), steps=steps,
                                          minimum_scaled_by_slot=np.min(errors, axis=0).tolist())
                            stream.write(json.dumps(record) + "\n")
                            records.append(record)
    api_client.EXPERIMENT = original
    summary = dict(endpoint_basis_constraints=len(endpoint_checks), all_endpoint_constraints_exact=True,
                   tensor_consequence="Same9 supplied Hermite nodal jets per corner; no source nodes changed",
                   native_AD=ad_summary, stencils=4 * len(records), derivative_entries=4 * 27 * len(records),
                   finite_difference_error="Minimum over four fixed steps per slot, then maximum over states/axes; no tolerance retuning",
                   additional_states_use_explicit_interior_offsets=interior_extra,
                   cross_cell_stencils=sum(step["crosses_native_cell"] for record in records for step in record["steps"]))
    for name in ("baseline", "variant"):
        for route in ("component", "final"):
            for group in ("declared_cases", "additional_blend_inactive"):
                selected = [record for record in records if record["model"] == name and record["route"] == route and
                            ((record["label"] in ("inactive", "blend")) == (group == "additional_blend_inactive"))]
                worst = np.max([record["minimum_scaled_by_slot"] for record in selected], axis=0)
                summary[name + "/" + route + "/" + group] = dict(maximum_best_scaled_by_slot=worst.tolist())
    runner.write_json(output / "results.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
