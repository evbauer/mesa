"""Broader blend and limiting-face diagnostics through the actual EOS API."""

import json

import numpy as np

from api_client import EOSClient, FIELDS, LN10, CRAD
from check_precision import ExtendedCache
from runtime_environment import WORKSPACE, EXPERIMENT, installed_inventory, verify_private_cache


def balances(result, point):
    values, density, thermal = result[4:31], result[31:58], result[58:85]
    rho, temperature = np.exp(np.asarray(point) * LN10)
    gas, energy, entropy = np.exp(values[:3])
    radiation = CRAD * temperature**4 / 3
    pressure, pressure_t = gas + radiation, gas * thermal[0] + 4 * radiation
    return [float((rho * energy * density[1] + pressure_t - pressure) / pressure),
            float(temperature * entropy * thermal[2] / (energy * thermal[1]) - 1),
            float(-rho * temperature * entropy * density[2] / pressure_t - 1),
            float(energy * thermal[1] / (temperature * values[FIELDS["Cv"]]) - 1)]


def run():
    before = installed_inventory()
    output = EXPERIMENT / "blend_face_checks"
    cache = ExtendedCache()
    blend_records, face_records = [], []
    with EOSClient(output) as client:
        for log_temp in (4.413, 5.513, 6.013, 6.453):
            for log_q in (-.413, -.213, -.013, .187, .387):
                point = np.array([log_q + 2 * log_temp - 12, log_temp])
                outputs = {(model, route): client.query(model, route, *point)
                           for model in ("off", "candidate") for route in ("component", "final")}
                statuses = {model + "_" + route: int(values[1]) for (model, route), values in outputs.items()}
                record = {"logRho": float(point[0]), "logT": float(point[1]), "statuses": statuses}
                if any(statuses.values()):
                    blend_records.append(record)
                    continue
                candidate, legacy = outputs["candidate", "final"], outputs["off", "final"]
                fraction = candidate[4 + FIELDS["frac_FreeEOS"]]
                record.update(fraction_FreeEOS=float(fraction), callback_count=int(candidate[2]),
                    candidate_component_balances=balances(outputs["candidate", "component"], point),
                    candidate_final_balances=balances(candidate, point), legacy_final_balances=balances(legacy, point))
                component_delta = outputs["candidate", "component"][4:31] - outputs["off", "component"][4:31]
                expected = legacy[4:31] + fraction * component_delta
                record["value_blend_error"] = float(np.max(abs(candidate[4:31] - expected) / np.maximum(1, abs(expected))))
                if not 0 < fraction < 1:
                    blend_records.append(record)
                    continue
                partial_errors = []
                physical = point * LN10
                for axis, start in ((0, 31), (1, 58)):
                    analytic = candidate[start:start + 27]
                    errors = []
                    for step in (1e-4, 3e-5, 1e-5):
                        offset = np.zeros(2)
                        offset[axis] = step
                        shifted = {factor: physical + factor * offset for factor in (-2, -1, 1, 2)}
                        if any(cache.cell(position) != cache.cell(physical) for position in shifted.values()):
                            raise ValueError("Blend finite difference leaves native cell")
                        sample = {factor: client.query("candidate", "final", *(position / LN10))
                                  for factor, position in shifted.items()}
                        if any(result[1] != 0 for result in sample.values()):
                            raise ValueError("Blend finite-difference query failed")
                        finite = (8 * (sample[1][4:31] - sample[-1][4:31]) -
                                  (sample[2][4:31] - sample[-2][4:31])) / (12 * step)
                        scale = np.maximum.reduce([abs(analytic), abs(candidate[4:31]), np.ones(27)])
                        errors.append(abs(finite - analytic) / scale)
                    partial_errors.append({"axis": axis, "best_scaled_error": np.min(errors, axis=0).tolist(),
                                           "errors_by_step": np.asarray(errors).tolist()})
                record["finite_difference"] = partial_errors
                blend_records.append(record)
        manifest = WORKSPACE / "freeeos_biquintic_progress/experiments/038_solver_manifest/queries_v2.json"
        for case in json.loads(manifest.read_text())["cases"]:
            face = np.array([case["logRho"], case["logT"]]) * LN10
            trials = []
            for distance in (1e-5, 1e-6, 1e-7, 1e-8):
                sides = []
                for direction, label in ((-1, "left"), (1, "right")):
                    point = face.copy()
                    if case["direction"] == "u":
                        point[0] += direction * distance
                    else:
                        point += direction * distance * np.array([2., 1.])
                    values = client.query("candidate", "final", *(point / LN10))
                    if values[1] != 0:
                        raise ValueError("Near-face actual API failure")
                    measured = values[[4 + 10, 31 + 10, 58 + 10]]
                    reference = np.array([case["split_" + label][field] for field in ("Cv", "Cv_x", "Cv_t")])
                    scale = np.maximum.reduce([abs(reference), np.full(3, abs(reference[0])), np.ones(3)])
                    sides.append({"side": label, "measured": measured.tolist(), "exact_manifest": reference.tolist(),
                        "scaled_error": (abs(measured - reference) / scale).tolist(),
                        "fraction_FreeEOS": float(values[4 + FIELDS["frac_FreeEOS"]]), "callback_count": int(values[2])})
                trials.append({"normal_distance_natural_log": distance, "sides": sides})
            face_records.append({"case_id": case["case_id"], "direction": case["direction"], "trials": trials})
    verify_private_cache(output)
    if before != installed_inventory():
        raise ValueError("Installed cache or library changed")
    blended = [record for record in blend_records if "finite_difference" in record]
    summary = {"scope": "preselected broader blend grid plus finite offsets from frozen faces; no stellar test",
        "blend_grid_count": len(blend_records), "strict_success_count": sum(not any(record["statuses"].values()) for record in blend_records),
        "blend_count": len(blended), "face_count": len(face_records), "installed_unchanged": True,
        "blend_value_error_max": max(record.get("value_blend_error", 0.) for record in blend_records),
        "blend_fd_error_max": max(max(check["best_scaled_error"]) for record in blended for check in record["finite_difference"]),
        "face_smallest_offset_error_max": max(max(side["scaled_error"]) for record in face_records for side in record["trials"][-1]["sides"]),
        "blend_records": blend_records, "face_records": face_records}
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if not key.endswith("records")}, indent=2))
    if len(blended) < 5 or summary["strict_success_count"] != 20 or summary["blend_value_error_max"] > 1e-12:
        raise AssertionError("Blend value/status gate failed")
    if summary["blend_fd_error_max"] > 2e-5 or summary["face_smallest_offset_error_max"] > 3e-6:
        raise AssertionError("Blend/face derivative gate failed")


if __name__ == "__main__":
    run()
