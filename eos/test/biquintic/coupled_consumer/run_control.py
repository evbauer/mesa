"""Frozen-face pressure plus standard-eps-grav standalone consumer."""

import argparse
import csv
import json
import os
from pathlib import Path
import sys

import numpy as np


CONTROL = Path(__file__).resolve().parent
CODE = CONTROL.parent
WORKSPACE = CONTROL.parents[4]
sys.path.insert(0, str(CODE))

from freeeos_potential import LN10, checksum, query_cells
from gravity_control import gravity_residual
from models import CoupledModels, ideal_state
from solver import Problem, reference, residual, solve
from thermal_consumer.models import CACHE, CORE_SHA, SOURCE, SOURCE_SHA


MANIFEST = WORKSPACE / "freeeos_biquintic_progress/experiments/038_solver_manifest/queries_v2.json"
MANIFEST_SHA = "afc6fe676ba13bf858fc0a611268402b3663f4ad1e0ff991049713fd78788f19"
OUTPUT = WORKSPACE / "freeeos_biquintic_progress/experiments/042_coupled_consumer"
SCALES = (1e-4, 1e-3, 1e-2)
SIGNS = (-1, 1)
RATIOS = (-1, 1)
X_HALF_WIDTH = .08
T_HALF_WIDTH = .015
TOLERANCE = 1e-10
MAX_ITERATIONS = 25
PILOT_IDS = ("u_ordinary_01", "u_tail_01", "tau_ordinary_01", "tau_tail_01")


def verify_manifest(models, manifest):
    if checksum(MANIFEST) != MANIFEST_SHA or manifest["schema"] != "FREEEOS_CONSUMER_MANIFEST_V2":
        raise ValueError("Frozen V2 manifest identity changed")
    if len(manifest["cases"]) != 24 or manifest["candidate"]["X"] != .7 or manifest["candidate"]["Z"] != .02:
        raise ValueError("Unexpected frozen case count or composition")
    if manifest["provenance"]["source_sha256"] != SOURCE_SHA or manifest["provenance"]["core_sha256"] != CORE_SHA:
        raise ValueError("Frozen source or core identity changed")
    for case in manifest["cases"]:
        if not all(models.tables.cells[tuple(case[side + "_cell"])] for side in ("left", "right")):
            raise ValueError("Incomplete case " + case["case_id"])


def face_point(case):
    return np.array([case["logRho"] * LN10, case["logT"] * LN10])


def box_at(face):
    return ((float(face[0] - X_HALF_WIDTH), float(face[0] + X_HALF_WIDTH)),
            (float(face[1] - T_HALF_WIDTH), float(face[1] + T_HALF_WIDTH)))


def nominal_points(case, sign, scale, ratio):
    face = face_point(case)
    shift = np.array([ratio * sign * scale, sign * scale]) / 2
    return face - shift, face + shift


def case_forcing(models, initial, target):
    states = {name: models.at(name, *target) for name in ("candidate", "legacy")}
    pressures = [states[name]["P"] for name in ("candidate", "legacy")]
    pressure_target = float(np.sqrt(pressures[0] * pressures[1]))
    gravities = [gravity_residual(states[name], *(target - initial))["value"]
                 for name in ("candidate", "legacy")]
    heating = float(-np.mean(gravities))
    initial_cv = [models.at(name, *initial)["Cv"] for name in ("candidate", "legacy")]
    energy_scale = float(np.exp(initial[1]) * np.sqrt(initial_cv[0] * initial_cv[1]))
    return pressure_target, heating, energy_scale


def model_problem(models, name, initial, pressure_target, heating, energy_scale, box):
    callback = lambda x, t, cell=None: models.at(name, x, t, cell)
    return Problem(callback, float(initial[0]), float(initial[1]),
                   pressure_target, heating, energy_scale, *box)


def native_coordinates(point):
    density_log, temperature_log = point
    return density_log - 2 * temperature_log + 12 * LN10, temperature_log


def count_crossings(first, second, faces):
    return int(np.count_nonzero((faces > min(first, second) + 1e-12)
                                  & (faces < max(first, second) - 1e-12)))


def crossing_record(models, case, initial, outcome):
    if outcome.get("root") is None:
        return {}
    root = np.array(outcome["root"])
    u_faces, tau_faces = models.tables.logq * LN10, models.tables.logt * LN10
    first, last = native_coordinates(initial), native_coordinates(root)
    accepted = [np.array([item["x"], item["t"]]) for item in outcome["history"]]
    accepted_counts = {"u": 0, "tau": 0}
    for start, stop in zip(accepted[:-1], accepted[1:]):
        first_native, second_native = native_coordinates(start), native_coordinates(stop)
        accepted_counts["u"] += count_crossings(first_native[0], second_native[0], u_faces)
        accepted_counts["tau"] += count_crossings(first_native[1], second_native[1], tau_faces)
    selected_axis = 0 if case["direction"] == "u" else 1
    selected_face = case["logQ"] * LN10 if selected_axis == 0 else case["logT"] * LN10
    return {"straight_initial_root_u_crossings": count_crossings(first[0], last[0], u_faces),
            "straight_initial_root_tau_crossings": count_crossings(first[1], last[1], tau_faces),
            "accepted_u_crossings": accepted_counts["u"],
            "accepted_tau_crossings": accepted_counts["tau"],
            "selected_face_crossed_straight": bool((first[selected_axis] - selected_face)
                                                  * (last[selected_axis] - selected_face) < 0),
            "selected_face_crossed_accepted": any(
                (native_coordinates(start)[selected_axis] - selected_face)
                * (native_coordinates(stop)[selected_axis] - selected_face) < 0
                for start, stop in zip(accepted[:-1], accepted[1:]))}


def run_trial(models, case, family, name, sign, scale, ratio, initial, nominal,
              forcing, manufactured=False):
    box = box_at(face_point(case))
    pressure_target, heating, energy_scale = forcing
    if manufactured:
        state = models.at(name, *nominal)
        pressure_target = state["P"]
        heating = -gravity_residual(state, *(nominal - initial))["value"]
    problem = model_problem(models, name, initial, pressure_target, heating, energy_scale, box)
    record = {"case_id": case["case_id"], "family": family, "model": name,
              "direction": case["direction"], "category": case["category"],
              "sign": sign, "scale": scale, "density_over_temperature_increment": ratio,
              "initial": initial.tolist(), "nominal_target": nominal.tolist(),
              "pressure_target": float(pressure_target), "heating": float(heating),
              "energy_scale": float(energy_scale), "box": box,
              "manufactured_model_specific": manufactured}
    independent = reference(problem)
    newton = solve(problem, initial, tolerance=TOLERANCE, max_iterations=MAX_ITERATIONS)
    record["reference"] = independent
    record["newton"] = newton
    record["status"] = newton["status"] if independent["status"] == "bracketed" else "unverified_" + independent["status"]
    if newton.get("root") is not None and independent.get("root") is not None:
        difference = np.asarray(newton["root"]) - np.asarray(independent["root"])
        record["root_error_inf"] = float(np.max(np.abs(difference)))
        record["root_difference"] = difference.tolist()
    if newton.get("root") is not None and manufactured:
        record["manufactured_target_error_inf"] = float(np.max(np.abs(np.asarray(newton["root"]) - nominal)))
    record["crossings"] = crossing_record(models, case, initial, newton)
    return record


def native_cell(models, point):
    native_u, native_t = native_coordinates(point)
    density_index, temperature_index, inside = query_cells(
        models.tables.logq, models.tables.logt,
        np.array([native_u / LN10]), np.array([native_t / LN10]))
    if not inside[0]:
        return None
    return int(density_index[0]), int(temperature_index[0])


def jacobian_checks(models, case):
    initial, nominal = nominal_points(case, 1, 1e-3, 1)
    forcing = case_forcing(models, initial, nominal)
    face = face_point(case)
    box = box_at(face)
    checks = {"case_id": case["case_id"], "within_cell": [], "one_sided_face": []}
    for name in ("candidate", "legacy"):
        problem = model_problem(models, name, initial, *forcing, box)
        for side, direction in (("left", -1), ("right", 1)):
            center = face.copy()
            center[0 if case["direction"] == "u" else 1] += direction * .001
            center_cell = native_cell(models, center)
            for step in (1e-4, 1e-5, 1e-6):
                value, analytic, _ = residual(problem, *center)
                for axis in range(2):
                    offset = np.zeros(2)
                    offset[axis] = step
                    plus, minus = center + offset, center - offset
                    if native_cell(models, plus) != center_cell or native_cell(models, minus) != center_cell:
                        checks["within_cell"].append({"model": name, "side": side, "axis": axis,
                                                      "step": step, "status": "crosses_another_cell"})
                        continue
                    finite = (residual(problem, *plus)[0] - residual(problem, *minus)[0]) / (2 * step)
                    for equation in range(2):
                        checks["within_cell"].append({"model": name, "side": side, "axis": axis,
                            "equation": equation, "step": step, "status": "checked",
                            "analytic": float(analytic[equation, axis]),
                            "finite_difference": float(finite[equation]),
                            "scaled_error": float((finite[equation] - analytic[equation, axis])
                                                  / max(1, abs(analytic[equation, axis])))})
        if name != "candidate":
            continue
        for side in ("left", "right"):
            cell = case[side + "_cell"]
            face_value, analytic, _ = residual(problem, *face, cell=cell)
            for axis in range(2):
                if case["direction"] == "u":
                    step_sign = (-1 if side == "left" else 1) * (1 if axis == 0 else -1)
                else:
                    step_sign = (-1 if side == "left" else 1) if axis == 1 else 1
                for step in (1e-6, 1e-7):
                    shifted = face.copy()
                    shifted[axis] += step_sign * step
                    finite = (residual(problem, *shifted, cell=cell)[0] - face_value) / (step_sign * step)
                    for equation in range(2):
                        checks["one_sided_face"].append({"side": side, "axis": axis,
                            "equation": equation, "step": step, "analytic": float(analytic[equation, axis]),
                            "finite_difference": float(finite[equation]),
                            "scaled_error": float((finite[equation] - analytic[equation, axis])
                                                  / max(1, abs(analytic[equation, axis])))})
    return checks


def write_queries(records, output):
    fields = ("case_id", "trial_id", "model", "stage", "logRho", "logT")
    with (output / "path_queries.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        for record in records:
            points = [("initial", record["initial"]), ("nominal_target", record["nominal_target"])]
            if record["reference"].get("root") is not None:
                points.append(("reference_root", record["reference"]["root"]))
            points += [("accepted_iterate", [item["x"], item["t"]])
                       for item in record["newton"].get("history", [])]
            points += [("rejected_" + item["reason"], [item["x"], item["t"]])
                       for item in record["newton"].get("rejected_history", [])]
            for stage, point in points:
                if np.all(np.isfinite(point)):
                    writer.writerow((record["case_id"], record["trial_id"], record["model"], stage,
                                     point[0] / LN10, point[1] / LN10))
    with (output / "span_queries.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        for record in records:
            root = record["newton"].get("root")
            if root is None:
                continue
            for point in np.linspace(record["initial"], root, 9):
                if np.all(np.isfinite(point)):
                    writer.writerow((record["case_id"], record["trial_id"], record["model"],
                                     "straight_span_sample", point[0] / LN10, point[1] / LN10))


def make_figure(records, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    common = [record for record in records if record["family"] == "common"]
    by_key = {}
    for record in common:
        key = (record["case_id"], record["sign"], record["scale"],
               record["density_over_temperature_increment"])
        by_key.setdefault(key, {})[record["model"]] = record
    paired = [(key, pair) for key, pair in by_key.items()
              if pair.get("candidate", {}).get("newton", {}).get("iterations") is not None
              and pair.get("legacy", {}).get("newton", {}).get("iterations") is not None]
    figure, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for category, color in (("ordinary", "#277da1"), ("tail", "#d1495b")):
        selected = [(key, pair) for key, pair in paired if pair["candidate"]["category"] == category]
        axes[0].scatter([key[2] for key, _ in selected],
                        [pair["candidate"]["newton"]["iterations"] - pair["legacy"]["newton"]["iterations"]
                         for _, pair in selected], color=color, label=category, alpha=.65, s=22)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("nominal |change in ln T|")
    axes[0].set_ylabel("candidate minus legacy Newton iterations")
    axes[0].legend()
    axes[0].grid(alpha=.2)
    for name, color in (("candidate", "#277da1"), ("legacy", "#eb8a44")):
        selected = [record for record in common if record["model"] == name]
        axes[1].scatter([record["scale"] for record in selected],
                        [max((item["condition"] for item in record["newton"].get("history", [])), default=np.nan)
                         for record in selected], color=color, label=name, alpha=.5, s=20)
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("nominal |change in ln T|")
    axes[1].set_ylabel("max 2x2 Jacobian condition")
    axes[1].legend()
    axes[1].grid(alpha=.2)
    figure.suptitle("Frozen-face coupled pressure + standard eps-grav consumer")
    figure.savefig(output / "coupled_summary.pdf")
    plt.close(figure)


def run(pilot=False):
    models = CoupledModels()
    manifest = json.loads(MANIFEST.read_text())
    verify_manifest(models, manifest)
    selected = [case for case in manifest["cases"] if not pilot or case["case_id"] in PILOT_IDS]
    output = OUTPUT / "pilot" if pilot else OUTPUT
    output.mkdir(parents=True, exist_ok=True)
    records = []
    for case in selected:
        face = face_point(case)
        box = box_at(face)
        for density_log in box[0]:
            for temperature_log in box[1]:
                if not models.valid(density_log, temperature_log):
                    raise ValueError("Prescribed box corner outside split support: " + case["case_id"])
        for scale in SCALES:
            for sign in SIGNS:
                for ratio in RATIOS:
                    initial, nominal = nominal_points(case, sign, scale, ratio)
                    forcing = case_forcing(models, initial, nominal)
                    for name in ("candidate", "legacy"):
                        records.append(run_trial(models, case, "common", name, sign, scale, ratio,
                                                 initial, nominal, forcing))
        if case["case_id"] in PILOT_IDS:
            for ratio in RATIOS:
                initial, nominal = nominal_points(case, 1, 1e-3, ratio)
                forcing = case_forcing(models, initial, nominal)
                for name in ("candidate", "legacy"):
                    records.append(run_trial(models, case, "manufactured", name, 1, 1e-3, ratio,
                                             initial, nominal, forcing, manufactured=True))
        print(case["case_id"] + " complete", flush=True)
    for index, record in enumerate(records, start=1):
        record["trial_id"] = f"C{index:04d}"
    checks = [jacobian_checks(models, case) for case in selected]
    (output / "trials.jsonl").write_text("".join(json.dumps(record, allow_nan=False) + "\n" for record in records))
    (output / "jacobian_checks.json").write_text(json.dumps(checks, indent=2, allow_nan=False) + "\n")
    write_queries(records, output)
    status_counts = {}
    for record in records:
        status_counts[record["status"]] = status_counts.get(record["status"], 0) + 1
    result = {"scope": "standalone fixed-composition pressure plus standard-eps-grav residual; no MESA runtime/stellar call",
              "provenance": {"manifest_sha256": checksum(MANIFEST), "source_sha256": checksum(SOURCE),
                             "core_sha256": checksum(CODE / "freeeos_potential.py"),
                             "legacy_cache_sha256": checksum(CACHE), "script_sha256": checksum(__file__),
                             "node": os.uname().nodename, "job_id": os.environ.get("SLURM_JOB_ID"),
                             "command": sys.argv, "python": sys.version, "numpy": np.__version__},
              "pilot": pilot, "selected_case_ids": [case["case_id"] for case in selected],
              "nominal_scales": SCALES, "signs": SIGNS, "density_over_temperature_ratios": RATIOS,
              "box_half_widths_natural_log": {"x": X_HALF_WIDTH, "t": T_HALF_WIDTH},
              "solver_tolerance_residual_inf": TOLERANCE, "max_iterations": MAX_ITERATIONS,
              "forcing": "shared geometric total-P target and arithmetic mean of negative model epsgrav at symmetric nominal target; common previous state; energy_scale=T0*geometric_mean(Cv0), dt=1",
              "reference": "independent nested pressure-in-x Brent at fixed t, then reduced epsgrav-in-t Brent; 11-point sampled monotonicity diagnostic, not proof of uniqueness",
              "trial_count": len(records), "common_trial_count": sum(r["family"] == "common" for r in records),
              "manufactured_trial_count": sum(r["family"] == "manufactured" for r in records),
              "status_counts": status_counts,
              "path_limit": "CSV records finite sampled states; straight spans and seed purity do not prove continuum MESA selection"}
    (output / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    make_figure(records, output)
    print(json.dumps({"pilot": pilot, "case_count": len(selected), "trial_count": len(records),
                      "status_counts": status_counts}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", action="store_true")
    args = parser.parse_args()
    run(args.pilot)
