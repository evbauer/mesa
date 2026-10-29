"""Frozen-face prescribed-heating thermal consumer; no runtime EOS call."""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

from models import CODE, WORKSPACE, CACHE, SOURCE, CORE_SHA, SOURCE_SHA, SmoothModel, TableModels
from solver import reference, residual, solve
from freeeos_potential import LN10, checksum, evaluate, thermodynamic_partials


MANIFEST = WORKSPACE / "freeeos_biquintic_progress/experiments/038_solver_manifest/queries_v2.json"
MANIFEST_SHA = "afc6fe676ba13bf858fc0a611268402b3663f4ad1e0ff991049713fd78788f19"
OUTPUT = WORKSPACE / "freeeos_biquintic_progress/experiments/039_thermal_consumer"
SCALES = (1e-4, 1e-3, 1e-2)
SIGNS = (-1, 1)
TOLERANCE = 1e-11
MAX_ITERATIONS = 40


def face_log_temperature(case):
    return float(case["logT"] * LN10)


def path_geometry(tables, case):
    face = face_log_temperature(case)
    density_log = case["logRho"] * LN10
    outer_lower, outer_upper = face - .05, face + .05
    u_faces = (density_log + 12 * LN10 - tables.logq * LN10) / 2
    tau_faces = tables.logt * LN10
    native_events = np.r_[u_faces, tau_faces]
    native_events = native_events[(native_events > outer_lower) & (native_events < outer_upper)
                                  & (abs(native_events - face) > 1e-10)]
    events = np.unique(np.r_[outer_lower, outer_upper, native_events, face])
    valid = [tables.valid(density_log, float((left + right) / 2))
             for left, right in zip(events[:-1], events[1:])]
    left_index = max(index for index, (left, right) in enumerate(zip(events[:-1], events[1:]))
                     if right <= face + 1e-12)
    right_index = min(index for index, (left, right) in enumerate(zip(events[:-1], events[1:]))
                      if left >= face - 1e-12)
    if not valid[left_index] or not valid[right_index]:
        raise ValueError("Seed-adjacent cells are invalid on fixed-rho path")
    while left_index > 0 and valid[left_index - 1]:
        left_index -= 1
    while right_index + 1 < len(valid) and valid[right_index + 1]:
        right_index += 1
    lower, upper = float(events[left_index] + 1e-10), float(events[right_index + 1] - 1e-10)
    if not lower < face < upper:
        raise ValueError("No positive-width valid bracket around face")
    samples = [lower, upper]
    for left, right in zip(events[left_index:right_index + 1], events[left_index + 1:right_index + 2]):
        samples.extend(np.linspace(left, right, 7)[1:-1])
    samples.extend(events[left_index + 1:right_index + 1])
    samples.append(face)
    return lower, upper, np.unique(samples), u_faces, tau_faces


def crossings(first, second, faces):
    return int(np.count_nonzero((faces > min(first, second) + 1e-12)
                                  & (faces < max(first, second) - 1e-12)))


def model_callback(tables, name, density_log, smooth_cv):
    if name == "smooth":
        return SmoothModel(smooth_cv)
    return lambda temperature_log: tables.at(name, density_log, temperature_log)


def verify_manifest(tables, manifest):
    if checksum(MANIFEST) != MANIFEST_SHA:
        raise ValueError("Frozen V2 manifest identity changed")
    if manifest["schema"] != "FREEEOS_CONSUMER_MANIFEST_V2" or len(manifest["cases"]) != 24:
        raise ValueError("Unexpected frozen manifest schema/count")
    if manifest["candidate"]["X"] != .7 or manifest["candidate"]["Z"] != .02:
        raise ValueError("Frozen V2 manifest X/Z metadata disagree with audited source header")
    if manifest["provenance"]["source_sha256"] != SOURCE_SHA or manifest["provenance"]["core_sha256"] != CORE_SHA:
        raise ValueError("Frozen manifest source/core identity changed")
    for case in manifest["cases"]:
        if not (tables.cells[tuple(case["left_cell"])] and tables.cells[tuple(case["right_cell"])]):
            raise ValueError(case["case_id"] + " lacks split candidate support")
        for side in ("left", "right"):
            cell = tuple(np.array([index]) for index in case[side + "_cell"])
            sample_q = np.array([case["logQ"]])
            sample_t = np.array([case["logT"]])
            native = evaluate(tables.nodes, tables.complete, tables.logq, tables.logt,
                              sample_q, sample_t, cell_indices=cell, add_radiation=True)
            state = thermodynamic_partials(native, sample_q, sample_t)
            for field in ("Cv", "Cv_t"):
                expected = case["split_" + side][field]
                if abs(float(state[field][0]) - expected) > 1e-8 * max(1, abs(expected)):
                    raise ValueError(case["case_id"] + " exact side mismatch: " + field)


def trial_record(tables, case, family, name, sign, scale, initial, heating,
                 root_target=None, solver_start=None, tolerance=TOLERANCE):
    density_log = case["logRho"] * LN10
    face = face_log_temperature(case)
    lower, upper, samples, u_faces, tau_faces = path_geometry(tables, case)
    callback = model_callback(tables, name, density_log,
                              np.sqrt(tables.at("candidate", density_log, initial)[0]
                                      * tables.at("legacy", density_log, initial)[0]))
    start = initial if solver_start is None else solver_start
    record = {"case_id": case["case_id"], "family": family, "model": name,
              "direction": case["direction"], "category": case["category"],
              "sign": sign, "scale": scale, "initial_y": initial,
              "solver_start_y": start, "tolerance_absolute_R": tolerance, "face_y": face,
              "dt_seconds": 1., "H_erg_per_g_s": heating,
              "bracket": [lower, upper], "target_y": root_target}
    try:
        for point in samples:
            cv, cv_t = callback(float(point))
            if not np.isfinite(cv) or cv <= 0 or not np.isfinite(cv_t):
                raise ValueError("Invalid Cv/Cv_t on sampled valid bracket")
        bracket_root, diagnostic = reference(callback, initial, heating, lower, upper, samples)
        record.update(reference_root=bracket_root, reference_diagnostic=diagnostic)
        if bracket_root is None:
            record["status"] = "reference_bracket_failure"
            return record
        solved = solve(callback, initial, heating, lower, upper, start,
                       tolerance=tolerance, max_iterations=MAX_ITERATIONS)
        record.update(status=solved.status, newton_root=solved.root,
                      iterations=solved.iterations, rejected=solved.rejected,
                      bisected=solved.bisected, history=solved.history,
                      rejected_history=solved.rejected_history,
                      root_error=None if solved.root is None else solved.root - bracket_root,
                      root_residual=None if solved.root is None else residual(callback, solved.root, initial, heating)[0],
                      target_error=None if root_target is None or solved.root is None else solved.root - root_target,
                      initial_to_root_u_crossings=crossings(initial, bracket_root, u_faces),
                      initial_to_root_tau_crossings=crossings(initial, bracket_root, tau_faces),
                      crossed_selected_face=bool((initial - face) * (bracket_root - face) < 0
                                                 or abs(bracket_root - face) < 1e-10))
        if name == "smooth":
            smooth_cv = callback(initial)[0]
            exact_temperature = np.exp(initial) + heating / smooth_cv
            record["continuous_ode_exact_root"] = float(np.log(exact_temperature))
            record["backward_euler_temporal_error"] = float(bracket_root - np.log(exact_temperature))
        history_y = [item["y"] for item in solved.history]
        record["accepted_u_crossings"] = sum(crossings(first, second, u_faces)
                                              for first, second in zip(history_y[:-1], history_y[1:]))
        record["accepted_tau_crossings"] = sum(crossings(first, second, tau_faces)
                                                for first, second in zip(history_y[:-1], history_y[1:]))
        return record
    except (ValueError, OverflowError, FloatingPointError) as error:
        record["status"] = "invalid_or_exception: " + str(error)
        return record


def jacobian_checks(tables, case):
    density_log = case["logRho"] * LN10
    face = face_log_temperature(case)
    _, _, _, u_faces, tau_faces = path_geometry(tables, case)
    neighbors = np.sort(np.unique(np.r_[u_faces, tau_faces]))
    below = float(face - neighbors[neighbors < face - 1e-11][-1])
    above = float(neighbors[neighbors > face + 1e-11][0] - face)
    cv_face = tables.at("candidate", density_log, face)[0]
    initial = face - .005
    heating = .005 * np.exp(face) * cv_face
    data = {"case_id": case["case_id"], "within_cell": [], "exact_face": {}}
    for name in ("candidate", "legacy", "smooth"):
        callback = model_callback(tables, name, density_log, cv_face)
        for side, distance in (("minus", below), ("plus", above)):
            center = face + (-1 if side == "minus" else 1) * min(.002, distance / 3)
            analytic = residual(callback, center, initial, heating)[1]
            for relative_step in (.1, .01, .001):
                step = relative_step * min(.002, distance / 3)
                finite = (residual(callback, center + step, initial, heating)[0]
                          - residual(callback, center - step, initial, heating)[0]) / (2 * step)
                data["within_cell"].append({"model": name, "side": side, "step": step,
                                            "analytic": analytic, "finite_difference": finite,
                                            "scaled_error": (finite - analytic) / max(1, abs(analytic))})
    for side, cell in (("left", case["left_cell"]), ("right", case["right_cell"])):
        sample_q = np.array([case["logQ"]])
        sample_t = np.array([case["logT"]])
        native = evaluate(tables.nodes, tables.complete, tables.logq, tables.logt,
                          sample_q, sample_t,
                          cell_indices=tuple(np.array([index]) for index in cell),
                          add_radiation=True)
        state = thermodynamic_partials(native, sample_q, sample_t)
        cv, cv_t = float(state["Cv"][0]), float(state["Cv_t"][0])
        increment = heating / (np.exp(face) * cv)
        data["exact_face"][side] = {"Cv": cv, "Cv_t": cv_t,
                                    "residual_jacobian": 1 + increment * (1 + cv_t / cv)}
    face_value = residual(lambda temperature_log: (cv_face, 0.), face, initial, heating)[0]
    for side, distance in (("minus", below), ("plus", above)):
        for step in (min(1e-5, distance / 100), min(1e-6, distance / 1000)):
            callback = model_callback(tables, "candidate", density_log, cv_face)
            shifted = face + (-step if side == "minus" else step)
            shifted_value = residual(callback, shifted, initial, heating)[0]
            finite = ((face_value - shifted_value) / step if side == "minus"
                      else (shifted_value - face_value) / step)
            physical_side = ("right" if side == "minus" else "left") if case["direction"] == "u" else side.replace("minus", "left").replace("plus", "right")
            data["exact_face"].setdefault("one_sided", []).append(
                {"y_side": side, "native_side": physical_side, "step": step,
                 "finite_difference": finite,
                 "exact_jacobian": data["exact_face"][physical_side]["residual_jacobian"]})
    return data


def make_plots(records, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    common = [record for record in records if record["family"] == "common"
              and record["status"] == "converged"]
    figure, axis = plt.subplots(figsize=(8, 4.7))
    for name, color in (("candidate", "tab:blue"), ("legacy", "tab:orange"), ("smooth", "tab:green")):
        selected = [record for record in common if record["model"] == name]
        axis.scatter([record["scale"] for record in selected], [record["iterations"] for record in selected],
                     c=color, label=name, alpha=.6, s=26)
    axis.set_xscale("log")
    axis.set_xlabel("Nominal |change in ln T|")
    axis.set_ylabel("Newton iterations")
    axis.grid(alpha=.2)
    axis.legend()
    figure.tight_layout()
    figure.savefig(output / "iterations.pdf")
    plt.close(figure)
    tail = [record for record in common if record["model"] == "candidate" and record["category"] == "tail"]
    if tail:
        worst = max(tail, key=lambda record: record["iterations"])
        figure, axis = plt.subplots(figsize=(7, 4))
        axis.semilogy(range(len(worst["history"])),
                      [max(abs(item["residual"]), 1e-16) for item in worst["history"]], marker="o")
        axis.set_xlabel("Accepted iterate")
        axis.set_ylabel("|R|")
        axis.set_title(worst["case_id"] + ", sign=" + str(worst["sign"]) + ", scale=" + str(worst["scale"]))
        axis.grid(alpha=.2)
        figure.tight_layout()
        figure.savefig(output / "tail_residual_history.pdf")
        plt.close(figure)


def run():
    output = OUTPUT
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(MANIFEST.read_text())
    tables = TableModels()
    verify_manifest(tables, manifest)
    records, checks = [], []
    for case in manifest["cases"]:
        checks.append(jacobian_checks(tables, case))
        face = face_log_temperature(case)
        density_log = case["logRho"] * LN10
        for sign in SIGNS:
            for scale in SCALES:
                initial = face - sign * scale / 2
                candidate_cv = tables.at("candidate", density_log, initial)[0]
                legacy_cv = tables.at("legacy", density_log, initial)[0]
                common_cv = np.sqrt(candidate_cv * legacy_cv)
                heating = sign * scale * np.exp(initial) * common_cv
                for name in ("candidate", "legacy", "smooth"):
                    records.append(trial_record(tables, case, "common", name, sign,
                                                scale, initial, heating))
                off_face_initial = face + sign * 2 * scale
                off_candidate_cv = tables.at("candidate", density_log, off_face_initial)[0]
                off_legacy_cv = tables.at("legacy", density_log, off_face_initial)[0]
                off_cv = np.sqrt(off_candidate_cv * off_legacy_cv)
                off_heating = sign * scale * np.exp(off_face_initial) * off_cv
                for name in ("candidate", "legacy", "smooth"):
                    records.append(trial_record(tables, case, "off_face_control", name, sign,
                                                scale, off_face_initial, off_heating))
        if case["category"] == "tail" and case["case_id"].endswith("01"):
            for sign in SIGNS:
                initial = face - sign * .005
                for offset in (-1e-6, 0., 1e-6):
                    target = face + offset
                    for name in ("candidate", "legacy", "smooth"):
                        cv = (tables.at(name, density_log, target)[0] if name != "smooth"
                              else np.sqrt(tables.at("candidate", density_log, initial)[0]
                                           * tables.at("legacy", density_log, initial)[0]))
                        heating = (target - initial) * np.exp(target) * cv
                        records.append(trial_record(tables, case, "manufactured", name, sign,
                                                    .005, initial, heating, root_target=target))
    selected_controls = {"u_ordinary_01", "u_tail_01", "tau_ordinary_01", "tau_tail_01"}
    for case in manifest["cases"]:
        if case["case_id"] not in selected_controls:
            continue
        face = face_log_temperature(case)
        density_log = case["logRho"] * LN10
        scale = 1e-2
        for sign in SIGNS:
            initial = face - sign * scale / 2
            candidate_cv = tables.at("candidate", density_log, initial)[0]
            legacy_cv = tables.at("legacy", density_log, initial)[0]
            common_cv = np.sqrt(candidate_cv * legacy_cv)
            heating = sign * scale * np.exp(initial) * common_cv
            for name in ("candidate", "legacy", "smooth"):
                for offset in (-.02, .02):
                    records.append(trial_record(tables, case, "restart_guess", name, sign,
                                                scale, initial, heating,
                                                solver_start=face + offset))
                for tolerance in (1e-9, 1e-11, 1e-13):
                    records.append(trial_record(tables, case, "tolerance_sweep", name, sign,
                                                scale, initial, heating, tolerance=tolerance))
    for trial_index, record in enumerate(records, start=1):
        record["trial_id"] = f"T{trial_index:04d}"
    output.joinpath("trials.jsonl").write_text("".join(json.dumps(record, allow_nan=False) + "\n" for record in records))
    output.joinpath("jacobian_checks.json").write_text(json.dumps(checks, indent=2, allow_nan=False) + "\n")
    case_by_id = {case["case_id"]: case for case in manifest["cases"]}
    with output.joinpath("path_queries.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("case_id", "trial_id", "model", "stage", "logRho", "logT"))
        for record in records:
            points = [("initial", record["initial_y"]), ("solver_start", record["solver_start_y"])]
            if record.get("reference_root") is not None:
                points.append(("reference_root", record["reference_root"]))
            points += [("accepted_iterate", item["y"]) for item in record.get("history", [])]
            points += [("rejected_" + item["reason"], item["y"])
                       for item in record.get("rejected_history", [])]
            for stage, point in points:
                writer.writerow((record["case_id"], record["trial_id"], record["model"], stage,
                                 case_by_id[record["case_id"]]["logRho"], point / LN10))
    with output.joinpath("span_queries.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("case_id", "trial_id", "model", "stage", "logRho", "logT"))
        for record in records:
            if record["status"] != "converged":
                continue
            for point in np.linspace(record["initial_y"], record["newton_root"], 9):
                writer.writerow((record["case_id"], record["trial_id"], record["model"],
                                 "span_sample", case_by_id[record["case_id"]]["logRho"], point / LN10))
    statuses = {}
    for record in records:
        key = record["status"]
        statuses[key] = statuses.get(key, 0) + 1
    summary = {"claim": "Standalone fixed-rho prescribed-H table-model solver only; root agreement is not source accuracy or MESA runtime evidence",
               "command": [sys.executable, *sys.argv], "mesa_head": "e92a7401",
               "manifest_sha256": checksum(MANIFEST), "source_sha256": checksum(SOURCE),
               "core_sha256": checksum(CODE / "freeeos_potential.py"), "legacy_cache_sha256": checksum(CACHE),
               "candidate": manifest["candidate"], "case_count": len(manifest["cases"]),
               "source_header_X_Z": [.7, .02],
               "manifest_v1_superseded_sha256": "74f5079a7b5830abeae2947c4063d78e214a9dde2be43dac59f507d6b6007eb8",
               "manifest_v1_correction": "V1 X/Z metadata labels reversed; V2 changed metadata only, with identical 24 case records; source was always X=.70,Z=.02",
               "common_trial_count": sum(record["family"] == "common" for record in records),
               "off_face_control_trial_count": sum(record["family"] == "off_face_control" for record in records),
               "manufactured_trial_count": sum(record["family"] == "manufactured" for record in records),
               "restart_guess_trial_count": sum(record["family"] == "restart_guess" for record in records),
               "tolerance_sweep_trial_count": sum(record["family"] == "tolerance_sweep" for record in records),
               "status_counts": statuses, "tolerance_absolute_R": TOLERANCE,
               "max_iterations": MAX_ITERATIONS, "dt_seconds": 1.,
               "scales": SCALES, "signs": SIGNS,
               "uniqueness_limit": "Positive sampled Jacobians and sign-change counts do not prove global uniqueness between samples.",
               "MESA_fraction_limit": "Pure FreeEOS at frozen seed only; path_queries.csv and separate span_queries.csv need independent runtime probing; finite span samples do not prove continuum routing."}
    output.joinpath("results.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    make_plots(records, output)
    print(json.dumps({"status_counts": statuses, "trials": len(records),
                      "manifest_sha256": summary["manifest_sha256"]}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    run()
