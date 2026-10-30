"""Repeat frozen scalar/coupled forcing through actual component and final APIs."""

import argparse
from collections import Counter
from dataclasses import asdict
from datetime import datetime
import json
import os
from pathlib import Path
import resource
import socket
import sys
import time

import numpy as np


DIRECTORY = Path(__file__).resolve().parent
CODE = DIRECTORY.parent
WORKSPACE = DIRECTORY.parents[4]
OUTPUT = WORKSPACE / "freeeos_biquintic_progress/experiments/051_runtime_consumers"
EXPERIMENTS = OUTPUT.parent
sys.path.insert(0, str(CODE))
sys.path.insert(0, str(CODE / "runtime_bridge"))

from api_client import EOSClient, LN10, IDENTITY
from runtime_environment import installed_inventory, verify_private_cache, sha256
from thermal_consumer import solver as scalar_solver
from coupled_consumer import solver as coupled_solver


PINNED = {
    EXPERIMENTS / "039_thermal_consumer/trials.jsonl": "41a2f4dcdb6caf6ee520edd5fa4aa3f9a2f087e61a78fb821b35bd62bcc1eb50",
    EXPERIMENTS / "042_coupled_consumer/trials.jsonl": "4a2dbc3dd80494a3438955e5d2e0271e6ea7f30e2a6902380aefbc87743d908d",
    EXPERIMENTS / "038_solver_manifest/queries_v2.json": "afc6fe676ba13bf858fc0a611268402b3663f4ad1e0ff991049713fd78788f19",
    CODE / "thermal_consumer/solver.py": "9782b8271245fa132638f4345334981603ea652ade256620ad18dae8b647768b",
    CODE / "coupled_consumer/solver.py": "bc23fd6c1900013ad20a3f64b95789c05b0fbdd55d1d246fd7d404026aeb6f20",
    CODE / "gravity_control.py": "3c06d8b5eea792fba523fe77a98306f74218ab82a2f8d4134bbdb50819f80f0c",
}
PILOT = {"u_ordinary_01", "u_tail_01", "tau_ordinary_01", "tau_tail_01"}


def serializable(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(type(value).__name__)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, default=serializable, allow_nan=False) + "\n")


def population(fractions):
    if fractions[2] >= 1 - 1e-12 and np.max(np.abs(np.delete(fractions, 2))) <= 1e-12:
        return "pure_FreeEOS"
    if fractions[2] <= 1e-12:
        return "inactive_FreeEOS"
    return "blend"


class LoggedClient(EOSClient):
    def __init__(self, output):
        super().__init__(output)
        self.evaluations = (self.output / "evaluations.jsonl").open("w")
        self.context = {}
        self.stage = "setup"
        self.stats = Counter()
        self.trial_stats = Counter()
        self.perturb_abundances = False

    def query(self, scenario, route, log_rho, log_temp):
        key = (scenario, route, float(log_rho), float(log_temp))
        cached = key in self.memo
        returned = super().query(scenario, route, log_rho, log_temp)
        status = int(returned[1])
        record = dict(self.context, stage=self.stage, scenario=scenario, route=route,
                      logRho=log_rho, logT=log_temp, actual_query_row=int(returned[0]), memoized=cached,
                      ierr=status, callback_count=int(returned[2]), hybrid_retained=int(returned[3]))
        self.stats["logical_evaluations"] += 1
        self.trial_stats["logical_evaluations"] += 1
        if status:
            label = f"error_{status}"
            if returned[2] > 0:
                self.stats["callback_errors"] += 1
        else:
            fractions = returned[23:31]
            record["fractions"] = fractions.tolist()
            label = population(fractions)
        record["population"] = label
        self.stats[f"{route}/{scenario}/{label}"] += 1
        self.trial_stats[label] += 1
        self.evaluations.write(json.dumps(record, default=serializable, allow_nan=False) + "\n")
        if self.perturb_abundances and status == 0:
            returned[85:] += 1e9 * np.arange(1, 41)
        return returned

    def begin(self, record):
        self.context = {key: record[key] for key in ("consumer", "trial_id", "case_id", "model", "route")}
        self.trial_stats = Counter()

    def close(self):
        try:
            super().close()
        finally:
            self.evaluations.close()


class Geometry:
    def __init__(self):
        lines = IDENTITY.read_text().splitlines()
        nq, nt = map(int, lines[1].split())
        umin, du, tmin, dt = map(float, lines[5].split())
        self.u = umin + du * np.arange(nq)
        self.t = tmin + dt * np.arange(nt)

    def cell(self, density, temperature):
        native = density - 2 * temperature + 12 * LN10
        if not self.u[0] <= native <= self.u[-1] or not self.t[0] <= temperature <= self.t[-1]:
            return None
        return (min(len(self.u) - 2, int(np.searchsorted(self.u, native, side="right") - 1)),
                min(len(self.t) - 2, int(np.searchsorted(self.t, temperature, side="right") - 1)))

    def scalar_samples(self, case, bracket):
        face = case["logT"] * LN10
        density = case["logRho"] * LN10
        outer_lower, outer_upper = face - .05, face + .05
        u_events = (density + 12 * LN10 - self.u) / 2
        events = np.r_[u_events, self.t]
        events = events[(events > outer_lower) & (events < outer_upper) & (abs(events - face) > 1e-10)]
        events = np.unique(np.r_[outer_lower, events, face, outer_upper])
        selected = events[(events >= bracket[0] - 2e-10) & (events <= bracket[1] + 2e-10)]
        points = [*bracket, face]
        points.extend(selected[1:-1])
        for lower, upper in zip(selected[:-1], selected[1:]):
            points.extend(np.linspace(lower, upper, 7)[1:-1])
        return np.unique(points)

    def crossings(self, case, history, scalar=False):
        if scalar:
            points = [(case["logRho"] * LN10, item["y"]) for item in history]
        else:
            points = [(item["x"], item["t"]) for item in history]
        native = [(density - 2 * temperature + 12 * LN10, temperature) for density, temperature in points]
        counts = {"u": 0, "tau": 0}
        for first, second in zip(native[:-1], native[1:]):
            for axis, (name, knots) in enumerate((("u", self.u), ("tau", self.t))):
                counts[name] += int(np.count_nonzero((knots > min(first[axis], second[axis]) + 1e-12) &
                                                   (knots < max(first[axis], second[axis]) - 1e-12)))
        axis = 0 if case["direction"] == "u" else 1
        face = case["logQ" if axis == 0 else "logT"] * LN10
        counts["selected_face_crossed"] = any((first[axis] - face) * (second[axis] - face) < 0
                                               for first, second in zip(native[:-1], native[1:]))
        return counts


def scalar_trial(client, original, case, route, geometry, settings):
    record = dict(original, consumer="scalar", route=route, source_trial_id=original["trial_id"])
    client.begin(record)
    density = case["logRho"] * LN10

    def callback(temperature):
        state = client.state(original["model"], route, density, temperature)
        return state["Cv"], state["Cv_t"]

    heating = original["H_erg_per_g_s"] * original["dt_seconds"]
    client.stage = "independent_reference"
    reference = None
    diagnostic = {}
    try:
        samples = geometry.scalar_samples(case, original["bracket"])
        reference, diagnostic = scalar_solver.reference(callback, original["initial_y"], heating,
                                                        *original["bracket"], samples)
    except (ValueError, OverflowError, FloatingPointError) as error:
        diagnostic = {"error": str(error)}
    client.stage = "newton_evaluation"
    outcome = scalar_solver.solve(callback, original["initial_y"], heating, *original["bracket"],
                                  original["solver_start_y"], tolerance=original["tolerance_absolute_R"],
                                  max_iterations=settings["max_iterations"])
    record = {key: record[key] for key in ("consumer", "route", "case_id", "trial_id", "source_trial_id", "family",
                                          "model", "direction", "category", "sign", "scale", "initial_y",
                                          "solver_start_y", "tolerance_absolute_R", "dt_seconds", "H_erg_per_g_s",
                                          "bracket", "target_y")}
    record.update(newton=asdict(outcome), reference_root=reference, reference_diagnostic=diagnostic,
                  status=outcome.status if reference is not None else "unverified_" + outcome.status,
                  original_status=original["status"], original_iterations=original["iterations"],
                  iteration_difference=outcome.iterations - original["iterations"],
                  original_reference_root=original["reference_root"],
                  root_error=None if outcome.root is None or reference is None else outcome.root - reference,
                  reference_shift=None if reference is None else reference - original["reference_root"],
                  newton_shift=None if outcome.root is None else outcome.root - original["newton_root"],
                  crossings=geometry.crossings(case, outcome.history, scalar=True),
                  evaluation_counts=dict(client.trial_stats), fallback_count=0)
    return record


def coupled_problem(client, original, route):
    def callback(density, temperature, cell=None):
        if cell is not None:
            raise ValueError("Actual API has no forced-cell evaluation")
        return client.state(original["model"], route, density, temperature)
    return coupled_solver.Problem(callback, *original["initial"], original["pressure_target"],
                                  original["heating"], original["energy_scale"], *original["box"])


def coupled_trial(client, original, case, route, geometry, settings):
    record = {key: original[key] for key in ("case_id", "trial_id", "family", "model", "direction", "category", "sign",
                                            "scale", "density_over_temperature_increment", "initial", "nominal_target",
                                            "pressure_target", "heating", "energy_scale", "box", "manufactured_model_specific")}
    record.update(consumer="coupled", route=route, source_trial_id=original["trial_id"])
    client.begin(record)
    problem = coupled_problem(client, original, route)
    client.stage = "independent_reference"
    reference = coupled_solver.reference(problem)
    client.stage = "newton_evaluation"
    start = [original["newton"]["history"][0][name] for name in ("x", "t")]
    newton = coupled_solver.solve(problem, start, tolerance=settings["solver_tolerance_residual_inf"],
                                 max_iterations=settings["max_iterations"])
    record.update(solver_start=start, reference=reference, newton=newton,
                  status=newton["status"] if reference["status"] == "bracketed" else "unverified_" + newton["status"],
                  original_status=original["status"], original_iterations=original["newton"]["iterations"],
                  iteration_difference=newton.get("iterations", 0) - original["newton"]["iterations"],
                  original_reference=original["reference"], crossings=geometry.crossings(case, newton["history"]),
                  evaluation_counts=dict(client.trial_stats), fallback_count=0)
    if "root" in reference:
        record["reference_shift_inf"] = float(np.max(np.abs(np.array(reference["root"]) - original["reference"]["root"])))
    if "root" in newton and "root" in reference:
        record["root_error_inf"] = float(np.max(np.abs(np.array(newton["root"]) - reference["root"])))
        record["newton_shift_inf"] = float(np.max(np.abs(np.array(newton["root"]) - original["newton"]["root"])))
    return record


def write_paths(stream, record, case):
    base = {key: record[key] for key in ("consumer", "route", "model", "case_id", "trial_id")}
    if record["consumer"] == "scalar":
        density = case["logRho"] * LN10
        entries = [("physical_initial", density, record["initial_y"]),
                   ("solver_start", density, record["solver_start_y"])]
        entries += [("accepted_iterate", density, item["y"]) for item in record["newton"]["history"]]
        entries += [("rejected_proposal:" + item["reason"], density, item["y"])
                    for item in record["newton"]["rejected_history"]]
        if record["reference_root"] is not None:
            entries.append(("independent_reference_root", density, record["reference_root"]))
    else:
        entries = [("physical_initial", *record["initial"]), ("solver_start", *record["solver_start"])]
        entries += [("accepted_iterate", item["x"], item["t"]) for item in record["newton"]["history"]]
        entries += [("rejected_proposal:" + item["reason"], item["x"], item["t"])
                    for item in record["newton"]["rejected_history"]]
        if "root" in record["reference"]:
            entries.append(("independent_reference_root", *record["reference"]["root"]))
    for stage, density, temperature in entries:
        stream.write(json.dumps(dict(base, stage=stage, logRho=density / LN10, logT=temperature / LN10)) + "\n")


def jacobian_checks(client, scalar, coupled, cases, geometry):
    records = []
    for consumer, originals in (("scalar", scalar), ("coupled", coupled)):
        representative = {}
        for original in originals:
            if original["family"] == "common" and original["sign"] == 1 and original["scale"] == .001:
                representative.setdefault((original["case_id"], original["model"]), original)
        for route in ("component", "final"):
            for (case_id, model), original in representative.items():
                case = cases[case_id]
                client.begin(dict(consumer=consumer, trial_id="FD_" + original["trial_id"], case_id=case_id,
                                  model=model, route=route))
                client.stage = "jacobian_check"
                face = np.array([case["logRho"] * LN10, case["logT"] * LN10])
                if consumer == "scalar":
                    def evaluate(point):
                        def callback(temperature):
                            state = client.state(model, route, face[0], temperature)
                            return state["Cv"], state["Cv_t"]
                        value, derivative = scalar_solver.residual(callback, point[1], original["initial_y"],
                                                                  original["H_erg_per_g_s"] * original["dt_seconds"])
                        return np.array([value]), np.array([[0., derivative]])
                    axes = (1,)
                else:
                    problem = coupled_problem(client, original, route)
                    def evaluate(point):
                        value, jacobian, _ = coupled_solver.residual(problem, *point)
                        return value, jacobian
                    axes = (0, 1)
                for kind in ("within_cell", "one_sided_interior"):
                    for side in (-1, 1):
                        center = face.copy()
                        normal_axis = 1 if consumer == "scalar" or case["direction"] == "tau" else 0
                        center[normal_axis] += side * (.001 if kind == "within_cell" else 1e-8)
                        for axis in axes:
                            sign = side if axis == normal_axis else (1 if case["direction"] == "tau" else -side)
                            for step in ((1e-5, 1e-6) if kind == "within_cell" else (1e-6, 1e-7)):
                                offset = np.zeros(2)
                                offset[axis] = step
                                points = ([center - offset, center + offset] if kind == "within_cell"
                                          else [center + sign * offset, center + 2 * sign * offset])
                                record = dict(consumer=consumer, route=route, case_id=case_id, model=model, kind=kind,
                                              side=side, axis=axis, step=step, center=center.tolist(),
                                              native_cell=geometry.cell(*center), interior_offset=1e-8 if kind != "within_cell" else .001)
                                if any(geometry.cell(*point) != geometry.cell(*center) for point in points):
                                    record["status"] = "crosses_other_native_cell"
                                    records.append(record)
                                    continue
                                try:
                                    value, analytic = evaluate(center)
                                    first, second = (evaluate(point)[0] for point in points)
                                    finite = ((second - first) / (2 * step) if kind == "within_cell"
                                              else sign * (-3 * value + 4 * first - second) / (2 * step))
                                    error = (finite - analytic[:, axis]) / np.maximum(1., np.abs(analytic[:, axis]))
                                    record.update(status="checked", analytic=analytic[:, axis].tolist(),
                                                  finite_difference=finite.tolist(), scaled_error=error.tolist(),
                                                  max_scaled_error=float(np.max(np.abs(error))))
                                except (ValueError, OverflowError, FloatingPointError) as error:
                                    record.update(status="api_or_arithmetic_error", error=str(error))
                                records.append(record)
    return records


def hybrid_checks(client, scalar, coupled, cases):
    checks = []
    for route in ("component", "final"):
        for model in ("candidate", "legacy"):
            thermal = next(record for record in scalar if record["model"] == model)
            paired = next(record for record in coupled if record["model"] == model)
            case = cases[thermal["case_id"]]
            client.begin(dict(consumer="hybrid_control", route=route, model=model, case_id=case["case_id"], trial_id="HYBRID"))
            client.stage = "finite_abundance_perturbation"
            def thermal_value():
                def callback(temperature):
                    state = client.state(model, route, case["logRho"] * LN10, temperature)
                    return state["Cv"], state["Cv_t"]
                return scalar_solver.residual(callback, thermal["initial_y"], thermal["initial_y"], thermal["H_erg_per_g_s"])
            problem = coupled_problem(client, paired, route)
            baseline_scalar = thermal_value()
            baseline_coupled = coupled_solver.residual(problem, *paired["initial"])
            client.perturb_abundances = True
            try:
                changed_scalar = thermal_value()
                changed_coupled = coupled_solver.residual(problem, *paired["initial"])
            finally:
                client.perturb_abundances = False
            checks.append(dict(route=route, model=model,
                               scalar_unchanged=baseline_scalar == changed_scalar,
                               coupled_unchanged=np.array_equal(baseline_coupled[0], changed_coupled[0]) and
                                                 np.array_equal(baseline_coupled[1], changed_coupled[1]),
                               scope="Perturb all40 exported abundance partials after API return before consumer formatting; fixed state"))
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", action="store_true")
    arguments = parser.parse_args()
    if os.environ.get("SLURMD_NODENAME") != socket.gethostname().split(".")[0]:
        raise RuntimeError("Inherited compute-node placement required")
    if int(os.environ.get("SLURM_JOB_END_TIME", "0")) - time.time() < 1800:
        raise RuntimeError("Insufficient remaining allocation for consumer run and checkpoint")
    output = OUTPUT / ("pilot" if arguments.pilot else "full")
    output.mkdir(parents=True, exist_ok=True)
    for path, expected in PINNED.items():
        if sha256(path) != expected:
            raise ValueError("Frozen input changed: " + str(path))
    manifest = json.loads((EXPERIMENTS / "038_solver_manifest/queries_v2.json").read_text())
    cases = {case["case_id"]: case for case in manifest["cases"]}
    scalar = [json.loads(line) for line in (EXPERIMENTS / "039_thermal_consumer/trials.jsonl").read_text().splitlines()]
    scalar = [record for record in scalar if record["model"] in ("candidate", "legacy")]
    coupled = [json.loads(line) for line in (EXPERIMENTS / "042_coupled_consumer/trials.jsonl").read_text().splitlines()]
    if arguments.pilot:
        scalar = [record for record in scalar if record["case_id"] in PILOT]
        coupled = [record for record in coupled if record["case_id"] in PILOT]
    settings_scalar = json.loads((EXPERIMENTS / "039_thermal_consumer/results.json").read_text())
    settings_coupled = json.loads((EXPERIMENTS / "042_coupled_consumer/results.json").read_text())
    geometry = Geometry()
    before = installed_inventory()
    provenance = dict(started=datetime.now().isoformat(), hostname=socket.gethostname(), job=os.environ.get("SLURM_JOB_ID"),
                      planned_trials=2 * (len(scalar) + len(coupled)), scalar_per_route=len(scalar), coupled_per_route=len(coupled),
                      input_sha256={str(path.relative_to(WORKSPACE)): expected for path, expected in PINNED.items()},
                      source_sha256=sha256(Path(__file__)), client_sha256=sha256(CODE / "runtime_bridge/api_client.py"),
                      identity_sha256=sha256(IDENTITY), installed_before=before, gate="lead cleared 18:52 PDT, bounded consumers only",
                      fallback="none; strict statuses retained", virtual_memory_limit_bytes=resource.getrlimit(resource.RLIMIT_AS)[0])
    write_json(output / "provenance.json", provenance)
    started = time.time()
    records = []
    with LoggedClient(output / "api") as client, (output / "trials.jsonl").open("w") as stream, (output / "states.jsonl").open("w") as paths:
        provenance["driver_sha256"] = sha256(client.command[0])
        provenance["driver_command"] = client.command
        write_json(output / "provenance.json", provenance)
        for route in ("component", "final"):
            for consumer, originals, settings, run_trial in (("scalar", scalar, settings_scalar, scalar_trial),
                                                            ("coupled", coupled, settings_coupled, coupled_trial)):
                for original in originals:
                    if int(os.environ["SLURM_JOB_END_TIME"]) - time.time() < 1500:
                        raise RuntimeError("Checkpoint deadline reached; completed JSONL records are saved")
                    record = run_trial(client, original, cases[original["case_id"]], route, geometry, settings)
                    records.append(record)
                    stream.write(json.dumps(record, default=serializable, allow_nan=False) + "\n")
                    stream.flush()
                    write_paths(paths, record, cases[original["case_id"]])
                    if len(records) % 25 == 0:
                        paths.flush()
                        client.evaluations.flush()
                        checkpoint = dict(completed=len(records), planned=provenance["planned_trials"], elapsed_seconds=time.time() - started,
                                          latest={key: record[key] for key in ("consumer", "route", "trial_id", "status")},
                                          statuses=dict(Counter(item["status"] for item in records)), query_counts=dict(client.stats))
                        write_json(output / "checkpoint.json", checkpoint)
                        print(f"{len(records)}/{provenance['planned_trials']} {consumer} {route} {record['trial_id']} {record['status']}", flush=True)
        derivatives = jacobian_checks(client, scalar, coupled, cases, geometry)
        write_json(output / "jacobian_checks.json", derivatives)
        hybrids = hybrid_checks(client, scalar, coupled, cases)
        write_json(output / "hybrid_checks.json", hybrids)
        stats = dict(client.stats)
    verify_private_cache(output / "api")
    after = installed_inventory()
    frozen_api_unchanged = (sha256(provenance["driver_command"][0]) == provenance["driver_sha256"] and
                            sha256(CODE / "runtime_bridge/api_client.py") == provenance["client_sha256"] and
                            sha256(IDENTITY) == provenance["identity_sha256"])
    results = dict(completed=len(records), planned=provenance["planned_trials"], statuses=dict(Counter(record["status"] for record in records)),
                   elapsed_seconds=time.time() - started, query_counts=stats, installed_unchanged=before == after,
                   fd_statuses=dict(Counter(record["status"] for record in derivatives)),
                   fd_max_scaled_error=max(record["max_scaled_error"] for record in derivatives if record["status"] == "checked"),
                   hybrid_unchanged=all(record["scalar_unchanged"] and record["coupled_unchanged"] for record in hybrids),
                   maxrss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                   child_maxrss_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss, fallback_count=0,
                   finished=datetime.now().isoformat(), installed_after=after, frozen_api_unchanged=frozen_api_unchanged)
    write_json(output / "results.json", results)
    write_json(output / "checkpoint.json", dict(completed=len(records), planned=provenance["planned_trials"], status="finished"))
    print(json.dumps({key: results[key] for key in ("completed", "statuses", "elapsed_seconds", "fd_max_scaled_error", "hybrid_unchanged")}), flush=True)
    if before != after or not frozen_api_unchanged:
        raise RuntimeError("Installed EOS/cache inventory or frozen API changed")


if __name__ == "__main__":
    main()
