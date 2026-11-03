"""Create the 070 survey in explicit, restartable stages on the allocated node."""

import argparse
import json
import os
import resource
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np

from common import (BIQUINTIC, CODE, HELM048, JOB, MASSFILE, OLD, OUTPUT, SOURCE,
                    WORKSPACE, digest, placement, save_json, scaled_jump, summary)

sys.path.insert(0, str(BIQUINTIC))
sys.path.insert(0, str(HELM048))
from freeeos_potential import (LN10, construct_twists, evaluate, read_table,
                               source_constraints, thermodynamic_partials)


STOP = datetime(2026, 11, 3, 0, 25, tzinfo=ZoneInfo("America/Los_Angeles"))
PROBE = BIQUINTIC / "mesa_probe/run_probe.py"
BATCH = 4096


def time_left():
    return max(0, int((STOP - datetime.now(STOP.tzinfo)).total_seconds()))


def environment(directory):
    sdk = WORKSPACE / "mesasdk-x86_64-linux-26.6.1"
    result = dict(os.environ)
    result.update(MESA_DIR=str(WORKSPACE / "mesa"), MESASDK_ROOT=str(sdk),
                  MESA_CACHES_DIR=str(directory / "runtime_cache"),
                  MESA_TEMP_CACHES_DIR=str(directory / "runtime_temp_cache"),
                  LD_LIBRARY_PATH=f"{sdk}/lib:{sdk}/lib64:" + result.get("LD_LIBRARY_PATH", ""),
                  PYTHONDONTWRITEBYTECODE="1")
    result.pop("MESA_TEMP_CACHES_DISABLE", None)
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        result[key] = "1"
    private_cache = directory / "runtime_cache/eosDT_cache"
    private_cache.mkdir(parents=True, exist_ok=True)
    (directory / "runtime_temp_cache").mkdir(parents=True, exist_ok=True)
    installed = WORKSPACE / "mesa/data/eosDT_data/cache/helm_table.bin"
    link = private_cache / installed.name
    if not link.exists() and not link.is_symlink():
        link.symlink_to(installed)
    if not link.is_symlink() or link.resolve() != installed.resolve():
        raise ValueError("private HELM cache link differs from installed cache")
    return result


def run_logged(command, stdin_path, stdout_path, stderr_path, cwd, env):
    remaining = time_left()
    if remaining < 60:
        raise RuntimeError("TASK070 numerical stop time reached")
    with stdin_path.open() if stdin_path else open(os.devnull) as source, \
         stdout_path.open("x") as target, stderr_path.open("x") as errors:
        completed = subprocess.run(command, cwd=cwd, env=env, stdin=source, stdout=target,
                                   stderr=errors, timeout=remaining, check=False)
    save_json(stdout_path.with_suffix(".command.json"), {"command": command,
              "returncode": completed.returncode, "job": JOB})
    if completed.returncode:
        raise RuntimeError(f"command returned {completed.returncode}: {command}")


def initialize():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    binary = OUTPUT / "build/helm_map_probe"
    if not binary.exists():
        remaining = time_left()
        if remaining < 60:
            raise RuntimeError("TASK070 numerical stop time reached")
        with (OUTPUT / "build.log").open("x") as log:
            subprocess.run(["make", "-f", str(CODE / "Makefile"), f"BUILD_DIR={OUTPUT / 'build'}"],
                           cwd=WORKSPACE, stdout=log, stderr=subprocess.STDOUT,
                           timeout=remaining, check=True, env=environment(OUTPUT))
    if not (OUTPUT / "temperature_axis.txt").exists():
        run_logged([str(binary), str(WORKSPACE / "mesa"), str(MASSFILE), str(OUTPUT)],
                   None, OUTPUT / "helm_initialize.tsv", OUTPUT / "helm_initialize.stderr",
                   OUTPUT, environment(OUTPUT))
    save_json(OUTPUT / "input_hashes.json", {"source": digest(SOURCE), "massfile": digest(MASSFILE),
              "helm_cache": digest(WORKSPACE / "mesa/data/eosDT_data/cache/helm_table.bin"),
              "old_manifest": digest(OLD / "face_manifest.npz"),
              "old_results": digest(OLD / "face_results.npz"),
              "old_probe": digest(OLD / "face_probe.tsv"),
              "new_driver": digest(CODE / "helm_probe.f90"),
              "new_binary": digest(binary), "new_script": digest(__file__)})


def old_bounds():
    with np.load(OLD / "face_manifest.npz") as old:
        coordinates = [(old[f"{direction}_logQ"] + 2 * old[f"{direction}_logT"] - 12,
                        old[f"{direction}_logT"]) for direction in ("u", "tau")]
    density_limits = [float(min(np.min(item[0]) for item in coordinates)),
                      float(max(np.max(item[0]) for item in coordinates))]
    temperature_limits = [float(min(np.min(item[1]) for item in coordinates)),
                          float(max(np.max(item[1]) for item in coordinates))]
    return {"old_logRho": density_limits, "old_logT": temperature_limits,
            "new_logRho": [density_limits[0], density_limits[1] + 2],
            "new_logT": [temperature_limits[0], temperature_limits[1] + 2]}


def grid(direction, density_axis, temperature_axis, ye, bounds):
    if direction == "T":
        density_indices = np.arange(len(density_axis) - 1)
        temperature_indices = np.arange(1, len(temperature_axis) - 1)
        density = (density_axis[:-1].astype(np.longdouble) + density_axis[1:].astype(np.longdouble)) / 2
        temperature = temperature_axis[temperature_indices]
    else:
        density_indices = np.arange(1, len(density_axis) - 1)
        temperature_indices = np.arange(len(temperature_axis) - 1)
        density = density_axis[density_indices]
        temperature = (temperature_axis[:-1].astype(np.longdouble) +
                       temperature_axis[1:].astype(np.longdouble)) / 2
    density = np.asarray(density, dtype=np.float64)
    temperature = np.asarray(temperature, dtype=np.float64)
    density_keep = (np.log10(density / ye) >= bounds["new_logRho"][0]) & (
        np.log10(density / ye) <= bounds["new_logRho"][1])
    temperature_keep = (np.log10(temperature) >= bounds["new_logT"][0]) & (
        np.log10(temperature) <= bounds["new_logT"][1])
    return density_indices[density_keep], temperature_indices[temperature_keep], density[density_keep], temperature[temperature_keep]


def prepare():
    composition = np.loadtxt(OUTPUT / "composition.tsv", skiprows=1)
    ye = float(composition[5])
    density_axis = np.loadtxt(OUTPUT / "density_axis.txt")
    temperature_axis = np.loadtxt(OUTPUT / "temperature_axis.txt")
    bounds = old_bounds()
    bounds.update({"helm_native_logD": [float(np.log10(density_axis[0])), float(np.log10(density_axis[-1]))],
                   "helm_native_logT": [float(np.log10(temperature_axis[0])), float(np.log10(temperature_axis[-1]))],
                   "Ye": ye})
    save_json(OUTPUT / "bounds.json", bounds)
    for direction in ("T", "d"):
        density_indices, temperature_indices, density, temperature = grid(
            direction, density_axis, temperature_axis, ye, bounds)
        manifest = {"density_index": density_indices, "temperature_index": temperature_indices,
                    "din": density, "temperature": temperature}
        np.savez_compressed(OUTPUT / f"helm_{direction}_grid.npz", **manifest)
        if (OUTPUT / f"helm_{direction}_queries.txt").exists():
            raise FileExistsError("Refusing to overwrite existing native queries")
        with (OUTPUT / f"helm_{direction}_queries.txt").open("x") as stream:
            for din in density:
                rho = din / ye
                for temp in temperature:
                    stream.write(f"{rho:.17e} {temp:.17e}\n")
        print(f"{direction}: {len(density)} x {len(temperature)} = {len(density)*len(temperature):,}", flush=True)


def prepare_routes():
    ye = float(np.loadtxt(OUTPUT / "composition.tsv", skiprows=1)[5])
    for direction in ("T", "d"):
        with np.load(OUTPUT / f"helm_{direction}_grid.npz") as grid, \
             (OUTPUT / f"helm_{direction}_route_queries.txt").open("x") as stream:
            for din in grid["din"]:
                log_density = np.log10(din / ye)
                for temperature in grid["temperature"]:
                    stream.write(f"{log_density:.17e} {np.log10(temperature):.17e}\n")


def route(direction):
    route_dir = OUTPUT / f"routing_{direction}_corrected"
    route_dir.mkdir(exist_ok=True)
    command = [sys.executable, "-B", str(PROBE), "--mesa-dir", str(WORKSPACE / "mesa"),
               "--output", str(route_dir)]
    run_logged(command, OUTPUT / f"helm_{direction}_route_queries.txt", route_dir / "probe.tsv",
               route_dir / "probe.stderr", route_dir, environment(route_dir))


def kernel(direction):
    command = [str(OUTPUT / "build/helm_map_probe"), str(WORKSPACE / "mesa"), str(MASSFILE), str(OUTPUT)]
    run_logged(command, OUTPUT / f"helm_{direction}_queries.txt", OUTPUT / f"helm_{direction}_kernel.tsv",
               OUTPUT / f"helm_{direction}_kernel.stderr", OUTPUT, environment(OUTPUT))


def freeeos():
    _, logq, logt, rows = read_table(SOURCE)
    low, source_valid, _ = source_constraints(rows, logq, logt, split_radiation=False,
                                              minimum_logT=3.2)
    nodes, complete, _ = construct_twists(low, source_valid, (logq[1] - logq[0]) * LN10,
                                          (logt[1] - logt[0]) * LN10, radius=2, route="mean")
    with np.load(OLD / "face_manifest.npz") as manifest, np.load(OLD / "face_results.npz") as old:
        results = {}
        checks = {}
        for axis, direction in ((0, "u"), (1, "tau")):
            face_q = manifest[f"{direction}_logQ"]
            face_t = manifest[f"{direction}_logT"]
            density_index = manifest[f"{direction}_density_index"]
            temperature_index = manifest[f"{direction}_temperature_index"]
            pure = old[f"{direction}_pure"]
            cause = old[f"{direction}_support_code"]
            chosen = np.flatnonzero(pure & (cause == 0))
            trace = {name: np.full(len(face_q), np.nan) for name in
                     ("Cv_left", "Cv_right", "Cv_t_left", "Cv_t_right", "denominator", "jump")}
            for start in range(0, len(chosen), BATCH):
                positions = chosen[start:start + BATCH]
                left_cell = density_index[positions], temperature_index[positions]
                right_cell = (density_index[positions] + (axis == 0),
                              temperature_index[positions] + (axis == 1))
                left = thermodynamic_partials(evaluate(nodes, complete, logq, logt,
                    face_q[positions], face_t[positions], cell_indices=left_cell),
                    face_q[positions], face_t[positions])
                right = thermodynamic_partials(evaluate(nodes, complete, logq, logt,
                    face_q[positions], face_t[positions], cell_indices=right_cell),
                    face_q[positions], face_t[positions])
                jump, denominator = scaled_jump(left["Cv"], right["Cv"], left["Cv_t"], right["Cv_t"])
                for side, values in (("left", left), ("right", right)):
                    trace[f"Cv_{side}"][positions] = values["Cv"]
                    trace[f"Cv_t_{side}"][positions] = values["Cv_t"]
                trace["denominator"][positions] = denominator
                trace["jump"][positions] = jump
                if time_left() < 60:
                    raise RuntimeError("TASK070 numerical stop time reached")
            archived = old[f"{direction}_Cv_t_jump"][chosen]
            current_old = (trace["Cv_t_left"][chosen] - trace["Cv_t_right"][chosen]) / np.maximum(1, trace["denominator"][chosen])
            difference = np.abs(current_old - archived)
            checks[direction] = {"selected": int(len(chosen)), "old_max_absolute_difference": float(np.max(difference)),
                                 "old_changed_count_over_1e-12": int(np.count_nonzero(difference > 1e-12)),
                                 "denominator_below_one": int(np.count_nonzero(trace["denominator"][chosen] < 1)),
                                 "new_vs_old_changed_over_1e-15": int(np.count_nonzero(np.abs(trace["jump"][chosen] - archived) > 1e-15)),
                                 "minimum_denominator": float(np.min(trace["denominator"][chosen])),
                                 "new_jump": summary(trace["jump"][chosen])}
            if checks[direction]["old_changed_count_over_1e-12"]:
                raise ValueError("old060 exact-face reproduction failed")
            results.update({f"{direction}_{key}": value for key, value in
                            {"logQ": face_q, "logT": face_t, "pure": pure,
                             "support_code": cause, **trace}.items()})
    np.savez_compressed(OUTPUT / "freeeos_traces.npz", **results)
    save_json(OUTPUT / "freeeos_summary.json", checks)
    print(json.dumps(checks, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("init", "prepare", "prepare_routes", "route_T", "route_d", "kernel", "kernel_d", "freeeos"))
    args = parser.parse_args()
    placement()
    resource.setrlimit(resource.RLIMIT_AS, (12 * 1024**3, 12 * 1024**3))
    if args.stage == "init":
        initialize()
    elif args.stage == "prepare":
        prepare()
    elif args.stage == "prepare_routes":
        prepare_routes()
    elif args.stage.startswith("route_"):
        route(args.stage.removeprefix("route_"))
    elif args.stage == "kernel":
        kernel("T")
    elif args.stage == "kernel_d":
        kernel("d")
    else:
        freeeos()


if __name__ == "__main__":
    main()
