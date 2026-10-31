"""Matched longer-burning runs using the frozen TASK064 native stellar driver."""

import argparse
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time


CODE = Path(__file__).resolve().parent
WORKSPACE = CODE.parents[4]
ORIGINAL = WORKSPACE / "mesa/eos/eosFreeEOS_builder/biquintic/stellar_experiment"
EXPERIMENT = WORKSPACE / "freeeos_biquintic_progress/experiments/067_stellar_burning"
FROZEN_EXPERIMENT = WORKSPACE / "freeeos_biquintic_progress/experiments/064_stellar_native"
FAMILY = WORKSPACE / "freeeos_biquintic_progress/experiments/062_family_cache/family"
DRIVER = FROZEN_EXPERIMENT / "driver_build/star"
CLEARANCE = "065_NATIVE_API_PASS_20261030"
MAX_MODELS = 800
MAX_TIMESTEP_YEARS = 1_000_000
MAX_RUNTIME_SECONDS = 1200

sys.dont_write_bytecode = True
specification = importlib.util.spec_from_file_location("frozen_stellar_harness", ORIGINAL / "harness.py")
frozen = importlib.util.module_from_spec(specification)
specification.loader.exec_module(frozen)


def target(backend, pair):
    if not re.fullmatch(r"[A-Za-z0-9_]+", pair):
        raise ValueError("pair must contain only letters, digits or underscore")
    return EXPERIMENT / f"{backend}_{pair}"


def frozen_build():
    freeze = json.loads((FROZEN_EXPERIMENT / "build_freeze.json").read_text())
    if frozen.sha256(DRIVER) != freeze["executable_sha256"]:
        raise RuntimeError("frozen driver hash changed")
    if frozen.archives() != freeze["private_archives"]:
        raise RuntimeError("coherent private archive hashes changed")
    mapping = (FROZEN_EXPERIMENT / "driver_build/link.map").read_text()
    linked = [str(Path(line[5:].strip()).resolve()) for line in mapping.splitlines()
              if line.startswith("LOAD ") and line.endswith(".a")]
    for package, record in freeze["private_archives"].items():
        if [path for path in linked if Path(path).name == f"lib{package}.a"] != [record["path"]]:
            raise RuntimeError(f"frozen driver linked unexpected {package} archive")
    return freeze


def prepare(pair):
    execution = frozen.placement()
    freeze = frozen_build()
    if frozen.sha256(frozen.SEED / "seed.mod") != frozen.SEED_SHA:
        raise RuntimeError("original seed checksum changed")
    if frozen.sha256(FAMILY / "family.manifest") != "bd62aaaee40c82804b54423c89faabc4bfb560da70e7184f7659bade8cc75158":
        raise RuntimeError("frozen family manifest checksum changed")
    template = (ORIGINAL / "inlist.template").read_text()
    if template.count("max_model_number = 8") != 1 or template.count("max_years_for_timestep = 1000") != 1:
        raise RuntimeError("original inlist controls are not as expected")
    common = template.replace("@TABLE_DIR@", str(FAMILY))
    common = common.replace("max_model_number = 8", f"max_model_number = {MAX_MODELS}")
    common = common.replace("max_years_for_timestep = 1000", f"max_years_for_timestep = {MAX_TIMESTEP_YEARS}")
    inputs = {}
    for backend, switch in (("baseline", ".false."), ("candidate", ".true.")):
        run_dir = target(backend, pair)
        run_dir.mkdir(parents=True, exist_ok=False)
        for filename in ("seed.mod", "frozen20.net", "history_columns.list", "profile_columns.list"):
            shutil.copyfile(frozen.SEED / filename, run_dir / filename)
        shutil.copyfile(ORIGINAL / "burning21.net", run_dir / "burning21.net")
        inlist = common.replace("@BIQUINTIC@", switch)
        (run_dir / "inlist").write_text(inlist)
        inputs[backend] = inlist
        frozen.record(run_dir / "manifest.json", {
            "execution": execution, "backend": backend, "pair": pair,
            "max_model_number": MAX_MODELS, "max_years_for_timestep": MAX_TIMESTEP_YEARS,
            "driver_sha256": freeze["executable_sha256"],
            "family_manifest_sha256": frozen.sha256(FAMILY / "family.manifest"),
            "inputs": {path.name: frozen.sha256(path) for path in run_dir.iterdir() if path.is_file()},
        })
    if inputs["baseline"].replace("FreeEOS_use_biquintic_interpolation = .false.", "FreeEOS_use_biquintic_interpolation = SWITCH") != inputs["candidate"].replace("FreeEOS_use_biquintic_interpolation = .true.", "FreeEOS_use_biquintic_interpolation = SWITCH"):
        raise RuntimeError("paired inlists differ beyond native switch")
    frozen.record(EXPERIMENT / f"pair_{pair}.json", {
        "execution": execution, "clearance": CLEARANCE,
        "driver_sha256": freeze["executable_sha256"],
        "private_archives": freeze["private_archives"],
        "seed_sha256": frozen.SEED_SHA,
        "family_manifest_sha256": frozen.sha256(FAMILY / "family.manifest"),
        "max_model_number": MAX_MODELS, "max_years_for_timestep": MAX_TIMESTEP_YEARS,
        "initial_timestep_years": 100, "runtime_seconds_per_backend": MAX_RUNTIME_SECONDS,
        "inlist_difference": "FreeEOS_use_biquintic_interpolation only",
    })


def completion(run_dir):
    path = run_dir / "completion.txt"
    if not path.exists():
        return {}
    return dict(line.split(maxsplit=1) for line in path.read_text().splitlines() if line.strip())


def run(backend, pair):
    execution = frozen.placement()
    freeze = frozen_build()
    run_dir = target(backend, pair)
    manifest = json.loads((run_dir / "manifest.json").read_text())
    if (run_dir / "command.json").exists():
        raise RuntimeError("attempt already exists; preserve its evidence")
    if manifest["driver_sha256"] != freeze["executable_sha256"]:
        raise RuntimeError("prepared driver checksum mismatch")
    for filename, expected in manifest["inputs"].items():
        if frozen.sha256(run_dir / filename) != expected:
            raise RuntimeError(f"prepared input changed: {filename}")
    if frozen.sha256(FAMILY / "family.manifest") != manifest["family_manifest_sha256"]:
        raise RuntimeError("family manifest changed")
    if backend == "candidate":
        baseline_dir = target("baseline", pair)
        baseline_outcome = json.loads((baseline_dir / "outcome.json").read_text())
        baseline_completion = completion(baseline_dir)
        if (baseline_outcome["returncode"] != 0 or baseline_outcome["timed_out"]
                or baseline_completion.get("model_number") != str(MAX_MODELS)
                or baseline_completion.get("native_failures") != "0"):
            raise RuntimeError("baseline did not pass; preserve and diagnose before candidate")
    installed_before = frozen.installed_inventory()
    frozen.record(run_dir / "installed_before.json", installed_before)
    cache_dir = run_dir / "runtime_cache/eosDT_cache"
    cache_dir.mkdir(parents=True)
    for source in (frozen.MESA / "data/eosDT_data/cache").glob("*.bin"):
        shutil.copy2(source, cache_dir / source.name)
    temp_dir = run_dir / "runtime_temp_cache"
    temp_dir.mkdir()
    runtime_env = frozen.environment()
    runtime_env.update(MESA_CACHES_DIR=str(cache_dir.parent), MESA_TEMP_CACHES_DIR=str(temp_dir))
    timeout = min(MAX_RUNTIME_SECONDS, int(frozen.STOP_TIME - time.time() - 90))
    if timeout < 30:
        raise RuntimeError("insufficient time before numerical stop")
    command = [str(DRIVER)]
    frozen.record(run_dir / "command.json", {
        "execution": execution, "command": command, "cwd": str(run_dir),
        "backend": backend, "pair": pair, "clearance": CLEARANCE,
        "threads": 1, "address_space_bytes": 16 * 1024**3,
        "timeout_seconds": timeout,
    })
    started = time.monotonic()
    timed_out = False
    with (run_dir / "run.log").open("w") as stream:
        process = subprocess.Popen(command, cwd=run_dir, env=runtime_env,
                                   stdout=stream, stderr=subprocess.STDOUT,
                                   preexec_fn=frozen.memory_limit)
        frozen.record(run_dir / "process.json", {"pid": process.pid, "started": time.time()})
        try:
            returncode = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            process.kill()
            returncode = process.wait()
    frozen.record(run_dir / "outcome.json", {
        "returncode": returncode, "timed_out": timed_out,
        "seconds": time.monotonic() - started,
        "completion_written": (run_dir / "completion.txt").exists(),
    })
    installed_after = frozen.installed_inventory()
    frozen.record(run_dir / "installed_after.json", installed_after)
    if installed_before != installed_after:
        raise RuntimeError("installed libraries/modules/EOS caches changed during run")
    if int(completion(run_dir).get("native_failures", "0")) != 0:
        raise RuntimeError("native EOS failure; stop and diagnose")


def analyze(pair):
    results = {"pair": pair, "backends": {}}
    for backend in ("baseline", "candidate"):
        run_dir = target(backend, pair)
        if not (run_dir / "outcome.json").exists():
            continue
        events = (run_dir / "events.tsv").read_text().splitlines() if (run_dir / "events.tsv").exists() else []
        accepted = [line.split() for line in events if line.strip().startswith("accepted")]
        startup = next((line.split() for line in events if line.strip().startswith("startup")), None)
        iterations = (run_dir / "iterations.tsv").read_text().splitlines() if (run_dir / "iterations.tsv").exists() else []
        iteration_numbers = [int(line.split()[2]) for line in iterations if line and not line.startswith("#")]
        state = completion(run_dir)
        outcome = json.loads((run_dir / "outcome.json").read_text())
        results["backends"][backend] = {
            "outcome": outcome, "completion": state,
            "accepted_models": len(accepted),
            "initial_center_h1": float(startup[10]) if startup else None,
            "last_accepted_age_yr": float(accepted[-1][4]) if accepted else None,
            "last_accepted_center_h1": float(accepted[-1][10]) if accepted else None,
            "center_h1_depletion": float(startup[10]) - float(accepted[-1][10]) if startup and accepted else None,
            "first_model_with_h1_depletion_at_least_1e_4": next(
                (int(row[1]) for row in accepted if startup and float(startup[10]) - float(row[10]) >= 1e-4), None),
            "last_accepted_center_he4": float(accepted[-1][11]) if accepted else None,
            "last_accepted_center_he3": float(accepted[-1][12]) if accepted else None,
            "accepted_age_yr": [float(row[4]) for row in accepted],
            "accepted_dt_s": [float(row[5]) for row in accepted],
            "maximum_retry_index": max((int(row[3]) for row in accepted), default=None),
            "maximum_newton_iteration": max(iteration_numbers, default=None),
            "eos_failure_log_lines": [line for line in (run_dir / "run.log").read_text().splitlines()
                                      if re.search(r"-7950|FreeEOS.*(error|fail)|eos.*ierr", line, re.I)][:20],
        }
    if len(results["backends"]) == 2:
        baseline = results["backends"]["baseline"]
        candidate = results["backends"]["candidate"]
        both = (baseline, candidate)
        results["checks"] = {
            "both_completed": all(item["outcome"]["returncode"] == 0
                                  and item["completion"].get("model_number") == str(MAX_MODELS)
                                  and item["completion"].get("termination_code") == "4"
                                  and item["accepted_models"] == MAX_MODELS for item in both),
            "burning_mixing_abundance_active": all(
                item["completion"].get(key) == "T" for item in both
                for key in ("saw_burn", "saw_mix", "saw_abundance_equations")),
            "native_only_candidate": baseline["completion"].get("native_successes") == "0"
                                     and int(candidate["completion"].get("native_successes", "0")) > 0,
            "native_zero_failures": all(item["completion"].get("native_failures") == "0" for item in both),
            "both_reach_h1_depletion_target": all(item["center_h1_depletion"] is not None
                                                   and item["center_h1_depletion"] >= 1e-4 for item in both),
        }
        if baseline["accepted_models"] == candidate["accepted_models"] and baseline["accepted_models"]:
            results["paired"] = {
                "max_abs_age_yr": max(abs(left - right) for left, right in zip(
                    baseline["accepted_age_yr"], candidate["accepted_age_yr"])),
                "max_abs_dt_s": max(abs(left - right) for left, right in zip(
                    baseline["accepted_dt_s"], candidate["accepted_dt_s"])),
                "final_center_h1_difference": candidate["last_accepted_center_h1"] - baseline["last_accepted_center_h1"],
            }
    frozen.record(EXPERIMENT / f"analysis_{pair}.json", results)
    print(json.dumps({"checks": results.get("checks"), "backends": {
        backend: {"accepted_models": item["accepted_models"],
                  "last_accepted_center_h1": item["last_accepted_center_h1"],
                  "outcome": item["outcome"]}
        for backend, item in results["backends"].items()}}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    commands.add_parser("freeze-check")
    preparation = commands.add_parser("prepare")
    preparation.add_argument("--pair", required=True)
    execution = commands.add_parser("run")
    execution.add_argument("backend", choices=("baseline", "candidate"))
    execution.add_argument("--pair", required=True)
    analysis = commands.add_parser("analyze")
    analysis.add_argument("--pair", required=True)
    args = parser.parse_args()
    if args.action == "freeze-check":
        print(json.dumps(frozen_build(), indent=2))
    elif args.action == "prepare":
        prepare(args.pair)
    elif args.action == "run":
        run(args.backend, args.pair)
    else:
        analyze(args.pair)


if __name__ == "__main__":
    main()
