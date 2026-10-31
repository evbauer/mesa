"""Run TASK069 matched-age native and legacy stellar models."""

import argparse
import importlib.util
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import time


CODE = Path(__file__).resolve().parent
WORKSPACE = CODE.parents[4]
ORIGINAL = WORKSPACE / "mesa/eos/eosFreeEOS_builder/biquintic/stellar_experiment"
EXPERIMENT = WORKSPACE / "freeeos_biquintic_progress/experiments/069_native_timestep"
FROZEN = WORKSPACE / "freeeos_biquintic_progress/experiments/064_stellar_native"
FAMILY = WORKSPACE / "freeeos_biquintic_progress/experiments/062_family_cache/family"
DRIVER = FROZEN / "driver_build/star"
STOP_TIME = 1793424600
CAPS = {"coarse": 1_000_000, "half": 500_000}
BACKENDS = {"legacy": ".false.", "native": ".true."}
MAX_MODELS = 2500
MAX_AGE = 600_000_000
MAX_SECONDS = 1800

sys.dont_write_bytecode = True
specification = importlib.util.spec_from_file_location("original_stellar_harness", ORIGINAL / "harness.py")
frozen = importlib.util.module_from_spec(specification)
specification.loader.exec_module(frozen)


def placement():
    hostname = socket.gethostname().split(".")[0]
    job = os.environ.get("SLURM_JOB_ID")
    if (hostname, job) != ("rzwhippet4", "5532121"):
        raise RuntimeError(f"TASK069 placement mismatch: {hostname}, {job}")
    if time.time() >= STOP_TIME - 90:
        raise RuntimeError("TASK069 numerical deadline is too near")
    return {"hostname": hostname, "job": job, "timestamp": time.time()}


def identities():
    freeze = json.loads((FROZEN / "build_freeze.json").read_text())
    if frozen.sha256(DRIVER) != freeze["executable_sha256"]:
        raise RuntimeError("frozen driver changed")
    if frozen.archives() != freeze["private_archives"]:
        raise RuntimeError("private archives changed")
    if frozen.sha256(frozen.SEED / "seed.mod") != frozen.SEED_SHA:
        raise RuntimeError("original seed changed")
    family_hash = frozen.sha256(FAMILY / "family.manifest")
    if family_hash != "bd62aaaee40c82804b54423c89faabc4bfb560da70e7184f7659bade8cc75158":
        raise RuntimeError("family manifest changed")
    runtime_hash = frozen.sha256(CODE.parent / "family_runtime/freeeos_family.f90")
    if runtime_hash != "03b32156c9113cdbf63f789cfe7298f2f1862999bcd024cd99704024a5bfa70f":
        raise RuntimeError("native runtime source changed")
    for filename, expected in freeze["sources"].items():
        if filename in ("Makefile", "inlist.template", "run.f90", "run_star_extras.f90", "burning21.net"):
            path = ORIGINAL / filename
            if frozen.sha256(path) != expected:
                raise RuntimeError(f"frozen source changed: {filename}")
    linked = [str(Path(line[5:].strip()).resolve())
              for line in (FROZEN / "driver_build/link.map").read_text().splitlines()
              if line.startswith("LOAD ") and line.endswith(".a")]
    for package, record in freeze["private_archives"].items():
        if [path for path in linked if Path(path).name == f"lib{package}.a"] != [record["path"]]:
            raise RuntimeError(f"driver linked unexpected {package} archive")
    return {"driver_sha256": freeze["executable_sha256"],
            "private_archives": freeze["private_archives"],
            "source_sha256": freeze["sources"],
            "current_original_harness_sha256": frozen.sha256(ORIGINAL / "harness.py"),
            "seed_sha256": frozen.SEED_SHA,
            "family_manifest_sha256": family_hash,
            "runtime_source_sha256": runtime_hash,
            "installed": frozen.installed_inventory()}


def directory(backend, cap):
    return EXPERIMENT / f"{backend}_{cap}"


def prepare():
    execution = placement()
    identity = identities()
    template = (ORIGINAL / "inlist.template").read_text()
    assert template.count("max_model_number = 8") == 1
    assert template.count("max_years_for_timestep = 1000") == 1
    common = template.replace("@TABLE_DIR@", str(FAMILY))
    common = common.replace("max_model_number = 8", f"max_model_number = {MAX_MODELS}")
    common = common.replace("   max_years_for_timestep = 1000", "   max_age = 600000000\n   max_num_profile_models = 2600\n   max_years_for_timestep = @CAP@")
    inputs = {}
    for backend, switch in BACKENDS.items():
        for cap, years in CAPS.items():
            run_dir = directory(backend, cap)
            run_dir.mkdir(parents=True, exist_ok=False)
            for filename in ("seed.mod", "frozen20.net", "history_columns.list", "profile_columns.list"):
                shutil.copyfile(frozen.SEED / filename, run_dir / filename)
            shutil.copyfile(ORIGINAL / "burning21.net", run_dir / "burning21.net")
            inlist = common.replace("@BIQUINTIC@", switch).replace("@CAP@", str(years))
            (run_dir / "inlist").write_text(inlist)
            inputs[f"{backend}_{cap}"] = inlist
            frozen.record(run_dir / "manifest.json", {
                "execution": execution, "backend": backend, "cap": cap,
                "clearance": "065_NATIVE_API_PASS_20261030",
                "identity": {key: value for key, value in identity.items() if key != "installed"},
                "inputs": {path.name: frozen.sha256(path) for path in run_dir.iterdir() if path.is_file()},
            })
    normalized = {name: value.replace(".false.", "@SWITCH@").replace(".true.", "@SWITCH@").replace("1000000", "@CAP@").replace("500000", "@CAP@")
                  for name, value in inputs.items()}
    if len(set(normalized.values())) != 1:
        raise RuntimeError("inlists differ beyond backend and cap")
    frozen.record(EXPERIMENT / "freeze_before.json", identity)
    frozen.record(EXPERIMENT / "design.json", {
        "execution": execution, "deadline_epoch": STOP_TIME,
        "max_age_years": MAX_AGE, "max_model_number": MAX_MODELS,
        "initial_timestep_years": 100, "max_num_profile_models": 2600,
        "caps_years": CAPS, "run_wall_limit_seconds": MAX_SECONDS,
        "common_changes_from_original": ["max_age", "max_model_number", "max_num_profile_models", "max_years_for_timestep"],
        "between_backend_difference": "FreeEOS_use_biquintic_interpolation only",
        "between_cap_difference": "max_years_for_timestep only"})


def run(backend, cap):
    execution = placement()
    before = identities()
    run_dir = directory(backend, cap)
    manifest = json.loads((run_dir / "manifest.json").read_text())
    if (run_dir / "command.json").exists():
        raise RuntimeError("attempt already exists")
    if {key: before[key] for key in manifest["identity"]} != manifest["identity"]:
        raise RuntimeError("frozen identity differs from preparation")
    for filename, expected in manifest["inputs"].items():
        if frozen.sha256(run_dir / filename) != expected:
            raise RuntimeError(f"prepared input changed: {filename}")
    frozen.record(run_dir / "installed_before.json", before["installed"])
    cache_dir = run_dir / "runtime_cache/eosDT_cache"
    cache_dir.mkdir(parents=True)
    for source in (frozen.MESA / "data/eosDT_data/cache").glob("*.bin"):
        shutil.copy2(source, cache_dir / source.name)
    temp_dir = run_dir / "runtime_temp_cache"
    temp_dir.mkdir()
    environment = frozen.environment()
    environment.update(MESA_CACHES_DIR=str(cache_dir.parent), MESA_TEMP_CACHES_DIR=str(temp_dir))
    timeout = min(MAX_SECONDS, int(STOP_TIME - time.time() - 60))
    if timeout < 30:
        raise RuntimeError("insufficient time before numerical stop")
    command = [str(DRIVER)]
    frozen.record(run_dir / "command.json", {"execution": execution, "command": command,
        "cwd": str(run_dir), "timeout_seconds": timeout, "threads": 1,
        "address_space_bytes": 16 * 1024**3})
    started = time.monotonic()
    timed_out = False
    with (run_dir / "run.log").open("w") as stream:
        process = subprocess.Popen(command, cwd=run_dir, env=environment,
                                   stdout=stream, stderr=subprocess.STDOUT,
                                   preexec_fn=frozen.memory_limit)
        frozen.record(run_dir / "process.json", {"pid": process.pid, "started": time.time()})
        try:
            returncode = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            process.kill()
            returncode = process.wait()
    frozen.record(run_dir / "outcome.json", {"returncode": returncode,
        "timed_out": timed_out, "seconds": time.monotonic() - started,
        "completion_written": (run_dir / "completion.txt").exists()})
    after = identities()
    frozen.record(run_dir / "installed_after.json", after["installed"])
    if before != after:
        raise RuntimeError("frozen or installed identities changed during run")
    if returncode or timed_out:
        raise RuntimeError(f"model failed: returncode={returncode}, timed_out={timed_out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run", "freeze-check", "finish"))
    parser.add_argument("backend", nargs="?", choices=BACKENDS)
    parser.add_argument("cap", nargs="?", choices=CAPS)
    args = parser.parse_args()
    if args.action == "prepare":
        prepare()
    elif args.action == "freeze-check":
        print(json.dumps(identities(), indent=2))
    elif args.action == "finish":
        after = identities()
        before = json.loads((EXPERIMENT / "freeze_before.json").read_text())
        if any(after[key] != value for key, value in before.items()):
            raise RuntimeError("frozen identities changed since preparation")
        frozen.record(EXPERIMENT / "freeze_after.json", after)
    elif args.backend and args.cap:
        run(args.backend, args.cap)
    else:
        parser.error("run requires backend and cap")


if __name__ == "__main__":
    main()
