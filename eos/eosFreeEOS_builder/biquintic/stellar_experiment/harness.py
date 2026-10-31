"""Bounded ordinary MESA stellar comparison for native FreeEOS interpolation."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import shutil
import socket
import subprocess
import time


CODE = Path(__file__).resolve().parent
WORKSPACE = CODE.parents[4]
MESA = WORKSPACE / "mesa"
SDK = WORKSPACE / "mesasdk-x86_64-linux-26.6.1"
EXPERIMENT = WORKSPACE / "freeeos_biquintic_progress/experiments/064_stellar_native"
PRIVATE = WORKSPACE / "freeeos_biquintic_progress/experiments/065_native_integration/build"
SEED = WORKSPACE / "freeeos_biquintic_progress/experiments/049_stellar_preparation/control_legacy"
SEED_SHA = "cd5b74e7ccd55d8253d045abc5bb95cba1506ba42fba70482bccf5f29cf9eb3e"
PACKAGES = ("eos", "rates", "kap", "net", "atm", "star_data", "star")
STOP_TIME = 1793410431


def suffix(identifier):
    if identifier and not re.fullmatch(r"[a-zA-Z0-9_]+", identifier):
        raise ValueError("pair/build identifier must contain only letters, digits or underscore")
    return f"_{identifier}" if identifier else ""


def build_path(build_id):
    return EXPERIMENT / f"driver_build{suffix(build_id)}"


def run_path(backend, pair_id):
    return EXPERIMENT / f"{backend}{suffix(pair_id)}"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def record(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def placement():
    hostname = socket.gethostname().split(".")[0]
    job = os.environ.get("SLURM_JOB_ID")
    if hostname != "rzwhippet13" or job != "5530692":
        raise RuntimeError(f"compute placement mismatch: {hostname} job {job}")
    if time.time() >= STOP_TIME - 120:
        raise RuntimeError("TASK064 numerical stop is too near")
    return {"hostname": hostname, "job": job, "timestamp": time.time()}


def environment():
    result = dict(os.environ)
    result.update(MESA_DIR=str(MESA), MESASDK_ROOT=str(SDK), OMP_NUM_THREADS="1",
                  OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1",
                  BLIS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    result["PATH"] = str(SDK / "bin") + ":" + result["PATH"]
    result["LD_LIBRARY_PATH"] = f"{SDK}/lib:{SDK}/lib64:" + result.get("LD_LIBRARY_PATH", "")
    return result


def memory_limit():
    resource.setrlimit(resource.RLIMIT_AS, (16 * 1024**3, 16 * 1024**3))


def archives():
    found = {}
    for package in PACKAGES:
        path = PRIVATE / package / "lib" / f"lib{package}.a"
        if not path.is_file():
            raise RuntimeError(f"coherent private archive absent: {path}")
        found[package] = {"path": str(path.resolve()), "sha256": sha256(path)}
    return found


def installed_inventory():
    paths = list((MESA / "build").glob("*/lib/*.a"))
    paths += list((MESA / "build").glob("*/include/*.mod"))
    paths += list((MESA / "data/eosDT_data/cache").glob("*.bin"))
    return {str(path.relative_to(MESA)): sha256(path) for path in sorted(paths)}


def build(build_id):
    execution = placement()
    dependency_record = archives()
    EXPERIMENT.mkdir(parents=True, exist_ok=True)
    driver_build = build_path(build_id)
    name_suffix = suffix(build_id)
    if driver_build.exists():
        raise RuntimeError("driver build already exists; preserve it")
    command = ["make", "-C", str(CODE), "-j4", f"BUILD={driver_build}"]
    record(EXPERIMENT / f"build_command{name_suffix}.json", {"execution": execution, "command": command,
           "private_archives": dependency_record})
    with (EXPERIMENT / f"build{name_suffix}.log").open("w") as stream:
        process = subprocess.run(command, env=environment(), stdout=stream,
                                 stderr=subprocess.STDOUT, timeout=600)
    if process.returncode:
        raise RuntimeError(f"driver build failed; see build{name_suffix}.log")
    if archives() != dependency_record:
        raise RuntimeError("private archives changed during driver build")
    mapping = (driver_build / "link.map").read_text()
    loaded = [Path(line[5:].strip()).resolve() for line in mapping.splitlines()
              if line.startswith("LOAD ") and line.endswith(".a")]
    for package, expected in dependency_record.items():
        matches = [str(path) for path in loaded if path.name == f"lib{package}.a"]
        if matches != [expected["path"]]:
            raise RuntimeError(f"driver linked unexpected {package}: {matches}")
    record(EXPERIMENT / f"build_freeze{name_suffix}.json", {"execution": execution,
           "executable_sha256": sha256(driver_build / "star"), "private_archives": dependency_record,
           "sources": {path.name: sha256(path) for path in CODE.iterdir() if path.is_file()}})


def prepare(table_dir, pair_id, max_model_number):
    execution = placement()
    if max_model_number not in (8, 100):
        raise ValueError("only eight or 100 accepted models are authorized")
    if sha256(SEED / "seed.mod") != SEED_SHA:
        raise RuntimeError("source seed checksum mismatch")
    table_dir = Path(table_dir).resolve()
    if not (table_dir / "family.manifest").is_file():
        raise RuntimeError(f"family manifest absent: {table_dir}")
    template = (CODE / "inlist.template").read_text()
    if "@BIQUINTIC@" not in template or "@TABLE_DIR@" not in template:
        raise RuntimeError("inlist template markers absent")
    if template.count("max_model_number = 8") != 1:
        raise RuntimeError("bounded model-number control missing")
    if "'" in str(table_dir):
        raise RuntimeError("table directory cannot contain a quote")
    for backend, value in (("baseline", ".false."), ("candidate", ".true.")):
        target = run_path(backend, pair_id)
        target.mkdir(parents=True, exist_ok=False)
        for filename in ("seed.mod", "frozen20.net", "history_columns.list", "profile_columns.list"):
            shutil.copyfile(SEED / filename, target / filename)
        shutil.copyfile(CODE / "burning21.net", target / "burning21.net")
        inlist = template.replace("@BIQUINTIC@", value).replace("@TABLE_DIR@", str(table_dir))
        inlist = inlist.replace("max_model_number = 8", f"max_model_number = {max_model_number}")
        (target / "inlist").write_text(inlist)
        record(target / "manifest.json", {"execution": execution, "backend": backend,
               "pair_id": pair_id, "max_model_number": max_model_number,
               "table_dir": str(table_dir), "family_manifest_sha256": sha256(table_dir / "family.manifest"),
               "inputs": {path.name: sha256(path) for path in target.iterdir() if path.is_file()}})


def run(backend, clearance, pair_id, build_id):
    execution = placement()
    if backend == "candidate" and not clearance:
        raise RuntimeError("candidate needs explicit lead clearance record")
    target = run_path(backend, pair_id)
    manifest = json.loads((target / "manifest.json").read_text())
    if (target / "command.json").exists():
        raise RuntimeError("attempt already exists; preserve its artifacts")
    for filename, expected in manifest["inputs"].items():
        if sha256(target / filename) != expected:
            raise RuntimeError(f"prepared input changed: {filename}")
    table_dir = Path(manifest["table_dir"])
    if sha256(table_dir / "family.manifest") != manifest["family_manifest_sha256"]:
        raise RuntimeError("family manifest changed")
    driver_build = build_path(build_id)
    freeze = json.loads((EXPERIMENT / f"build_freeze{suffix(build_id)}.json").read_text())
    if sha256(driver_build / "star") != freeze["executable_sha256"] or archives() != freeze["private_archives"]:
        raise RuntimeError("driver or private archives changed")
    installed_before = installed_inventory()
    record(target / "installed_before.json", installed_before)
    cache_dir = target / "runtime_cache" / "eosDT_cache"
    cache_dir.mkdir(parents=True)
    for source in (MESA / "data/eosDT_data/cache").glob("*.bin"):
        shutil.copy2(source, cache_dir / source.name)
    temp_dir = target / "runtime_temp_cache"
    temp_dir.mkdir()
    runtime_env = environment()
    runtime_env.update(MESA_CACHES_DIR=str(cache_dir.parent), MESA_TEMP_CACHES_DIR=str(temp_dir))
    timeout = min(900, int(STOP_TIME - time.time() - 90))
    if timeout < 30:
        raise RuntimeError("insufficient time before numerical stop")
    command = [str(driver_build / "star")]
    record(target / "command.json", {"execution": execution, "command": command,
           "cwd": str(target), "backend": backend, "pair_id": pair_id,
           "build_id": build_id, "lead_clearance": clearance,
           "timeout_seconds": timeout, "threads": 1, "address_space_bytes": 16 * 1024**3})
    start = time.monotonic()
    timed_out = False
    with (target / "run.log").open("w") as stream:
        process = subprocess.Popen(command, cwd=target, env=runtime_env,
                                   stdout=stream, stderr=subprocess.STDOUT, preexec_fn=memory_limit)
        record(target / "process.json", {"pid": process.pid, "started": time.time()})
        try:
            returncode = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            process.kill()
            returncode = process.wait()
    record(target / "outcome.json", {"returncode": returncode, "timed_out": timed_out,
           "seconds": time.monotonic() - start, "completion_written": (target / "completion.txt").exists()})
    installed_after = installed_inventory()
    record(target / "installed_after.json", installed_after)
    if installed_before != installed_after:
        raise RuntimeError("installed libraries, modules or EOS caches changed during run")


def analyze(pair_id):
    results = {}
    for backend in ("baseline", "candidate"):
        target = run_path(backend, pair_id)
        outcome = json.loads((target / "outcome.json").read_text())
        completion = {}
        if (target / "completion.txt").exists():
            for line in (target / "completion.txt").read_text().splitlines():
                key, value = line.split(maxsplit=1)
                completion[key] = value.strip()
        events = (target / "events.tsv").read_text().splitlines() if (target / "events.tsv").exists() else []
        accepted = [line.split() for line in events if line.strip().startswith("accepted")]
        iterations = ((target / "iterations.tsv").read_text().splitlines()
                      if (target / "iterations.tsv").exists() else [])
        iteration_values = [int(line.split()[2]) for line in iterations if line and not line.startswith("#")]
        manifest = json.loads((target / "manifest.json").read_text())
        results[backend] = {"outcome": outcome, "completion": completion,
                            "max_model_number": manifest.get("max_model_number", 8),
                            "accepted_events": len(accepted),
                            "accepted_ages_yr": [float(row[4]) for row in accepted],
                            "accepted_dt_s": [float(row[5]) for row in accepted],
                            "final_center_h1": float(accepted[-1][10]) if accepted else None,
                            "final_center_he3": float(accepted[-1][12]) if accepted else None,
                            "max_retries_at_acceptance": max((int(row[3]) for row in accepted), default=None),
                            "max_newton_iteration": max(iteration_values, default=None)}
    baseline = results["baseline"]
    candidate = results["candidate"]
    requested = baseline["max_model_number"]
    if baseline["accepted_events"] == candidate["accepted_events"] and baseline["accepted_events"]:
        results["paired_differences"] = {
            "max_abs_age_yr": max(abs(left - right) for left, right in zip(
                baseline["accepted_ages_yr"], candidate["accepted_ages_yr"])),
            "max_abs_dt_s": max(abs(left - right) for left, right in zip(
                baseline["accepted_dt_s"], candidate["accepted_dt_s"])),
            "final_center_h1": candidate["final_center_h1"] - baseline["final_center_h1"],
            "final_center_he3": candidate["final_center_he3"] - baseline["final_center_he3"],
        }
    results["checks"] = {
        "both_exited_zero": all(entry["outcome"]["returncode"] == 0 for entry in (baseline, candidate)),
        "both_completed_requested": candidate["max_model_number"] == requested and all(
                                    entry["completion"].get("model_number") == str(requested)
                                    and entry["completion"].get("termination_code") == "4"
                                    and entry["accepted_events"] == requested for entry in (baseline, candidate)),
        "native_only_candidate": baseline["completion"].get("native_successes") == "0"
                                 and int(candidate["completion"].get("native_successes", "0")) > 0,
        "native_zero_failures": all(entry["completion"].get("native_failures") == "0"
                                    for entry in (baseline, candidate)),
        "native_loaded": int(candidate["completion"].get("native_loads", "0")) > 0,
        "burning_and_mixing": all(entry["completion"].get("saw_burn") == "T"
                                  and entry["completion"].get("saw_mix") == "T"
                                  and entry["completion"].get("saw_abundance_equations") == "T"
                                  for entry in (baseline, candidate)),
        "composition_changed": all(float(entry["completion"].get("delta_center_he3", "0")) != 0
                                   or float(entry["completion"].get("delta_center_h1", "0")) != 0
                                   for entry in (baseline, candidate)),
    }
    record(EXPERIMENT / f"analysis{suffix(pair_id)}.json", results)
    print(json.dumps(results["checks"], indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    building = commands.add_parser("build")
    building.add_argument("--build-id", default="")
    analysis = commands.add_parser("analyze")
    analysis.add_argument("--pair", default="")
    preparation = commands.add_parser("prepare")
    preparation.add_argument("table_dir")
    preparation.add_argument("--pair", default="")
    preparation.add_argument("--max-model-number", type=int, default=8)
    execution = commands.add_parser("run")
    execution.add_argument("backend", choices=("baseline", "candidate"))
    execution.add_argument("--lead-clearance", default="")
    execution.add_argument("--pair", default="")
    execution.add_argument("--build-id", default="")
    args = parser.parse_args()
    if args.action == "build":
        build(args.build_id)
    elif args.action == "prepare":
        prepare(args.table_dir, args.pair, args.max_model_number)
    elif args.action == "analyze":
        analyze(args.pair)
    else:
        run(args.backend, args.lead_clearance, args.pair, args.build_id)


if __name__ == "__main__":
    main()
