"""Explicitly cleared TASK052_same_seed_100; unchanged frozen stellar executable."""

import argparse
import json
import shutil
import subprocess
import time

from harness import (CODE, CONTROL, IDENTITY, INPUT_NAMES, MESA, RUNS, SEED_SHA, digest, environment,
                     inventory, limit_memory, placement, verify_private_closure, write_json)


EXTENSION = "TASK052_same_seed_100"
NAMES = {"baseline_100": "BASELINE", "candidate_100": "CANDIDATE"}
RAW_NAMES = ("queries.tsv", "component_queries.tsv", "events.tsv", "iterations.tsv", "matrices.bin", "LOGS", "photos", "run.log")


def raw_bytes(directory):
    total = 0
    for name in RAW_NAMES:
        path = directory / name
        if path.is_dir():
            total += sum(item.stat().st_size for item in path.rglob("*") if item.is_file())
        elif path.exists():
            total += path.stat().st_size
    return total


def prepare(name):
    execution = placement()
    verify_private_closure()
    target = RUNS / name
    target.mkdir(exist_ok=False)
    for filename in INPUT_NAMES:
        shutil.copyfile(CONTROL / filename, target / filename)
    shutil.copyfile(IDENTITY, target / "solar.identity")
    original = (target / "inlist").read_text()
    old = "   max_model_number = 8\n"
    new = "   max_model_number = 100\n"
    if original.count(old) != 1:
        raise ValueError("original eight-step limit is not unique")
    (target / "inlist").write_text(original.replace(old, new))
    write_json(target / "manifest.json", {"execution": execution, "source": str(CONTROL), "extension": EXTENSION,
               "scope": "explicitly cleared100-step exploratory frozen-composition relaxation",
               "approval": "lead/user explicit clearance at19:31PDT on2026-10-29",
               "inputs": {filename: digest(target / filename) for filename in INPUT_NAMES + ("solar.identity",)},
               "only_input_change_from_first_pair": {"max_model_number": [8, 100]},
               "original_inlist_sha256": digest(CONTROL / "inlist"),
               "force_retry": True, "steps": 100, "source_seed_sha256": SEED_SHA})
    print(target, flush=True)


def run(name):
    execution = placement()
    verify_private_closure()
    target = RUNS / name
    manifest = json.loads((target / "manifest.json").read_text())
    if manifest["extension"] != EXTENSION or manifest["steps"] != 100:
        raise ValueError("wrong extension identity")
    for filename, expected in manifest["inputs"].items():
        if digest(target / filename) != expected:
            raise ValueError("input changed: " + filename)
    original = (CONTROL / "inlist").read_text()
    if (target / "inlist").read_text() != original.replace("   max_model_number = 8\n", "   max_model_number = 100\n"):
        raise ValueError("unapproved control change")
    if (target / "command.json").exists():
        raise RuntimeError("preserve existing attempt")
    freeze = json.loads((RUNS / "driver_build/freeze.json").read_text())
    for filename, expected in freeze["sources"].items():
        if digest(CODE / filename) != expected:
            raise ValueError("frozen source changed: " + filename)
    before = inventory()
    executable = RUNS / "driver_build/star"
    if before != freeze["dependencies"] or digest(executable) != freeze["executable_sha256"]:
        raise RuntimeError("frozen dependency/executable changed")
    if shutil.disk_usage(RUNS).free < 6 * 1024**3:
        raise RuntimeError("less than6GiB free for this private-cache run")
    cache_root = target / "runtime_cache"
    eos_cache = cache_root / "eosDT_cache"
    eos_cache.mkdir(parents=True)
    for original_cache in (MESA / "data/eosDT_data/cache").glob("*.bin"):
        shutil.copy2(original_cache, eos_cache / original_cache.name)
    temporary = target / "runtime_temp_cache"
    temporary.mkdir()
    runtime_environment = environment()
    runtime_environment.update(TASK052_BACKEND=NAMES[name], MESA_CACHES_DIR=str(cache_root), MESA_TEMP_CACHES_DIR=str(temporary))
    timeout = min(900, int(1793334000 - time.time() - 120))
    if timeout < 30:
        raise RuntimeError("insufficient time before hard stop")
    command = [str(executable)]
    write_json(target / "command.json", {"execution": execution, "command": command, "cwd": str(target),
               "extension": EXTENSION, "backend": NAMES[name], "accepted_step_limit": 100,
               "executable_sha256": digest(executable), "timeout_seconds": timeout,
               "numerical_threads": 1, "memory_limit_bytes": 8 * 1024**3, "force_retry": True,
               "raw_logging_cap_bytes": 1024**3})
    write_json(target / "before.json", before)
    start = time.monotonic()
    stop_reason = None
    with (target / "run.log").open("w") as stream:
        process = subprocess.Popen(command, cwd=target, env=runtime_environment, stdout=stream,
                                   stderr=subprocess.STDOUT, preexec_fn=limit_memory)
        write_json(target / "process.json", {"pid": process.pid, "started": time.time(), "backend": NAMES[name]})
        while process.poll() is None:
            elapsed = time.monotonic() - start
            if elapsed > timeout:
                stop_reason = "timeout"
            elif raw_bytes(target) > 1024**3:
                stop_reason = "raw_logging_cap"
            if stop_reason:
                process.kill()
                break
            time.sleep(.5)
        returncode = process.wait()
    seconds = time.monotonic() - start
    after = inventory()
    write_json(target / "after.json", after)
    outcome = {"returncode": returncode, "timed_out": stop_reason == "timeout", "stop_reason": stop_reason,
               "seconds": seconds, "protected_unchanged": before == after,
               "completion_written": (target / "completion.txt").exists(), "raw_logging_bytes": raw_bytes(target)}
    write_json(target / "outcome.json", outcome)
    print(json.dumps(outcome), flush=True)
    if returncode or stop_reason or before != after or not outcome["completion_written"]:
        raise RuntimeError("bounded extension stopped or failed; preserve attempt and report")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run"))
    parser.add_argument("name", choices=tuple(NAMES))
    arguments = parser.parse_args()
    (prepare if arguments.action == "prepare" else run)(arguments.name)


if __name__ == "__main__":
    main()
