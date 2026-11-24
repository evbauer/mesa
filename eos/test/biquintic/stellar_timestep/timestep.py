"""Cleared four-run TASK057 refinement with the unchanged052stellar executable."""

import argparse
from datetime import datetime
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

sys.dont_write_bytecode = True
CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE.parent / "stellar_candidate"))
import harness as frozen

EXPERIMENT = frozen.WORKSPACE / "freeeos_biquintic_progress/experiments/057_stellar_timestep"
RUNS = EXPERIMENT / "runs"
ANALYSIS = EXPERIMENT / "analysis"
TARGET_AGE = 1089056.4362468508
NAMES = ("coarse_baseline", "coarse_candidate", "half_baseline", "half_candidate")
SOURCE = frozen.RUNS / "baseline_100"
EXECUTABLE = frozen.RUNS / "driver_build/star"
LAUNCH_STOP = datetime.fromisoformat("2026-10-29T20:25:00-07:00").timestamp()
NUMERICAL_STOP = datetime.fromisoformat("2026-10-29T20:30:00-07:00").timestamp()
WHITELIST = {("star_job", "years_for_initial_dt"), ("controls", "max_years_for_timestep"),
             ("controls", "max_model_number"), ("controls", "max_age")}
RAW_NAMES = ("queries.tsv", "component_queries.tsv", "events.tsv", "iterations.tsv", "matrices.bin", "LOGS", "photos", "run.log")


def assignments(text):
    result = {}
    block = None
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("&"):
            block = stripped[1:]
        elif stripped == "/":
            block = None
        elif "=" in stripped:
            match = re.fullmatch(r"([a-zA-Z0-9_()]+)\s*=\s*(.+)", stripped)
            if not match or block is None:
                raise ValueError("unexpected inlist syntax")
            key = (block, match.group(1))
            if key in result:
                raise ValueError("duplicate inlist key")
            result[key] = match.group(2)
    return result


def expected_inlist(name):
    original = (SOURCE / "inlist").read_text()
    changes = {"   max_model_number = 100\n": "   max_model_number = 300\n   max_age = 1089056.4362468508d0\n"}
    if name.startswith("half_"):
        changes.update({"   years_for_initial_dt = 100\n": "   years_for_initial_dt = 50\n",
                        "   max_years_for_timestep = 1000\n": "   max_years_for_timestep = 500\n"})
    result = original
    for old, new in changes.items():
        if result.count(old) != 1:
            raise ValueError("nonunique whitelisted input line")
        result = result.replace(old, new)
    before, after = assignments(original), assignments(result)
    differences = {key for key in before.keys() | after.keys() if before.get(key) != after.get(key)}
    if not differences <= WHITELIST:
        raise ValueError("nonwhitelisted input change")
    return result, {"/".join(key): [before.get(key), after.get(key)] for key in sorted(differences)}


def verify_frozen(full=False):
    frozen.placement()
    frozen.verify_private_closure()
    freeze = json.loads((frozen.RUNS / "driver_build/freeze.json").read_text())
    if frozen.digest(EXECUTABLE) != freeze["executable_sha256"] or frozen.inventory() != freeze["dependencies"]:
        raise ValueError("052executable/private or installed dependency changed")
    for filename, expected in freeze["sources"].items():
        if frozen.digest(frozen.CODE / filename) != expected:
            raise ValueError("052frozen driver source changed: " + filename)
    if full:
        handoff = json.loads((frozen.SUMMARY / "FINAL_HANDOFF.json").read_text())
        for relative, expected in handoff["artifact_hashes"].items():
            if frozen.digest(frozen.WORKSPACE / relative) != expected:
                raise ValueError("052frozen artifact changed: " + relative)
        return len(handoff["artifact_hashes"])
    return freeze


def raw_bytes(directory):
    total = 0
    for name in RAW_NAMES:
        path = directory / name
        if path.is_dir():
            total += sum(item.stat().st_size for item in path.rglob("*") if item.is_file())
        elif path.exists():
            total += path.stat().st_size
    return total


def prepare():
    execution = frozen.placement()
    if RUNS.exists():
        raise RuntimeError("preserve existing057preparation")
    verified = verify_frozen(full=True)
    if shutil.disk_usage(EXPERIMENT).free < 28 * 1024**3:
        raise RuntimeError("less than authorized28GiB disk reserve")
    RUNS.mkdir()
    ANALYSIS.mkdir(exist_ok=True)
    source_inputs = json.loads((SOURCE / "manifest.json").read_text())["inputs"]
    for name in NAMES:
        target = RUNS / name
        target.mkdir()
        for filename, expected in source_inputs.items():
            if frozen.digest(SOURCE / filename) != expected:
                raise ValueError("052source input changed")
            shutil.copyfile(SOURCE / filename, target / filename)
        text, changes = expected_inlist(name)
        (target / "inlist").write_text(text)
        frozen.write_json(target / "manifest.json", {"execution": execution, "task": "TASK057", "name": name,
                          "approval": "explicit lead/user19:48PDT2026-10-29; exact frozen052proposal",
                          "source": str(SOURCE), "allowed_input_differences": changes,
                          "inputs": {filename: frozen.digest(target / filename) for filename in source_inputs},
                          "target_age_years": TARGET_AGE, "endpoint_ulp_tolerance": 2, "safety_model_cap": 300,
                          "requested_termination": "max_age", "force_retry": True})
    frozen.write_json(EXPERIMENT / "preparation.json", {"execution": execution, "verified052artifacts": verified,
                      "proposal_sha256": frozen.digest(frozen.SUMMARY / "TASK057_PROPOSAL.md"),
                      "freeze_sha256": frozen.digest(frozen.SUMMARY / "FINAL_HANDOFF.json"),
                      "executable_sha256": frozen.digest(EXECUTABLE), "target_age_years": TARGET_AGE})
    print("prepared four independent runs; verified", verified, "frozen052artifacts", flush=True)


def check_order(name):
    if any(EXPERIMENT.glob("CONTINUATION_STOP_*.json")) and name != "half_candidate":
        raise RuntimeError("decoder continuation permits only the remaining half candidate")
    continuation = name in ("half_baseline", "half_candidate") and continuation_authorized()
    if (EXPERIMENT / "STOP.json").exists() and not continuation:
        raise RuntimeError("sequence stopped; no automatic rerun")
    for prior in NAMES[:NAMES.index(name)]:
        record = ANALYSIS / (prior + ".json")
        if not record.exists() or (not json.loads(record.read_text())["audit_passed"] and not (continuation and prior == "coarse_candidate")):
            raise RuntimeError("prior run not audited: " + prior)
    if name.startswith("half_"):
        for backend in ("baseline", "candidate"):
            if not json.loads((ANALYSIS / ("coarse_" + backend + ".json")).read_text())["coarse_reproduction_gate"]["passed"] and not (continuation and backend == "candidate"):
                raise RuntimeError("coarse reproduction gate not passed")


def continuation_authorized():
    path = EXPERIMENT / "CONTINUATION.json"
    if not path.exists():
        return False
    stops = list(EXPERIMENT.glob("CONTINUATION_STOP_*.json"))
    if stops:
        decoder_path = EXPERIMENT / "DECODER_CONTINUATION.json"
        expected_stop = EXPERIMENT / "CONTINUATION_STOP_half_baseline.json"
        if stops != [expected_stop] or not decoder_path.exists():
            raise RuntimeError("exploratory continuation already stopped")
        decoder = json.loads(decoder_path.read_text())
        if decoder.get("authorization_id") != "LEAD_2026_10_29_2014_INDEX_DECODER_AND_LAST_HALF_CANDIDATE":
            raise ValueError("unrecognized decoder continuation")
        if decoder.get("authorized_run") != "half_candidate":
            raise ValueError("decoder continuation permits no other run")
        if frozen.digest(expected_stop) != "36ec3412d58530e412ac3f8f48fd5e5aa8453919e347bd68bfdaef3257791f08":
            raise ValueError("original decoder stop changed")
        decoded = json.loads((ANALYSIS / "half_baseline.json").read_text())
        if not decoded["audit_passed"] or not decoded["endpoint"]["passed"]:
            raise ValueError("corrected half-baseline audit did not pass")
    authorization = json.loads(path.read_text())
    if authorization.get("authorization_id") != "LEAD_2026_10_29_HALF_ONLY_AFTER_COARSE_GATE_STOP":
        raise ValueError("unrecognized continuation authorization")
    if authorization.get("authorized_runs") != ["half_baseline", "half_candidate"] or authorization.get("waived_blocker") != "final_selected_profile_ceiling_only":
        raise ValueError("continuation scope is not the specific reviewed half pair")
    expected = {
        EXPERIMENT / "STOP.json": "7163096f1ded72d1cba58b080e8fb6b86cd509064625d0e405bdb6381e3778d5",
        ANALYSIS / "coarse_candidate.json": "9a03970330e1034b9fb1d9f53361c773b1e97fe7b54c4323e80662f5d53841fd",
        ANALYSIS / "coarse_baseline.json": "5e2df7224d59d4821f0f2b844cfc9216d3874373bee8ebe5c5c03c70ff3093a8"}
    for record, checksum in expected.items():
        if frozen.digest(record) != checksum:
            raise ValueError("reviewed original stop/audit changed: " + str(record))
    candidate = json.loads((ANALYSIS / "coarse_candidate.json").read_text())
    gate = candidate["coarse_reproduction_gate"]
    required = ("accepted_count_matches052", "pre_endpoint_events_exact", "pre_endpoint_profile_bytes_exact", "pre_endpoint_matrix_bytes_exact")
    if not all(gate[key] for key in required) or gate["final_scaled_roundoff_ceiling"] != 5e-11:
        raise ValueError("continuation would waive an unreviewed coarse gate")
    if not candidate["endpoint"]["passed"] or candidate["outcome"]["returncode"] != 0 or not candidate["outcome"]["protected_unchanged"]:
        raise ValueError("coarse endpoint/outcome is not supported for continuation")
    if candidate["query_audit"]["component_errors"] != 0 or candidate["query_audit"]["hybrid_failures"] != 0:
        raise ValueError("coarse support/hybrid error cannot be waived")
    if gate["final_dt_difference_seconds"] > gate["final_dt_roundoff_ceiling_seconds"]:
        raise ValueError("coarse endpoint timestep gate cannot be waived")
    from timestep_audit import coarse_gate, mesa_profile, records
    accepted = [row for row in records(RUNS / "coarse_candidate/events.tsv") if row[0] == "accepted"]
    current_gate = coarse_gate("coarse_candidate", accepted, mesa_profile(RUNS / "coarse_candidate/LOGS/profile100.data"))
    if current_gate != gate:
        raise ValueError("reviewed candidate coarse prefix or endpoint evidence changed")
    return True


def run(name):
    check_order(name)
    execution = frozen.placement()
    if time.time() >= LAUNCH_STOP:
        raise RuntimeError("no new run after20:25")
    freeze = verify_frozen()
    target = RUNS / name
    manifest = json.loads((target / "manifest.json").read_text())
    for filename, expected in manifest["inputs"].items():
        if frozen.digest(target / filename) != expected:
            raise ValueError("prepared input changed: " + filename)
    if (target / "inlist").read_text() != expected_inlist(name)[0]:
        raise ValueError("nonwhitelisted control difference")
    if (target / "command.json").exists():
        raise RuntimeError("preserve existing attempt")
    before = frozen.inventory()
    cache_root = target / "runtime_cache"
    eos_cache = cache_root / "eosDT_cache"
    eos_cache.mkdir(parents=True)
    for original_cache in (frozen.MESA / "data/eosDT_data/cache").glob("*.bin"):
        shutil.copy2(original_cache, eos_cache / original_cache.name)
    temporary = target / "runtime_temp_cache"
    temporary.mkdir()
    runtime_environment = frozen.environment()
    backend = "CANDIDATE" if name.endswith("candidate") else "BASELINE"
    runtime_environment.update(TASK052_BACKEND=backend, MESA_CACHES_DIR=str(cache_root), MESA_TEMP_CACHES_DIR=str(temporary))
    if time.time() >= LAUNCH_STOP:
        raise RuntimeError("cache preparation passed20:25launch cutoff")
    timeout = min(300, int(NUMERICAL_STOP - time.time()))
    command = [str(EXECUTABLE)]
    frozen.write_json(target / "command.json", {"execution": execution, "command": command, "cwd": str(target),
                      "backend": backend, "executable_sha256": freeze["executable_sha256"],
                      "timeout_seconds": timeout, "numerical_threads": 1, "memory_limit_bytes": 8 * 1024**3,
                      "raw_logging_cap_bytes": 2 * 1024**3, "force_retry": True, "requested_termination": "max_age",
                      "exploratory_continuation": name.startswith("half_"),
                      "continuation_authorization_sha256": frozen.digest(EXPERIMENT / "CONTINUATION.json") if name.startswith("half_") else None})
    frozen.write_json(target / "before.json", before)
    start = time.monotonic()
    stop_reason = None
    with (target / "run.log").open("w") as stream:
        process = subprocess.Popen(command, cwd=target, env=runtime_environment, stdout=stream,
                                   stderr=subprocess.STDOUT, preexec_fn=frozen.limit_memory)
        frozen.write_json(target / "process.json", {"pid": process.pid, "started": time.time(), "backend": backend})
        while process.poll() is None:
            if time.monotonic() - start >= timeout or time.time() >= NUMERICAL_STOP:
                stop_reason = "timeout_or20:30cutoff"
            elif raw_bytes(target) > 2 * 1024**3:
                stop_reason = "raw_logging_cap"
            if stop_reason:
                process.kill()
                break
            time.sleep(.5)
        returncode = process.wait()
    seconds = time.monotonic() - start
    after = frozen.inventory()
    frozen.write_json(target / "after.json", after)
    outcome = {"returncode": returncode, "stop_reason": stop_reason, "seconds": seconds,
               "protected_unchanged": before == after, "completion_written": (target / "completion.txt").exists(),
               "raw_logging_bytes": raw_bytes(target)}
    frozen.write_json(target / "outcome.json", outcome)
    print(name, json.dumps(outcome), flush=True)
    if returncode or stop_reason or before != after or not outcome["completion_written"]:
        raise RuntimeError("run stopped or failed; preserve attempt")
    from timestep_audit import audit_run
    audit_run(name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run"))
    parser.add_argument("name", nargs="?", choices=NAMES)
    arguments = parser.parse_args()
    try:
        if arguments.action == "prepare":
            prepare()
        elif arguments.name:
            run(arguments.name)
        else:
            raise ValueError("run requires name")
    except Exception as error:
        stop_path = EXPERIMENT / ("CONTINUATION_STOP_" + (arguments.name or arguments.action) + ".json") if (EXPERIMENT / "STOP.json").exists() else EXPERIMENT / "STOP.json"
        if stop_path.exists():
            stop_path = stop_path.with_name(stop_path.stem + "_" + str(int(time.time())) + ".json")
        frozen.write_json(stop_path, {"time": time.time(), "action": arguments.action,
                          "run": arguments.name, "error_type": type(error).__name__, "error": str(error),
                          "policy": "preserve attempt; no automatic rerun, retune or further launch"})
        raise


if __name__ == "__main__":
    main()
