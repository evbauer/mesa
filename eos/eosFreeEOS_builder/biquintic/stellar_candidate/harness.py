"""Authorized, bounded TASK052 stellar A/B with one coherent private executable."""

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
EXPERIMENT = WORKSPACE / "freeeos_biquintic_progress/experiments/052_stellar_api"
RUNS = EXPERIMENT / "runs"
SUMMARY = EXPERIMENT / "candidate_summary"
PRIVATE = EXPERIMENT / "build"
PRIVATE_EOS = WORKSPACE / "freeeos_biquintic_progress/experiments/050_runtime_api/build/eos"
CONTROL = WORKSPACE / "freeeos_biquintic_progress/experiments/049_stellar_preparation/control_legacy"
IDENTITY = WORKSPACE / "freeeos_biquintic_progress/experiments/047_runtime_adapter/solar.identity"
CACHE = WORKSPACE / "freeeos_biquintic_progress/experiments/008_standalone_table/solar_split.cache"
PATTERN = MESA / "eos/eosFreeEOS_builder/mass_frac_solar.txt"
SEED_SHA = "cd5b74e7ccd55d8253d045abc5bb95cba1506ba42fba70482bccf5f29cf9eb3e"
IDENTITY_SHA = "ac3eef30bba5389bb8c04f9adf56a22ce053e4df47260f27e0bb00fd63f03ea3"
CACHE_SHA = "e2c72c6f02cb900649a6a629bfd4fc53c2e8423e5bbe6afd6a189c45cac3a204"
PATTERN_SHA = "bd62922804e073a2a410480994b8a28e85b0bc53e469729d6ed5136ea9da10a6"
MUTABLE_PACKAGES = ("star", "star_data", "net", "kap", "atm", "rates")
INPUT_NAMES = ("seed.mod", "inlist", "frozen20.net", "uniform_xa.data", "history_columns.list", "profile_columns.list")
CORE_SOURCES = ("eos/private/helm_polynomials.f90", "eos/eosFreeEOS_builder/biquintic/test_fortran/freeeos_biquintic.f90",
                "eos/eosFreeEOS_builder/biquintic/table_fortran/freeeos_table.f90",
                "eos/eosFreeEOS_builder/biquintic/runtime_adapter/runtime_sha256.f90",
                "eos/eosFreeEOS_builder/biquintic/runtime_adapter/freeeos_runtime_adapter.f90")


def digest(path):
    checksum = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            checksum.update(chunk)
    return checksum.hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def placement():
    host = socket.gethostname().split(".")[0]
    if host != "rzwhippet14" or os.environ.get("SLURM_JOB_ID") != "5528075":
        raise RuntimeError("TASK052 requires inherited job5528075 on rzwhippet14")
    if time.time() > 1793334000:
        raise RuntimeError("TASK052 hard stop at21:20PDT")
    return {"hostname": host, "job": "5528075", "time": time.time()}


def environment():
    result = dict(os.environ)
    result.update(MESA_DIR=str(MESA), MESASDK_ROOT=str(SDK), OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                  MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1", BLIS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
                  TASK052_FORCE_RETRY="1", TASK052_CACHE=str(CACHE), TASK052_PATTERN=str(PATTERN),
                  LD_LIBRARY_PATH=f"{SDK}/lib:{SDK}/lib64:" + result.get("LD_LIBRARY_PATH", ""))
    result["PATH"] = f"{SDK}/bin:" + result["PATH"]
    return result


def limit_memory():
    resource.setrlimit(resource.RLIMIT_AS, (8 * 1024**3, 8 * 1024**3))


def inventory():
    paths = list((MESA / "build").glob("*/lib/*.a")) + list((MESA / "build").glob("*/include/*.mod"))
    paths += list(PRIVATE.glob("*/lib/*.a")) + list(PRIVATE.glob("*/include/*.mod"))
    paths += [MESA / path for path in CORE_SOURCES]
    paths += [MESA / "eos/plotter/README.rst", MESA / "eos/plotter/src/eos_plotter.f90", PATTERN, IDENTITY]
    return {str(path.relative_to(WORKSPACE)): digest(path) for path in sorted(set(paths))}


def verify_private_closure():
    for package in MUTABLE_PACKAGES:
        directory = PRIVATE / package
        if directory.is_symlink() or not (directory / "lib" / f"lib{package}.a").is_file():
            raise RuntimeError("missing coherent private package: " + package)
    if (PRIVATE / "eos").resolve() != PRIVATE_EOS.resolve():
        raise RuntimeError("private EOS link does not resolve to frozen050build")
    for directory in PRIVATE.iterdir():
        if directory.name in MUTABLE_PACKAGES or directory.name == "eos":
            continue
        if not directory.is_symlink() or directory.resolve() != (MESA / "build" / directory.name).resolve():
            raise RuntimeError("unexpected dependency link: " + str(directory))
    for path, expected in ((CONTROL / "seed.mod", SEED_SHA), (IDENTITY, IDENTITY_SHA), (CACHE, CACHE_SHA), (PATTERN, PATTERN_SHA)):
        if digest(path) != expected:
            raise RuntimeError("frozen input identity changed: " + str(path))
    expected_adapter = {
        "freeeos_runtime_adapter.f90": "868985fec27b85c0e413944d005529c0944e548ef424823c0507f36465d33233",
        "runtime_sha256.f90": "d98e666d918317a1bbf11971d6f3013b10973ed27e24e1b5a3dc1c762a3b8236"}
    for name, expected in expected_adapter.items():
        if digest(CODE.parent / "runtime_adapter" / name) != expected:
            raise RuntimeError("frozen047source changed: " + name)


def build():
    execution = placement()
    verify_private_closure()
    RUNS.mkdir(exist_ok=True)
    output = RUNS / f"build_record_{int(time.time())}"
    output.mkdir()
    before = inventory()
    write_json(output / "before.json", before)
    command = ["make", "-C", str(CODE), "-j8"]
    write_json(output / "command.json", {"execution": execution, "command": command,
               "approval": "explicit bounded exploratory stellar A/B; private driver build only",
               "sources": {path.name: digest(path) for path in CODE.iterdir() if path.is_file()}})
    with (output / "build.log").open("w") as stream:
        process = subprocess.run(command, env=environment(), stdout=stream, stderr=subprocess.STDOUT,
                                 timeout=600, preexec_fn=limit_memory)
    after = inventory()
    write_json(output / "after.json", after)
    if before != after:
        raise RuntimeError("protected inputs or private/installed dependencies changed during build")
    if process.returncode:
        raise RuntimeError(f"build failed: {output / 'build.log'}")
    mapping = (RUNS / "driver_build/link.map").read_text()
    archives = [Path(line[5:].strip()).resolve() for line in mapping.splitlines() if line.startswith("LOAD ") and line.endswith(".a")]
    linked = {}
    for package in MUTABLE_PACKAGES + ("eos",):
        expected = (PRIVATE / package / "lib" / f"lib{package}.a").resolve()
        actual = [path for path in archives if path.name == f"lib{package}.a"]
        if actual != [expected]:
            raise RuntimeError(f"unexpected linked archive for {package}: {actual}")
        linked[package] = {"archive": str(expected), "sha256": digest(expected)}
    write_json(output / "linked_private.json", linked)
    write_json(RUNS / "driver_build/freeze.json", {"executable_sha256": digest(RUNS / "driver_build/star"),
               "dependencies": after, "linked_private": linked,
               "sources": {path.name: digest(path) for path in CODE.iterdir() if path.is_file()}})
    print(f"coherent private executable built; records: {output}")


def prepare(name):
    execution = placement()
    verify_private_closure()
    if not re.fullmatch(r"[a-z][a-z0-9_]*", name):
        raise ValueError("use a simple run name")
    target = RUNS / name
    target.mkdir(exist_ok=False)
    for filename in INPUT_NAMES:
        shutil.copyfile(CONTROL / filename, target / filename)
    shutil.copyfile(IDENTITY, target / "solar.identity")
    write_json(target / "manifest.json", {"execution": execution, "source": str(CONTROL),
               "scope": "bounded exploratory frozen-composition stellar relaxation; no burning or global jump-stress claim",
               "inputs": {filename: digest(target / filename) for filename in INPUT_NAMES + ("solar.identity",)},
               "force_retry": True, "steps": 8, "source_seed_sha256": SEED_SHA})
    print(target)


def run(name, backend):
    execution = placement()
    verify_private_closure()
    target = RUNS / name
    manifest = json.loads((target / "manifest.json").read_text())
    for filename, expected in manifest["inputs"].items():
        if digest(target / filename) != expected:
            raise RuntimeError("prepared input changed: " + filename)
    if (target / "command.json").exists():
        raise RuntimeError("preserve old attempts; prepare a fresh directory")
    freeze = json.loads((RUNS / "driver_build/freeze.json").read_text())
    executable = RUNS / "driver_build/star"
    before = inventory()
    if digest(executable) != freeze["executable_sha256"] or before != freeze["dependencies"]:
        raise RuntimeError("frozen executable/dependencies changed")
    cache_root = target / "runtime_cache"
    eos_cache = cache_root / "eosDT_cache"
    eos_cache.mkdir(parents=True)
    for original in (MESA / "data/eosDT_data/cache").glob("*.bin"):
        shutil.copy2(original, eos_cache / original.name)
    temporary = target / "runtime_temp_cache"
    temporary.mkdir()
    runtime_environment = environment()
    runtime_environment.update(TASK052_BACKEND=backend, MESA_CACHES_DIR=str(cache_root), MESA_TEMP_CACHES_DIR=str(temporary))
    timeout = min(900, int(1793334000 - time.time() - 120))
    if timeout < 30:
        raise RuntimeError("insufficient time before hard stop")
    command = [str(executable)]
    write_json(target / "command.json", {"execution": execution, "command": command, "cwd": str(target),
               "backend": backend, "executable_sha256": digest(executable), "timeout_seconds": timeout,
               "numerical_threads": 1, "memory_limit_bytes": 8 * 1024**3, "force_retry": True})
    write_json(target / "before.json", before)
    start = time.monotonic()
    timed_out = False
    with (target / "run.log").open("w") as stream:
        process = subprocess.Popen(command, cwd=target, env=runtime_environment, stdout=stream,
                                   stderr=subprocess.STDOUT, preexec_fn=limit_memory)
        write_json(target / "process.json", {"pid": process.pid, "started": time.time(), "backend": backend})
        try:
            returncode = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            process.kill()
            returncode = process.wait()
    after = inventory()
    write_json(target / "after.json", after)
    result = {"returncode": returncode, "timed_out": timed_out, "seconds": time.monotonic() - start,
              "protected_unchanged": before == after, "completion_written": (target / "completion.txt").exists()}
    write_json(target / "outcome.json", result)
    print(json.dumps(result))
    if before != after:
        raise RuntimeError("protected dependencies changed during run")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    commands.add_parser("build")
    preparation = commands.add_parser("prepare")
    preparation.add_argument("name")
    execution = commands.add_parser("run")
    execution.add_argument("name")
    execution.add_argument("backend", choices=("BASELINE", "CANDIDATE"))
    arguments = parser.parse_args()
    if arguments.action == "build":
        build()
    elif arguments.action == "prepare":
        prepare(arguments.name)
    else:
        run(arguments.name, arguments.backend)


if __name__ == "__main__":
    main()
