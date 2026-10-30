"""Prepare and run bounded installed-LEGACY frozen-composition controls."""

import argparse
import hashlib
import json
import math
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
EXPERIMENT = WORKSPACE / "freeeos_biquintic_progress/experiments/049_stellar_preparation"
SEED = MESA / "star/test_suite/simplex_solar_calibration/standard_late_pre_zams_1.0M.mod"
PATTERN = MESA / "eos/eosFreeEOS_builder/mass_frac_solar.txt"
PATTERN_SHA = "bd62922804e073a2a410480994b8a28e85b0bc53e469729d6ed5136ea9da10a6"
NAMES = "h1 he4 c12 n14 o16 ne20 na23 mg24 al27 si28 p31 s32 cl35 ar40 ca40 ti48 cr52 mn55 fe56 ni58".split()
ELEMENTS = "C N O Ne Na Mg Al Si P S Cl Ar Ca Ti Cr Mn Fe Ni".split()


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2) + "\n")


def placement():
    host = socket.gethostname().split(".")[0]
    if host != "rzwhippet14" or os.environ.get("SLURM_JOB_ID") != "5528075":
        raise RuntimeError("TASK049 requires inherited job 5528075 on rzwhippet14")
    if time.time() > 1793334000:
        raise RuntimeError("TASK049 hard stop: 2026-10-29 21:20 PDT")
    return {"hostname": host, "job": os.environ["SLURM_JOB_ID"], "time": time.time()}


def environment():
    result = dict(os.environ)
    result.update(MESA_DIR=str(MESA), MESASDK_ROOT=str(SDK), OMP_NUM_THREADS="1",
                  OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1",
                  VECLIB_MAXIMUM_THREADS="1", BLIS_NUM_THREADS="1",
                  LD_LIBRARY_PATH=f"{SDK}/lib:{SDK}/lib64:" + result.get("LD_LIBRARY_PATH", ""))
    result["PATH"] = f"{SDK}/bin:" + result["PATH"]
    return result


def limit_memory():
    resource.setrlimit(resource.RLIMIT_AS, (8 * 1024**3, 8 * 1024**3))


def inventory():
    paths = sorted((MESA / "build").glob("*/lib/*.a"))
    paths += sorted((MESA / "build").glob("*/include/*.mod"))
    paths += [MESA / "eos/plotter/README.rst", MESA / "eos/plotter/src/eos_plotter.f90", PATTERN, SEED]
    return {str(path.relative_to(WORKSPACE)): digest(path) for path in paths}


def pinned_mixture():
    if digest(PATTERN) != PATTERN_SHA:
        raise RuntimeError("audited metal-pattern identity changed")
    rows = [line.split() for line in PATTERN.read_text().splitlines()[1:] if line.strip()]
    if [row[0] for row in rows] != ELEMENTS:
        raise RuntimeError("metal-pattern element order changed")
    weights = [float(row[1]) for row in rows]
    total = sum(weights)
    if abs(total - 1) > 1e-12:
        raise RuntimeError("metal-pattern normalization changed")
    return [0.70, 0.28] + [0.02 * value / total for value in weights]


def prepare(name, seed=SEED):
    execution = placement()
    if not re.fullmatch(r"[a-z][a-z0-9_]*", name):
        raise ValueError("run name must be a simple lowercase identifier")
    target = EXPERIMENT / name
    target.mkdir(exist_ok=False)
    masses = pinned_mixture()
    seed = Path(seed).resolve()
    original = seed.read_text().splitlines()
    header_index = next(index for index, line in enumerate(original) if line.split()[:2] == ["lnd", "lnT"])
    shells = int(next(line.split()[-1] for line in original[:header_index] if line.split()[:1] == ["n_shells"]))
    header = original[:header_index]
    header = [re.sub(r"('basic.net')", "'frozen20.net'", line) for line in header]
    header = [re.sub(r"(\s+species\s+)\d+", r"\g<1>20", line) for line in header]
    columns = original[header_index].split()
    first_abundance = columns.index("h1")
    data = original[header_index + 1:header_index + 1 + shells]
    transformed = []
    states = []
    for expected, line in enumerate(data, 1):
        fields = line.split()
        if int(fields[0]) != expected or len(fields) != len(columns) + 1:
            raise ValueError("unexpected saved-model row layout")
        structural = fields[:1 + first_abundance]
        transformed.append(" ".join(structural + [f"{value:.17E}".replace("E", "D") for value in masses]))
        log_rho = float(fields[1].replace("D", "E")) / math.log(10)
        log_temp = float(fields[2].replace("D", "E")) / math.log(10)
        states.append({"zone": expected, "logRho": log_rho, "logT": log_temp, "q": log_rho - 2 * log_temp + 12})
    output = header + [" ".join(columns[:first_abundance] + NAMES)] + transformed + original[header_index + 1 + shells:]
    (target / "seed.mod").write_text("\n".join(output) + "\n")
    (target / "frozen20.net").write_text("add_isos(\n" + ",\n".join(NAMES) + "\n)\n")
    (target / "uniform_xa.data").write_text("".join(f"{name} {mass:.17e}\n" for name, mass in zip(NAMES, masses)))
    shutil.copyfile(CODE / "inlist.template", target / "inlist")
    (target / "history_columns.list").write_text("model_number\nstar_age\nlog_dt\nstar_mass\nlog_L\nlog_R\nlog_Teff\nnum_retries\nnum_iters\n")
    (target / "profile_columns.list").write_text("zone\nmass\nlogRho\nlogT\nlogP\nlogR\nluminosity\ncp\ncv\ngrada\neps_grav\neps_nuc\neos_frac_FreeEOS\neos_frac_HELM\neos_frac_OPAL_SCVH\neos_frac_PC\neos_frac_Skye\neos_frac_CMS\neos_frac_ideal\n" + "\n".join(NAMES) + "\n")
    write_json(target / "seed_states.json", states)
    write_json(target / "manifest.json", {"execution": execution, "backend": "LEGACY", "candidate_gate": "closed",
               "source_seed": str(seed.relative_to(WORKSPACE)), "source_seed_sha256": digest(seed),
               "derived_seed_sha256": digest(target / "seed.mod"), "pattern_sha256": PATTERN_SHA,
               "species": NAMES, "mass_fractions": masses, "shells": shells,
               "derivation": "structural columns unchanged; network and every cell composition repinned; requires legacy warmup",
               "logT_range": [min(state["logT"] for state in states), max(state["logT"] for state in states)],
               "logRho_range": [min(state["logRho"] for state in states), max(state["logRho"] for state in states)],
               "source_code_sha256": {path.name: digest(path) for path in CODE.iterdir() if path.is_file()}})
    print(target)


def build():
    execution = placement()
    before = inventory()
    write_json(EXPERIMENT / "installed_before_build.json", before)
    command = ["make", "-C", str(CODE), "-j8"]
    write_json(EXPERIMENT / "build_command.json", {"command": command, "execution": execution,
               "abi_policy": "all module includes and linked archives from installed mesa/build; no private EOS overlay"})
    log = EXPERIMENT / f"build_{int(time.time())}.log"
    with log.open("w") as stream:
        process = subprocess.run(command, env=environment(), stdout=stream, stderr=subprocess.STDOUT,
                                 timeout=600, preexec_fn=limit_memory)
    after = inventory()
    write_json(EXPERIMENT / "installed_after_build.json", after)
    if before != after:
        raise RuntimeError("installed dependencies or protected inputs changed during build")
    print(f"build exit={process.returncode}; log={log}")
    if process.returncode:
        raise SystemExit(process.returncode)


def run(name, abundance_mode, force_retry):
    execution = placement()
    target = EXPERIMENT / name
    manifest = json.loads((target / "manifest.json").read_text())
    if manifest["backend"] != "LEGACY" or manifest["candidate_gate"] != "closed":
        raise RuntimeError("This installed-dependency runner is LEGACY-only")
    if digest(target / "seed.mod") != manifest["derived_seed_sha256"]:
        raise RuntimeError("derived seed changed after preparation")
    expected_xa = "".join(f"{isotope} {mass:.17e}\n" for isotope, mass in zip(NAMES, pinned_mixture()))
    if (target / "uniform_xa.data").read_text() != expected_xa:
        raise RuntimeError("runtime pinned-mixture file differs from audited preparation")
    if (target / "frozen20.net").read_text() != "add_isos(\n" + ",\n".join(NAMES) + "\n)\n":
        raise RuntimeError("reaction-free network changed after preparation")
    if (target / "run_command.json").exists():
        raise RuntimeError("use a fresh prepared directory; do not overwrite a run")
    before = inventory()
    write_json(target / "installed_before.json", before)
    runtime_cache = target / "runtime_cache"
    eos_cache = runtime_cache / "eosDT_cache"
    eos_cache.mkdir(parents=True)
    for path in (MESA / "data/eosDT_data/cache").glob("*.bin"):
        shutil.copy2(path, eos_cache / path.name)
    temporary = target / "runtime_temp_cache"
    temporary.mkdir()
    runtime_environment = environment()
    runtime_environment.update(TASK049_BACKEND="LEGACY", TASK049_ABUNDANCE_MODE=str(abundance_mode),
                               TASK049_FORCE_RETRY=str(int(force_retry)),
                               MESA_CACHES_DIR=str(runtime_cache), MESA_TEMP_CACHES_DIR=str(temporary))
    executable = EXPERIMENT / "build/star"
    command = [str(executable)]
    remaining = 1793334000 - time.time() - 120
    timeout = min(900, int(remaining))
    if timeout < 30:
        raise RuntimeError("insufficient time before hard stop")
    write_json(target / "run_command.json", {"execution": execution, "command": command, "cwd": str(target),
               "abundance_mode": abundance_mode, "force_retry": force_retry,
               "executable_sha256": digest(executable), "timeout_seconds": timeout,
               "memory_limit_bytes": 8 * 1024**3, "numerical_threads": 1,
               "input_sha256": {filename: digest(target / filename) for filename in
                   ("seed.mod", "inlist", "frozen20.net", "uniform_xa.data", "history_columns.list", "profile_columns.list")}})
    start = time.monotonic()
    timed_out = False
    with (target / "run.log").open("w") as stream:
        try:
            process = subprocess.run(command, cwd=target, env=runtime_environment, stdout=stream,
                                     stderr=subprocess.STDOUT, timeout=timeout, preexec_fn=limit_memory)
            returncode = process.returncode
        except subprocess.TimeoutExpired:
            timed_out, returncode = True, None
    after = inventory()
    write_json(target / "installed_after.json", after)
    outcome = {"returncode": returncode, "timed_out": timed_out, "seconds": time.monotonic() - start,
               "installed_unchanged": before == after, "completion_written": (target / "completion.txt").exists()}
    write_json(target / "run_outcome.json", outcome)
    print(json.dumps(outcome))
    if before != after:
        raise RuntimeError("installed dependencies or protected inputs changed during run")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subcommands = parser.add_subparsers(dest="action", required=True)
    subcommands.add_parser("build")
    preparation = subcommands.add_parser("prepare")
    preparation.add_argument("name")
    preparation.add_argument("--seed", type=Path, default=SEED)
    execution = subcommands.add_parser("run")
    execution.add_argument("name")
    execution.add_argument("--abundance-mode", type=int, choices=(-1, 0, 1), default=0)
    execution.add_argument("--force-retry", action="store_true")
    arguments = parser.parse_args()
    if arguments.action == "prepare":
        prepare(arguments.name, arguments.seed)
    elif arguments.action == "build":
        build()
    else:
        run(arguments.name, arguments.abundance_mode, arguments.force_retry)


if __name__ == "__main__":
    main()
