"""Run standalone regression suites with source/cache protection and provenance."""

import argparse
import datetime
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path


CODE = Path(__file__).resolve().parent
WORKSPACE = CODE.parents[3]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def code_inventory():
    selected = {}
    for path in sorted(CODE.rglob("*")):
        if not path.is_file() or any(part in ("build", "__pycache__", ".git") for part in path.relative_to(CODE).parts):
            continue
        if path.suffix in (".py", ".f90", ".md", ".sh", ".txt") or path.name in ("Makefile", ".gitignore", "inlist_probe", "inlist_benchmark"):
            selected[str(path.relative_to(WORKSPACE))] = sha256(path)
    return selected


def protected_inventory():
    source = WORKSPACE / "mesa/data/eosFreeEOS_data"
    tables = {str(path.relative_to(WORKSPACE)): sha256(path) for path in sorted(source.glob("mesa-FreeEOS_*.data"))}
    cache = WORKSPACE / "mesa/data/eosDT_data/cache"
    caches = {path.name: {"size": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns}
              for path in sorted(cache.glob("*.bin"))}
    plotter = {str(path.relative_to(WORKSPACE)): sha256(path) for path in
               (WORKSPACE / "mesa/eos/plotter/README.rst", WORKSPACE / "mesa/eos/plotter/src/eos_plotter.f90")}
    return {"source_table_sha256": tables, "installed_cache_size_mtime": caches,
            "preserved_plotter_sha256": plotter}


def suites(skip_fortran):
    command_list = []
    directories = [CODE, *sorted(path for path in CODE.iterdir() if path.is_dir() and any(path.glob("test_*.py"))
                                 and path.name not in ("mesa_probe", "table_fortran"))]
    for directory in directories:
        command_list.append(("python_" + (directory.name if directory != CODE else "lead"),
                             [sys.executable, "-B", "-m", "unittest", "discover", "-s", str(directory), "-p", "test_*.py", "-v"], WORKSPACE))
    if not skip_fortran:
        for mode, directory in (("optimized", CODE / "test_fortran"), ("debug", CODE / "test_fortran/build/debug")):
            command_list.append(("cell_" + mode, [str(directory / "test_biquintic")], WORKSPACE))
            command_list.append(("cell_protocol_" + mode, ["bash", str(CODE / "test_fortran/test_driver.sh"),
                                 str(directory / "cell_driver"), str(CODE / f"test_fortran/build/{mode}/final_validation")], CODE / "test_fortran"))
        for mode, directory in (("optimized", CODE / "table_fortran/build"), ("debug", CODE / "table_fortran/build/debug")):
            command_list.append(("table_" + mode, [str(directory / "test_table"), str(directory)], WORKSPACE))
            command_list.append(("table_protocol_" + mode, [sys.executable, "-B", str(CODE / "table_fortran/test_driver.py"),
                                 str(directory / "table_query"), str(directory)], WORKSPACE))
        command_list.append(("actual_MESA_probe", [sys.executable, "-B", str(CODE / "mesa_probe/test_smoke.py")], WORKSPACE))
    return command_list


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=WORKSPACE / "freeeos_biquintic_progress/experiments/022_final_validation")
    parser.add_argument("--skip-fortran", action="store_true")
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    environment = dict(os.environ)
    sdk = WORKSPACE / "mesasdk-x86_64-linux-26.6.1"
    environment.update(MESASDK_ROOT=str(sdk), OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                       PYTHONDONTWRITEBYTECODE="1", LD_LIBRARY_PATH=f"{sdk}/lib:{sdk}/lib64:" + environment.get("LD_LIBRARY_PATH", ""))
    before_code, before_protected = code_inventory(), protected_inventory()
    started = time.monotonic()
    result = {"started_at": datetime.datetime.now().astimezone().isoformat(),
              "scope": "Standalone regression and MESA probe smoke only; no stellar validation, installed table generation or dependency rebuild",
              "command": sys.argv, "python": sys.executable, "sdk": str(sdk),
              "mesa_commit": subprocess.check_output(["git", "-C", str(WORKSPACE / "mesa"), "rev-parse", "HEAD"], text=True).strip(),
              "mesa_worktree_status": subprocess.check_output(["git", "-C", str(WORKSPACE / "mesa"), "status", "--short"], text=True),
              "code_sha256": before_code, "protected_before": before_protected, "suites": []}
    for name, command, working_directory in suites(arguments.skip_fortran):
        suite_started = time.monotonic()
        try:
            completed = subprocess.run(command, cwd=working_directory, env=environment, capture_output=True, text=True, timeout=300)
            text = completed.stdout + completed.stderr
            returncode = completed.returncode
        except (OSError, subprocess.TimeoutExpired) as error:
            text, returncode = str(error), -1
        log = output / (name + ".txt")
        log.write_text(text)
        unit_count = re.search(r"Ran (\d+) tests? in", text)
        result["suites"].append({"name": name, "command": command, "cwd": str(working_directory),
                                 "returncode": returncode, "elapsed_seconds": time.monotonic() - suite_started,
                                 "unittest_count": int(unit_count.group(1)) if unit_count else None,
                                 "log": str(log.relative_to(WORKSPACE)), "log_sha256": sha256(log)})
        print(f"{name}: {'PASS' if returncode == 0 else 'FAIL'}", flush=True)
    after_code, after_protected = code_inventory(), protected_inventory()
    result.update(elapsed_seconds=time.monotonic() - started, code_inventory_unchanged=before_code == after_code,
                  protected_inventory_unchanged=before_protected == after_protected, protected_after=after_protected,
                  unittest_total=sum(suite["unittest_count"] or 0 for suite in result["suites"]))
    result["passed"] = all(suite["returncode"] == 0 for suite in result["suites"]) and result["code_inventory_unchanged"] and result["protected_inventory_unchanged"]
    (output / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"{'PASS' if result['passed'] else 'FAIL'}: {len(result['suites'])} suites; {result['unittest_total']} unittest cases; protected files unchanged={result['protected_inventory_unchanged']}")
    raise SystemExit(0 if result["passed"] else 1)


if __name__ == "__main__":
    main()
