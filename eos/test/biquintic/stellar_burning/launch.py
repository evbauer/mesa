"""Detach long TASK067 runs from the short-lived CLI command session."""

import argparse
import json
from pathlib import Path
import subprocess
import sys
import time


CODE = Path(__file__).resolve().parent
WORKSPACE = CODE.parents[4]
EXPERIMENT = WORKSPACE / "freeeos_biquintic_progress/experiments/067_stellar_burning"


def run_dir(backend, pair):
    return EXPERIMENT / f"{backend}_{pair}"


def supervise(backend, pair):
    destination = run_dir(backend, pair)
    command = [sys.executable, str(CODE / "harness.py"), "run", backend, "--pair", pair]
    started = time.time()
    with (destination / "launcher.log").open("w") as stream:
        result = subprocess.run(command, stdin=subprocess.DEVNULL,
                                stdout=stream, stderr=subprocess.STDOUT)
    (destination / "launcher_outcome.json").write_text(json.dumps({
        "command": command, "returncode": result.returncode,
        "started": started, "finished": time.time(),
    }, indent=2) + "\n")


def launch(backend, pair):
    destination = run_dir(backend, pair)
    if not (destination / "manifest.json").is_file():
        raise RuntimeError("prepared run directory absent")
    if (destination / "launcher.json").exists() or (destination / "command.json").exists():
        raise RuntimeError("attempt already launched; preserve it")
    command = [sys.executable, str(__file__), "supervise", backend, "--pair", pair]
    with (destination / "launcher_boot.log").open("w") as stream:
        process = subprocess.Popen(command, stdin=subprocess.DEVNULL,
                                   stdout=stream, stderr=subprocess.STDOUT,
                                   start_new_session=True)
    (destination / "launcher.json").write_text(json.dumps({
        "command": command, "pid": process.pid, "started": time.time(),
        "detached_session": True,
    }, indent=2) + "\n")
    print(json.dumps({"pid": process.pid, "backend": backend, "pair": pair}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("launch", "supervise"))
    parser.add_argument("backend", choices=("baseline", "candidate"))
    parser.add_argument("--pair", required=True)
    args = parser.parse_args()
    if args.action == "launch":
        launch(args.backend, args.pair)
    else:
        supervise(args.backend, args.pair)


if __name__ == "__main__":
    main()
