"""Exercise the actual installed EOS final callback with supported and rejected mixtures."""

import argparse
import json
import re
import subprocess

from harness import EXPERIMENT, PATTERN, environment, inventory, limit_memory, placement, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-name', default='guard_checks')
    arguments = parser.parse_args()
    if not re.fullmatch(r'[a-z][a-z0-9_]*', arguments.output_name):
        raise ValueError('use a simple output directory name')
    execution = placement()
    output = EXPERIMENT / arguments.output_name
    output.mkdir(exist_ok=False)
    before = inventory()
    runtime_environment = environment()
    temporary = output / "runtime_temp_cache"
    temporary.mkdir()
    runtime_environment.update(MESA_CACHES_DIR=str(EXPERIMENT / "control_legacy/runtime_cache"),
                               MESA_TEMP_CACHES_DIR=str(temporary))
    outcomes = []
    for name, expected in (("ordered", None), ("reordered", None),
                           ("metal_swap", "TASK049 composition drift"),
                           ("previous_swap", "TASK049 composition drift"),
                           ("retry_swap", "TASK049 composition drift"),
                           ("normalization", "TASK049 nonnormalized abundances"),
                           ("extra_zero", "TASK049 wrong species count")):
        case = output / name
        case.mkdir()
        command = [str(EXPERIMENT / "build/guard_driver"), name,
                   str(EXPERIMENT / "control_legacy/uniform_xa.data"), str(PATTERN)]
        with (case / "run.log").open("w") as stream:
            process = subprocess.run(command, cwd=case, env=runtime_environment, stdout=stream,
                                     stderr=subprocess.STDOUT, timeout=60, preexec_fn=limit_memory)
        log = (case / "run.log").read_text()
        passed = (process.returncode == 0 and "TASK049_GUARD_ACCEPTED" in log) if expected is None else (
            process.returncode != 0 and expected in log)
        relative_error = float(next(line.split()[-1] for line in log.splitlines()
                                    if line.startswith('TASK049_REFERENCE_RELATIVE_ERROR')))
        outcomes.append({"case": name, "command": command, "returncode": process.returncode,
                         "expected_rejection": expected, "reference_relative_error": relative_error,
                         "passed": passed and relative_error <= 5e-13})
    after = inventory()
    result = {"execution": execution, "tests": outcomes, "installed_unchanged": before == after,
              "passed": all(outcome["passed"] for outcome in outcomes) and before == after,
              "limitation": "previous/retry negative inputs exercise the same actual callback with phase labels; real retry/previous lifecycle is checked by stellar runs"}
    write_json(output / "results.json", result)
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
