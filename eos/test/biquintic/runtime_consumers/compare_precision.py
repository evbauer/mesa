"""Hold polynomial basis and summation order fixed while changing precision only."""

from dataclasses import asdict
import json
import time

import numpy as np

from diagnose_precision import ExtendedCache
from run_consumers import EXPERIMENTS, OUTPUT, LN10, scalar_solver, write_json
from runtime_environment import sha256


def main():
    output = OUTPUT / "precision_matched"
    output.mkdir(parents=True, exist_ok=True)
    failures = json.loads((OUTPUT / "precision_control/frozen_failures.json").read_text())
    failure = next(record for record in failures if record["route"] == "component")
    manifest = json.loads((EXPERIMENTS / "038_solver_manifest/queries_v2.json").read_text())
    case = next(record for record in manifest["cases"] if record["case_id"] == failure["case_id"])
    density_log = case["logRho"] * LN10
    samples = [json.loads(line) for line in (OUTPUT / "precision_control/representable_states.jsonl").read_text().splitlines()]
    samples = [record for record in samples if record["route"] == "component"]
    initial, heating = failure["initial_y"], failure["H_erg_per_g_s"] * failure["dt_seconds"]
    started = time.time()
    statistics, controls, values = {}, {}, {}
    for name, precision in (("binary64", np.float64), ("extended", np.longdouble)):
        evaluator = ExtendedCache(precision)
        residuals = []
        for sample in samples:
            cv, _ = evaluator.evaluate(density_log, sample["y"])
            residuals.append(float(np.longdouble(sample["y"]) - np.longdouble(initial) - np.longdouble(heating) / (
                np.exp(np.longdouble(sample["y"])) * np.longdouble(cv))))
        changes = np.diff(residuals)
        statistics[name] = dict(significand_bits=np.finfo(precision).nmant + 1,
                                max_adjacent_residual_jump=float(max(abs(changes))),
                                decreasing_residual_steps=int(np.count_nonzero(changes < 0)),
                                min_abs_residual=float(min(abs(np.array(residuals)))))
        def callback(temperature):
            return tuple(float(value) for value in evaluator.evaluate(density_log, temperature))
        solved = scalar_solver.solve(callback, initial, heating, *failure["bracket"], failure["solver_start_y"],
                                     tolerance=failure["tolerance_absolute_R"], max_iterations=40)
        controls[name] = asdict(solved)
        values[name] = residuals
    result = dict(trial_id=failure["trial_id"], case_id=case["case_id"], samples=len(samples),
                  scope="Offline diagnostic only; no EOS calls, source changes, or replacement of frozen actual-API failures",
                  held_fixed="Same binary64 cache nodes, effective API coordinates, polynomial basis, loop order, forcing, solver and tolerance",
                  precision_statistics=statistics,
                  controls={name: dict(status=control["status"], iterations=control["iterations"],
                                       last_residual=control["history"][-1]["residual"]) for name, control in controls.items()},
                  elapsed_seconds=time.time() - started,
                  helper_sha256=sha256(__import__("diagnose_precision").__file__))
    write_json(output / "results.json", result)
    write_json(output / "traces.json", controls)
    write_json(output / "samples.json", dict(y=[record["y"] for record in samples], residuals=values))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
