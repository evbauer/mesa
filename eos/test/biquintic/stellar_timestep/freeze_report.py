"""Freeze057completed exploratory continuation while preserving original failures."""

import json

import timestep as task


def main():
    execution = task.frozen.placement()
    verified = task.verify_frozen(full=True)
    comparison = json.loads((task.ANALYSIS / "comparison.json").read_text())
    if comparison["status"] != "TASK057_EXPLORATORY_CONTINUATION_COMPLETE_WITH_RETAINED_ORIGINAL_GATE_FAILURE" or comparison["original_gate_passed"]:
        raise ValueError("wrong final scientific disposition")
    pinned = {"STOP.json": "7163096f1ded72d1cba58b080e8fb6b86cd509064625d0e405bdb6381e3778d5",
              "analysis/coarse_candidate.json": "9a03970330e1034b9fb1d9f53361c773b1e97fe7b54c4323e80662f5d53841fd",
              "CONTINUATION_STOP_half_baseline.json": "36ec3412d58530e412ac3f8f48fd5e5aa8453919e347bd68bfdaef3257791f08"}
    for relative, expected in pinned.items():
        if task.frozen.digest(task.EXPERIMENT / relative) != expected:
            raise ValueError("original failure evidence changed")
    paths = [path for path in task.CODE.iterdir() if path.is_file()]
    paths += [path for path in task.EXPERIMENT.iterdir() if path.is_file() and path.name not in ("HANDOFF.json", "freeze_stdout.log", "CHECKPOINT.md")]
    paths += [path for path in task.ANALYSIS.iterdir() if path.is_file()]
    for name in task.NAMES:
        directory = task.RUNS / name
        paths += [path for path in directory.iterdir() if path.is_file()]
        for subdirectory in ("LOGS", "photos"):
            paths += [path for path in (directory / subdirectory).rglob("*") if path.is_file()]
    reports = [task.frozen.WORKSPACE / "freeeos_biquintic_progress/progress_reports" / ("2026-10-29_stellar_timestep." + suffix)
               for suffix in ("md", "tex", "pdf")]
    paths += reports
    hashes = {str(path.relative_to(task.frozen.WORKSPACE)): task.frozen.digest(path) for path in sorted(set(paths))}
    result = {"execution": execution, "status": comparison["status"], "original_gate_passed": False,
              "exploratory_continuation_completed": True, "original_failure_hashes": pinned,
              "verified_unchanged052artifacts": verified, "owned_stellar_processes_running": 0,
              "endpoint_ages": comparison["final_ages"], "endpoint_interpolation": False,
              "final_profiles": {name: {"path": result["final_profile"], "sha256": result["final_profile_sha256"],
                                         "accepted_steps": result["accepted_steps"], "offset_ulps": result["endpoint"]["offset_ulps"]}
                                 for name, result in comparison["runs"].items()},
              "reports": {str(path.relative_to(task.frozen.WORKSPACE)): task.frozen.digest(path) for path in reports},
              "artifact_hashes": hashes,
              "hash_scope": "source, stable metadata/authorizations/failures, all four retained run outputs/profiles/photos, analyses and report formats; mutable checkpoint, freeze stdout and runtime cache copies excluded",
              "limitations": ["original5e-11coarse gate remains failed", "half profiles2through89not retained by unchanged MESA policy",
                              "endpoint sensitivity is empirical, not a rigorous error bound", "two resolutions do not prove convergence order",
                              "all prior fixed-composition/hybrid/support and051stress qualifications remain", "no additional run authorized"]}
    task.frozen.write_json(task.EXPERIMENT / "HANDOFF.json", result)
    print(json.dumps({"status": result["status"], "original_gate_passed": False, "hashed_artifacts": len(hashes),
                      "reports": result["reports"], "final_profiles": result["final_profiles"]}, indent=2))


if __name__ == "__main__":
    main()
