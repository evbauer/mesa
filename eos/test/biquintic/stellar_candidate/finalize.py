"""Freeze completed TASK052 evidence and reports without another stellar run."""

import json

from harness import CODE, RUNS, SUMMARY, WORKSPACE, digest, inventory, placement, write_json


def main():
    execution = placement()
    first = json.loads((SUMMARY / "comparison.json").read_text())
    extension = json.loads((SUMMARY / "extension_comparison.json").read_text())
    freeze = json.loads((RUNS / "driver_build/freeze.json").read_text())
    if not first["comparison_completed"] or not extension["comparison_completed"]:
        raise ValueError("unfinished paired comparison")
    if inventory() != freeze["dependencies"] or digest(RUNS / "driver_build/star") != freeze["executable_sha256"]:
        raise ValueError("frozen dependencies/executable changed")
    for filename, expected in freeze["sources"].items():
        if digest(CODE / filename) != expected:
            raise ValueError("frozen source changed: " + filename)
    paths = [path for path in CODE.iterdir() if path.is_file()]
    paths += [RUNS / "driver_build/star", RUNS / "driver_build/freeze.json", RUNS / "driver_build/link.map"]
    for name in ("baseline_00", "candidate_00", "baseline_100", "candidate_100"):
        directory = RUNS / name
        paths += [path for path in directory.iterdir() if path.is_file()]
        for subdirectory in ("LOGS", "photos"):
            paths += [path for path in (directory / subdirectory).rglob("*") if path.is_file()]
    for name in ("comparison.json", "HANDOFF.json", "extension_comparison.json", "report_numbers.json",
                 "extension_report_numbers.json", "profile_comparison.pdf", "extension_profile_comparison.pdf"):
        paths.append(SUMMARY / name)
    reports = [WORKSPACE / "freeeos_biquintic_progress/progress_reports" / ("2026-10-29_stellar_candidate." + suffix)
               for suffix in ("md", "tex", "pdf")]
    paths += reports
    hashes = {str(path.relative_to(WORKSPACE)): digest(path) for path in sorted(set(paths))}
    profile_paths = {backend: RUNS / (backend + "_100") / "LOGS/profile100.data" for backend in ("baseline", "candidate")}
    result = {"execution": execution, "status": "TASK052_FROZEN_COMPLETE",
              "scope": "eight-step and100-step exploratory frozen-composition stellar relaxation, not production/tail-stress acceptance",
              "first_pair_completed": True, "same_seed_100_pair_completed": True,
              "active_stellar_processes_owned_by_task": 0, "further_evolution_authorized_by_this_handoff": False,
              "executable_sha256": freeze["executable_sha256"],
              "seed_sha256": "cd5b74e7ccd55d8253d045abc5bb95cba1506ba42fba70482bccf5f29cf9eb3e",
              "final_profile_age_years": extension["candidate"]["final_age_years"],
              "exact_common_profile_ages": extension["common_age_comparison"]["exact_age_pairs"],
              "interpolated_profile_ages": extension["common_age_comparison"]["interpolated_age_pairs"],
              "final_profiles": {backend: {"path": str(path.relative_to(WORKSPACE)), "sha256": digest(path)}
                                 for backend, path in profile_paths.items()},
              "reports": {str(path.relative_to(WORKSPACE)): digest(path) for path in reports},
              "artifact_hashes": hashes,
              "hash_scope": "all driver source, executable/link/freeze, four run direct files, all LOGS/photos, principal analyses and report sources/PDF; copied runtime-cache files excluded",
              "limitations": ["051retained1e-13stress failures remain failures", "hybrid abundance derivatives are not composition-consistent",
                              "crossing counts cover accepted profile endpoints, not global rough-face stress",
                              "no burning, mesh/time convergence, source-accuracy or production-evolution claim"]}
    write_json(SUMMARY / "FINAL_HANDOFF.json", result)
    print(json.dumps({"status": result["status"], "hashed_artifacts": len(hashes), "reports": result["reports"],
                      "final_profiles": result["final_profiles"]}, indent=2))


if __name__ == "__main__":
    main()
