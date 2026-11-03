"""Hash final TASK070 source, inputs, controls, and deliverables."""

import os
import socket
from datetime import datetime
from zoneinfo import ZoneInfo

from common import CODE, OLD, OUTPUT, WORKSPACE, digest, save_json


def main():
    source_paths = [
        WORKSPACE / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data",
        WORKSPACE / "mesa/data/eosDT_data/helm_table.dat",
        WORKSPACE / "mesa/data/eosDT_data/cache/helm_table.bin",
        WORKSPACE / "mesa/eos/eosFreeEOS_builder/mass_frac_solar.txt",
        WORKSPACE / "mesa/build/eos/lib/libeos.a",
        WORKSPACE / "mesa/eos/eosFreeEOS_builder/biquintic/mesa_probe/build/mesa_probe",
        OLD / "face_manifest.npz", OLD / "face_results.npz", OLD / "face_probe.tsv",
        WORKSPACE / "freeeos_biquintic_progress/experiments/048_helm_control/results.json",
    ]
    code_paths = [CODE / name for name in ("common.py", "run.py", "analyze.py", "plot.py",
                                           "display.py", "norm_audit.py", "validate.py",
                                           "diagnostics.py", "provenance.py", "test_map.py",
                                           "helm_probe.f90", "Makefile", "README.md")]
    artifact_paths = [OUTPUT / name for name in (
        "build/helm_map_probe", "bounds.json", "input_hashes.json", "commands.txt",
        "freeeos_traces.npz", "freeeos_summary.json",
        "helm_T_grid.npz", "helm_d_grid.npz", "helm_T_traces_corrected.npz",
        "helm_d_traces.npz", "helm_summary_corrected.json", "coverage_diagnostics.json",
        "decimal_control_corrected.json", "validation/summary.json",
        "display_bounds_corrected.json", "in_display_stats_corrected.json",
        "normalization_audit.json", "freeeos_helm_native_faces_corrected.png",
        "freeeos_helm_native_faces_corrected.pdf",
        "freeeos_helm_native_faces_corrected_presentation.json",
        "routing_T_corrected/probe.tsv", "routing_d_corrected/probe.tsv",
        "helm_T_kernel.tsv", "helm_d_kernel.tsv")]
    report = WORKSPACE / "freeeos_biquintic_progress/progress_reports"
    report_paths = [report / f"2026-11-02_freeeos_helm_map.{extension}" for extension in ("md", "tex", "pdf")]
    result = {"created_PST": datetime.now(ZoneInfo("America/Los_Angeles")).isoformat(),
              "job": os.environ.get("SLURM_JOB_ID"), "host": socket.gethostname(),
              "sources": {str(path.relative_to(WORKSPACE)): digest(path) for path in source_paths},
              "code": {str(path.relative_to(WORKSPACE)): digest(path) for path in code_paths},
              "artifacts": {str(path.relative_to(WORKSPACE)): digest(path) for path in artifact_paths},
              "report": {str(path.relative_to(WORKSPACE)): digest(path) for path in report_paths}}
    output = OUTPUT / "provenance_final.json"
    if output.exists():
        raise FileExistsError("provenance_final.json is immutable")
    save_json(output, result)
    print(f"Saved {len(result['sources'])} sources, {len(result['code'])} code files, "
          f"{len(result['artifacts'])} artifacts, {len(result['report'])} report files")


if __name__ == "__main__":
    main()
