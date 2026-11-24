"""Read-only actual-MESA selection check for the response-tail slices."""

import argparse
import json
import os
import subprocess
from pathlib import Path

import numpy as np

from freeeos_potential import checksum


def main():
    workspace = Path(__file__).resolve().parents[4]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    mesa = workspace / "mesa"
    executable = mesa / "build/eos-plotter/bin/eos-plotter"
    cache_directory = output / "cache/eosDT_cache"
    cache_directory.mkdir(parents=True, exist_ok=True)
    for original in (mesa / "data/eosDT_data/cache").glob("*.bin"):
        target = cache_directory / original.name
        if not target.exists():
            target.symlink_to(original)
    environment = os.environ.copy()
    sdk = workspace / "mesasdk-x86_64-linux-26.6.1"
    environment.update(MESA_DIR=str(mesa), MESASDK_ROOT=str(sdk), OMP_NUM_THREADS="1",
                       MESA_CACHES_DIR=str(output / "cache"), MESA_TEMP_CACHES_DIR=str(output / "temp_cache"))
    environment["LD_LIBRARY_PATH"] = str(sdk / "lib") + ":" + str(sdk / "lib64") + ":" + environment.get("LD_LIBRARY_PATH", "")
    results = {"purpose": "actual default-MESA fractions, not new potential backend or solver comparison",
               "composition": "plotter H1/He4/C12/O16=.70/.28/.01/.01; not independent scaled-solar direct-source mixture",
               "mesa_revision": subprocess.check_output(["git", "-C", str(mesa), "rev-parse", "HEAD"], text=True).strip(),
               "executable_sha256": checksum(executable), "WDM_AV": False, "slices": {}}
    for name, logtemperature, density_minimum, density_maximum in (
            ("u_face_tail", 6.5117106292763225, 3.283421258552645, 3.783421258552645),
            ("tau_face_tail", 6.6, 3.05, 3.85),
            ("chi_tail", 6.765258774036014, 2.15, 2.75)):
        directory = output / name
        directory.mkdir(exist_ok=True)
        inlist = directory / "inlist_plotter"
        inlist.write_text(f"""&eos
 use_WDM_AV = .false.
/
&plotter
 ignore_ierr = .true.
 export_audit = .true.
 i_eos = 0
 i_eos_other = -1
 only_blend_regions = .false.
 i_var = 11
 doing_partial = .false.
 doing_dfridr = .false.
 doing_d_dlnd = .false.
 doing_consistency = .false.
 i_cons = 2
 xname = 'Rho'
 yname = 'T'
 nT = 1
 nRho = 1001
 nX = 1
 nZ = 1
 X_center = 0.70
 Z_center = 0.02
 delta_X = 0.0
 delta_Z = 0.0
 logT_min = {logtemperature:.17e}
 logT_max = {logtemperature:.17e}
 logRho_min = {density_minimum:.17e}
 logRho_max = {density_maximum:.17e}
/
""")
        with (directory / "run.log").open("w") as stream:
            subprocess.run([str(executable), "inlist_plotter", "audit.dat"], cwd=directory, env=environment,
                           stdout=stream, stderr=subprocess.STDOUT, check=True)
        data = np.loadtxt(directory / "audit.dat", skiprows=4)
        if data.shape != (1001, 85):
            raise ValueError(f"Unexpected plotter export shape {data.shape}")
        fractions = data[:, 3 + 19:3 + 27]
        if not np.all(data[:, -1] == 0):
            raise ValueError("MESA selection query failed")
        target_q = {"u_face_tail": 2.51, "tau_face_tail": 2.2848309909032887, "chi_tail": .7999999999999989}[name]
        target_density = target_q + 2 * logtemperature - 12
        nearest = int(np.argmin(abs(data[:, 1] - target_density)))
        results["slices"][name] = {"calls": len(data), "errors": 0, "logT": logtemperature,
                                   "logRho": [density_minimum, density_maximum],
                                   "FreeEOS_fraction_minmax": [float(fractions[:, 2].min()), float(fractions[:, 2].max())],
                                   "Skye_fraction_minmax": [float(fractions[:, 4].min()), float(fractions[:, 4].max())],
                                   "fraction_columns": ["HELM", "OPAL_SCVH", "FreeEOS", "PC", "Skye", "CMS", "ideal", "WDM_AV"],
                                   "nearest_recorded_extremum": {"target_logQ": target_q, "target_logRho": target_density,
                                                                 "evaluated_logRho": float(data[nearest, 1]),
                                                                 "FreeEOS_fraction": float(fractions[nearest, 2]), "Skye_fraction": float(fractions[nearest, 4])}}
    results["failed_attempts"] = [{"failure": "SDK lib64 missing from LD_LIBRARY_PATH; GFORTRAN_10/GCC_14.0.0 loader versions not found",
                                   "resolution": "add local SDK lib and lib64, rerun all three slices successfully"}]
    (output / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results["slices"], indent=2))


if __name__ == "__main__":
    main()
