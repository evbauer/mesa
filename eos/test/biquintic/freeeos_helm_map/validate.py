"""Independent native-kernel offsets and archived 048 overlap controls."""

import json
import resource
import sys

import numpy as np

from common import CODE, OUTPUT, WORKSPACE, digest, placement, save_json, scaled_jump
from run import environment, run_logged

sys.path.insert(0, str(CODE.parent / "helm_control"))
from native import NativeHelm, decimal_cv


OLD048 = WORKSPACE / "freeeos_biquintic_progress/experiments/048_helm_control/results.json"
EPSILONS = (1e-4, 1e-6)


def selected_map_faces(grid, traces):
    shown = traces["shown"].ravel()
    positions = np.flatnonzero(shown)
    if not len(positions):
        return []
    metrics = np.abs(traces["jump"].ravel()[positions])
    chosen = np.unique(np.concatenate((positions[np.linspace(0, len(positions) - 1, 20, dtype=int)],
                                       positions[np.argsort(metrics)[-8:]])))
    width = len(grid["temperature"])
    return [(int(position), *divmod(int(position), width)) for position in chosen]


def queries(grid, traces, table, ye, frozen):
    records = []
    for position, density_column, temperature_column in selected_map_faces(grid, traces):
        density_index = int(grid["density_index"][density_column])
        temperature_index = int(grid["temperature_index"][temperature_column])
        density = float(grid["din"][density_column])
        face_temperature = float(grid["temperature"][temperature_column])
        for epsilon in EPSILONS:
            for side, sign, span in (("left", -1, table.temperature[temperature_index] - table.temperature[temperature_index - 1]),
                                     ("right", 1, table.temperature[temperature_index + 1] - table.temperature[temperature_index])):
                records.append({"kind": "map_offset", "flat_index": position, "side": side,
                                "epsilon_fraction_linear_cell": epsilon, "rho": density / ye,
                                "T": float(face_temperature + sign * epsilon * span)})
    for index, face in enumerate(frozen["faces"]):
        records.append({"kind": "048_center", "face_index": index,
                        "rho": face["rho"], "T": face["T"]})
    return records


def main():
    placement()
    resource.setrlimit(resource.RLIMIT_AS, (12 * 1024**3, 12 * 1024**3))
    frozen = json.loads(OLD048.read_text())
    ye = float(np.loadtxt(OUTPUT / "composition.tsv", skiprows=1)[5])
    table = NativeHelm(WORKSPACE / "mesa/data/eosDT_data/cache/helm_table.bin",
                       OUTPUT / "density_axis.txt", OUTPUT / "temperature_axis.txt")
    with np.load(OUTPUT / "helm_T_grid.npz") as grid, np.load(OUTPUT / "helm_T_traces_corrected.npz") as traces:
        records = queries(grid, traces, table, ye, frozen)
        validation_dir = OUTPUT / "validation"
        validation_dir.mkdir(exist_ok=True)
        query_path = validation_dir / "queries.txt"
        with query_path.open("x") as stream:
            for record in records:
                stream.write(f"{record['rho']:.17e} {record['T']:.17e}\n")
        command = [str(OUTPUT / "build/helm_map_probe"), str(WORKSPACE / "mesa"),
                   str(WORKSPACE / "mesa/eos/eosFreeEOS_builder/mass_frac_solar.txt"), str(validation_dir)]
        run_logged(command, query_path, validation_dir / "kernel.tsv", validation_dir / "kernel.stderr",
                   validation_dir, environment(validation_dir))
        kernel = np.loadtxt(validation_dir / "kernel.tsv", skiprows=1, ndmin=2)
        if kernel.shape != (len(records), 14) or np.any(kernel[:, 2] != 0) or not np.isfinite(kernel).all():
            raise ValueError("validation kernel failed")
        offsets = []
        overlap = []
        for record, row in zip(records, kernel):
            if record["kind"] == "map_offset":
                position = record["flat_index"]
                side = record["side"]
                denominator = float(traces["denominator"].ravel()[position])
                exact_cv = float(traces[f"Cv_{side}"].ravel()[position])
                exact_t = float(traces[f"Cv_t_{side}"].ravel()[position])
                smooth_cv = float(traces["smooth_Cv"].ravel()[position])
                smooth_t = float(traces["smooth_Cv_t"].ravel()[position])
                offsets.append({"flat_index": position, "side": side,
                                "epsilon_fraction_linear_cell": record["epsilon_fraction_linear_cell"],
                                "total_limit_error_over_scale": float(max(abs(row[11] - exact_cv),
                                                                          abs(row[12] - exact_t)) / denominator),
                                "smooth_limit_error_over_scale": float(max(abs(row[5] + row[7] + row[9] - smooth_cv),
                                                                           abs(row[6] + row[8] + row[10] - smooth_t)) / denominator)})
                continue
            face = frozen["faces"][record["face_index"]]
            additions = np.asarray(row[[5, 6]] + row[[7, 8]] + row[[9, 10]])
            old_components = face["center_components"]
            old_smooth = np.asarray(old_components["ion"] + old_components["radiation"] +
                                    old_components["coulomb_by_subtraction"]).reshape(3, 3).sum(axis=0)[[0, 2]]
            scale = max(1, abs(face["total_left"][0]), abs(face["total_left"][2]), abs(face["total_right"][2]))
            if face["direction"] == "T":
                left_cell, right_cell = tuple(face["left_cell"]), tuple(face["right_cell"])
                electron_left = decimal_cv(table, face["din"], face["T"], ye, left_cell)
                electron_right = decimal_cv(table, face["din"], face["T"], ye, right_cell)
                calculated_jump, _ = scaled_jump(float(electron_left[0]) + old_smooth[0],
                                                 float(electron_right[0]) + old_smooth[0],
                                                 float(electron_left[2]) + old_smooth[1],
                                                 float(electron_right[2]) + old_smooth[1])
                archived_jump = face["total_jump"]["Cv_t"]["signed_scaled"]
                jump_difference = abs(float(calculated_jump) - archived_jump)
            else:
                jump_difference = 0.0
            overlap.append({"face_id": face["face_id"], "direction": face["direction"],
                            "smooth_addition_error_over_scale": float(np.max(np.abs(additions - old_smooth)) / scale),
                            "recomputed_jump_absolute_difference": jump_difference,
                            "old_denominator_floor_inactive": face["total_jump"]["Cv_t"]["scale"] > 1})
    report = {"map_offset_count": len(offsets), "overlap048_faces": len(overlap),
              "offset_errors": {str(epsilon): {"maximum_total_limit_error_over_scale":
                  max(item["total_limit_error_over_scale"] for item in offsets if item["epsilon_fraction_linear_cell"] == epsilon),
                  "maximum_smooth_limit_error_over_scale":
                  max(item["smooth_limit_error_over_scale"] for item in offsets if item["epsilon_fraction_linear_cell"] == epsilon)}
                  for epsilon in EPSILONS},
              "overlap048_max_smooth_error_over_scale": max(item["smooth_addition_error_over_scale"] for item in overlap),
              "overlap048_max_jump_difference": max(item["recomputed_jump_absolute_difference"] for item in overlap),
              "overlap048_all_denominator_floors_inactive": all(item["old_denominator_floor_inactive"] for item in overlap),
              "old048_sha256": digest(OLD048), "validation_kernel_sha256": digest(validation_dir / "kernel.tsv")}
    save_json(validation_dir / "offsets.json", offsets)
    save_json(validation_dir / "overlap048.json", overlap)
    save_json(validation_dir / "summary.json", report)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
