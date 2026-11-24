"""Analyze stored 070 queries without changing EOS source or old artifacts."""

import json
import resource
import sys

import numpy as np

from common import CODE, OUTPUT, WORKSPACE, digest, placement, save_json, summary

sys.path.insert(0, str(CODE.parent / "helm_control"))
sys.path.insert(0, str(CODE.parent / "pure_freeeos_faces"))
from native import NativeHelm, REAL, basis, decimal_cv
from audit import classify_probe


def route(direction, grid):
    probe = np.loadtxt(OUTPUT / f"routing_{direction}_corrected/probe.tsv", skiprows=1, ndmin=2)
    count = len(grid["din"]) * len(grid["temperature"])
    if probe.shape != (count, 14):
        raise ValueError(f"{direction} probe shape {probe.shape}, expected ({count}, 14)")
    expected_density = np.repeat(grid["din"] / float(np.loadtxt(OUTPUT / "composition.tsv", skiprows=1)[5]),
                                 len(grid["temperature"]))
    expected_temperature = np.tile(grid["temperature"], len(grid["din"]))
    if np.max(np.abs(probe[:, 0] - np.log10(expected_density))) > 1e-12 or \
       np.max(np.abs(probe[:, 1] - np.log10(expected_temperature))) > 1e-12:
        raise ValueError("probe coordinates do not match native-face grid")
    successful, pure, good = classify_probe(probe)
    return successful, pure, good, probe


def density_interpolated_jets(table, density_indices, temperature_index):
    width = table.density[density_indices + 1] - table.density[density_indices]
    jets = []
    for slots in ((0, 1, 3), (2, 5, 6), (4, 7, 8)):
        lower = [np.asarray(table.fields[slot][density_indices, temperature_index], dtype=REAL)
                 for slot in slots]
        upper = [np.asarray(table.fields[slot][density_indices + 1, temperature_index], dtype=REAL)
                 for slot in slots]
        jets.append(REAL("0.5") * (lower[0] + upper[0]) + REAL("0.15625") * width *
                    (lower[1] - upper[1]) + REAL("0.015625") * width**2 *
                    (lower[2] + upper[2]))
    return jets


def temperature_endpoints(table, density_indices, temperature_index, ye):
    previous = density_interpolated_jets(table, density_indices, temperature_index - 1)
    current = density_interpolated_jets(table, density_indices, temperature_index)
    following = density_interpolated_jets(table, density_indices, temperature_index + 1)
    left_width = table.temperature[temperature_index] - table.temperature[temperature_index - 1]
    right_width = table.temperature[temperature_index + 1] - table.temperature[temperature_index]
    left_weights = basis(REAL(1), left_width, 3)
    right_weights = basis(REAL(0), right_width, 3)
    third_left = sum(left_weights[0, slot] * previous[slot] + left_weights[1, slot] * current[slot]
                     for slot in range(3))
    third_right = sum(right_weights[0, slot] * current[slot] + right_weights[1, slot] * following[slot]
                      for slot in range(3))
    temperature = table.temperature[temperature_index]
    electron_cv = -REAL(ye) * temperature * current[2]
    electron_left = -REAL(ye) * temperature * (current[2] + temperature * third_left)
    electron_right = -REAL(ye) * temperature * (current[2] + temperature * third_right)
    electron_delta = -REAL(ye) * temperature**2 * (third_left - third_right)
    return electron_cv, electron_left, electron_right, electron_delta


def temperature_analysis(table, grid, ye):
    successful, pure, good, probe = route("T", grid)
    kernel = np.loadtxt(OUTPUT / "helm_T_kernel.tsv", skiprows=1, ndmin=2)
    count = len(probe)
    if kernel.shape != (count, 14):
        raise ValueError(f"HELM kernel shape {kernel.shape}, expected ({count}, 14)")
    if np.max(np.abs(np.log10(kernel[:, 0]) - probe[:, 0])) > 1e-12 or \
       np.max(np.abs(np.log10(kernel[:, 1]) - probe[:, 1])) > 1e-12:
        raise ValueError("HELM kernel coordinates do not match probe")
    shape = (len(grid["din"]), len(grid["temperature"]))
    traces = {name: np.full(shape, np.nan) for name in
              ("Cv_left", "Cv_right", "Cv_t_left", "Cv_t_right", "electron_Cv",
               "electron_Cv_t_left", "electron_Cv_t_right", "electron_delta_Cv_t",
               "smooth_Cv", "smooth_Cv_t", "denominator", "jump", "precision_floor")}
    kernel_good = (kernel[:, 2] == 0) & np.isfinite(kernel).all(axis=1)
    kernel_grid = kernel.reshape(shape + (14,))
    good_grid = kernel_good.reshape(shape)
    for column, temperature_index in enumerate(grid["temperature_index"]):
        electron_cv, electron_left, electron_right, electron_delta = temperature_endpoints(
            table, grid["density_index"], int(temperature_index), ye)
        rows = kernel_grid[:, column]
        valid = good_grid[:, column]
        smooth_cv = np.asarray(rows[:, 5] + rows[:, 7] + rows[:, 9], dtype=REAL)
        smooth_t = np.asarray(rows[:, 6] + rows[:, 8] + rows[:, 10], dtype=REAL)
        total_cv = electron_cv + smooth_cv
        total_left = electron_left + smooth_t
        total_right = electron_right + smooth_t
        denominator = np.maximum.reduce((abs(total_cv), abs(total_left), abs(total_right)))
        jump = np.divide(electron_delta, denominator, out=np.zeros(len(denominator), dtype=REAL),
                         where=denominator != 0)
        floor = REAL("1e-12") * np.maximum.reduce((abs(electron_cv), abs(electron_left), abs(electron_right)))
        values = {"Cv_left": total_cv, "Cv_right": total_cv,
                  "Cv_t_left": total_left, "Cv_t_right": total_right,
                  "electron_Cv": electron_cv, "electron_Cv_t_left": electron_left,
                  "electron_Cv_t_right": electron_right, "electron_delta_Cv_t": electron_delta,
                  "smooth_Cv": smooth_cv, "smooth_Cv_t": smooth_t,
                  "denominator": denominator, "jump": jump, "precision_floor": floor}
        for name, value in values.items():
            traces[name][valid, column] = np.asarray(value[valid], dtype=float)
    finite = np.isfinite(traces["jump"].ravel())
    external = successful & good & ~pure & kernel_good
    delta = np.abs(traces["electron_delta_Cv_t"].ravel())
    floor = traces["precision_floor"].ravel()
    resolved = finite & (delta > floor)
    gamma = kernel[:, 13]
    branch_near = kernel_good & (np.abs(gamma - 1) < 1e-6)
    shown = external & resolved & ~branch_near
    archive = {key: value for key, value in traces.items()}
    archive.update({"route_success": successful.reshape(shape), "pure": pure.reshape(shape),
                    "fraction_good": good.reshape(shape), "kernel_good": kernel_good.reshape(shape),
                    "resolved": resolved.reshape(shape), "branch_near": branch_near.reshape(shape),
                    "shown": shown.reshape(shape), "Gamma": gamma.reshape(shape)})
    np.savez_compressed(OUTPUT / "helm_T_traces_corrected.npz", **archive)
    selected_denom = traces["denominator"].ravel()[shown]
    result = {"native_faces": count, "routing_success": int(successful.sum()),
              "routing_bad_fractions": int((successful & ~good).sum()), "pure_FreeEOS": int(pure.sum()),
              "kernel_success": int(kernel_good.sum()), "kernel_failure": int((~kernel_good).sum()),
              "external": int(external.sum()), "precision_unresolved_external": int((external & ~resolved).sum()),
              "coulomb_branch_near_external": int((external & branch_near).sum()),
              "shown": int(shown.sum()), "minimum_denominator": float(np.min(selected_denom)) if len(selected_denom) else None,
              "minimum_abs_Gamma_minus_one": float(np.min(np.abs(gamma[kernel_good] - 1))) if kernel_good.any() else None,
              "shown_jump": summary(traces["jump"].ravel()[shown]),
              "kernel_sha256": digest(OUTPUT / "helm_T_kernel.tsv"),
              "probe_sha256": digest(OUTPUT / "routing_T_corrected/probe.tsv")}
    return result, shown, traces


def density_analysis(grid):
    successful, pure, good, probe = route("d", grid)
    kernel = np.loadtxt(OUTPUT / "helm_d_kernel.tsv", skiprows=1, ndmin=2)
    if kernel.shape != (len(probe), 14):
        raise ValueError("HELM density kernel shape does not match direct routing probe")
    if np.max(np.abs(np.log10(kernel[:, 0]) - probe[:, 0])) > 1e-12 or \
       np.max(np.abs(np.log10(kernel[:, 1]) - probe[:, 1])) > 1e-12:
        raise ValueError("HELM density kernel coordinates do not match probe")
    shape = (len(grid["din"]), len(grid["temperature"]))
    kernel_good = (kernel[:, 2] == 0) & np.isfinite(kernel).all(axis=1)
    branch_near = kernel_good & (np.abs(kernel[:, 13] - 1) < 1e-6)
    external = successful & good & ~pure & kernel_good
    shown = external & ~branch_near
    denominator = np.maximum(abs(kernel[:, 11]), abs(kernel[:, 12]))
    np.savez_compressed(OUTPUT / "helm_d_traces.npz", route_success=successful.reshape(shape),
                        pure=pure.reshape(shape), fraction_good=good.reshape(shape),
                        kernel_good=kernel_good.reshape(shape), branch_near=branch_near.reshape(shape),
                        shown=shown.reshape(shape), jump=np.zeros(shape),
                        Cv_left=kernel[:, 11].reshape(shape), Cv_right=kernel[:, 11].reshape(shape),
                        Cv_t_left=kernel[:, 12].reshape(shape), Cv_t_right=kernel[:, 12].reshape(shape),
                        denominator=denominator.reshape(shape), Gamma=kernel[:, 13].reshape(shape))
    return {"native_faces": len(probe), "routing_success": int(successful.sum()),
            "routing_bad_fractions": int((successful & ~good).sum()),
            "pure_FreeEOS": int(pure.sum()), "kernel_success": int(kernel_good.sum()),
            "kernel_failure": int((~kernel_good).sum()), "external_zero_faces": int(shown.sum()),
            "coulomb_branch_near_external": int((external & branch_near).sum()),
            "minimum_denominator": float(np.min(denominator[shown])) if shown.any() else None,
            "minimum_abs_Gamma_minus_one": float(np.min(np.abs(kernel[kernel_good, 13] - 1))) if kernel_good.any() else None,
            "zero_derivation": "Cv_t has identical one-sided limits at native d faces for the tensor Hermite electron table; ion/radiation and same-branch Coulomb terms are smooth",
            "probe_sha256": digest(OUTPUT / "routing_d_corrected/probe.tsv"),
            "kernel_sha256": digest(OUTPUT / "helm_d_kernel.tsv")}


def decimal_control(table, grid, ye, shown, traces):
    positions = np.flatnonzero(shown)
    if not len(positions):
        return {"count": 0}
    selected = np.unique(np.concatenate((positions[np.linspace(0, len(positions) - 1, 24, dtype=int)],
                                         positions[np.argsort(np.abs(traces["jump"].ravel()[positions]))[-8:]])))
    records = []
    columns = len(grid["temperature"])
    for position in selected:
        density_column, temperature_column = divmod(int(position), columns)
        density_index = int(grid["density_index"][density_column])
        temperature_index = int(grid["temperature_index"][temperature_column])
        din = float(grid["din"][density_column])
        temperature = float(grid["temperature"][temperature_column])
        left = decimal_cv(table, din, temperature, ye, (density_index, temperature_index - 1))
        right = decimal_cv(table, din, temperature, ye, (density_index, temperature_index))
        decimal_difference = float(left[2] - right[2])
        calculated = float(traces["electron_delta_Cv_t"].ravel()[position])
        denominator = float(traces["denominator"].ravel()[position])
        records.append({"flat_index": int(position), "density_index": density_index,
                        "temperature_index": temperature_index, "logRho": float(np.log10(din / ye)),
                        "logT": float(np.log10(temperature)), "decimal_difference": decimal_difference,
                        "longdouble_difference": calculated,
                        "absolute_disagreement_over_total_scale": abs(decimal_difference - calculated) / denominator})
    save_json(OUTPUT / "decimal_control_corrected.json", records)
    return {"count": len(records), "maximum_disagreement_over_total_scale":
            max(item["absolute_disagreement_over_total_scale"] for item in records)}


def main():
    placement()
    resource.setrlimit(resource.RLIMIT_AS, (12 * 1024**3, 12 * 1024**3))
    ye = float(np.loadtxt(OUTPUT / "composition.tsv", skiprows=1)[5])
    table = NativeHelm(WORKSPACE / "mesa/data/eosDT_data/cache/helm_table.bin",
                       OUTPUT / "density_axis.txt", OUTPUT / "temperature_axis.txt")
    with np.load(OUTPUT / "helm_T_grid.npz") as grid:
        temperature_result, shown, traces = temperature_analysis(table, grid, ye)
        validation = decimal_control(table, grid, ye, shown, traces)
    if (OUTPUT / "helm_d_traces.npz").exists():
        density_result = json.loads((OUTPUT / "helm_summary.json").read_text())["density_faces"]
    else:
        with np.load(OUTPUT / "helm_d_grid.npz") as grid:
            density_result = density_analysis(grid)
    report = {"temperature_faces": temperature_result, "density_faces": density_result,
              "decimal_validation": validation, "normalization": "max(abs(Cv_L),abs(Cv_R),abs(Cv_t_L),abs(Cv_t_R)); all-zero returns 0",
              "display_floor": 1e-17, "coordinate_system": "native HELM linear d=Ye*rho and T; log axes only for display"}
    save_json(OUTPUT / "helm_summary_corrected.json", report)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
