"""Replay Figure 2 from saved 036 fits and add returned-entropy consistency."""

import argparse
import json
import os
from pathlib import Path
import sys

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"
sys.dont_write_bytecode = True

import numpy as np

CODE = Path(__file__).resolve().parents[1]
WORKSPACE = CODE.parents[3]
REPORT = WORKSPACE / "freeeos_biquintic_progress"
ARCHIVE = REPORT / "experiments/036_dense_derivatives"
OUTPUT = REPORT / "experiments/068_figure2_consistency"
sys.path.insert(0, str(CODE))
sys.path.insert(0, str(WORKSPACE / "freeeos_report"))
import audit_freeeos as legacy
from freeeos_potential import (EVAL_ORDERS, LN10, checksum, construct_twists, evaluate,
                               query_cells, read_table, source_constraints,
                               thermodynamic_partials, to_physical)
from smooth_patch import evaluate_fit

FLOOR = 1e-16
NAMES = ("biquintic", "native_1.0", "thermo_0.1", "legacy")
COLORS = {"biquintic": "tab:blue", "legacy": "tab:orange",
          "native_1.0": "tab:green", "thermo_0.1": "tab:red"}
LABELS = {"biquintic": "Exact nodal biquintic potential",
          "legacy": "Legacy bicubic $C_V$ column",
          "native_1.0": "Smooth potential: native-derivative fit",
          "thermo_0.1": "Smooth potential: thermodynamic fit"}


def summary(values):
    values = np.asarray(values)
    if not values.size or not np.isfinite(values).all():
        raise ValueError("Empty or nonfinite diagnostic")
    return {"count": int(values.size), "signed_min": float(values.min()),
            "signed_max": float(values.max()), "abs_median": float(np.median(abs(values))),
            "abs_p99": float(np.quantile(abs(values), .99)), "abs_max": float(abs(values).max())}


def legacy_state(cache, logq, logt, sample_q, sample_t):
    values = legacy.evaluate_spline(cache, logq, logt, sample_q, sample_t)
    derivative_q = legacy.evaluate_spline(cache, logq, logt, sample_q, sample_t, (1, 0))
    derivative_t = legacy.evaluate_spline(cache, logq, logt, sample_q, sample_t, (0, 1))
    cv_column, entropy_column = legacy.FIELDS["Cv"], legacy.FIELDS["lnS"]
    entropy = np.exp(values[:, entropy_column])
    temperature_entropy_partial = entropy * (derivative_t[:, entropy_column] -
                                               2 * derivative_q[:, entropy_column]) / LN10
    cv = values[:, cv_column]
    return {"Cv": cv, "Cv_x": derivative_q[:, cv_column] / LN10,
            "Cv_t": (derivative_t[:, cv_column] - 2 * derivative_q[:, cv_column]) / LN10,
            "S": entropy, "T_s_T": temperature_entropy_partial,
            "dse_Cv": temperature_entropy_partial / cv - 1}


def potential_state(native, sample_q, sample_t):
    state = thermodynamic_partials(native, sample_q, sample_t)
    physical = to_physical(native, EVAL_ORDERS)
    temperature_entropy_partial = (physical[0, 1] - physical[0, 2]) / (10.0**sample_t)
    return {"Cv": state["Cv"], "Cv_x": state["Cv_x"], "Cv_t": state["Cv_t"],
            "S": state["S"], "T_s_T": temperature_entropy_partial,
            "dse_Cv": temperature_entropy_partial / state["Cv"] - 1}


def finite_differences(evaluators, logq, logt, sample_q, sample_t, faces):
    indices = np.arange(3, len(sample_q), 37)
    indices = indices[np.min(abs(sample_q[indices, None] - faces[None, :]), axis=1) > 2e-4]
    center_q, center_t = sample_q[indices], sample_t[indices]
    center_cell_q, center_cell_t, _ = query_cells(logq, logt, center_q, center_t)
    report = {}
    for name, evaluator in evaluators.items():
        center = evaluator(center_q, center_t)
        report[name] = {}
        for step in (1e-4, 3e-5, 1e-5):
            shift_q, shift_t = -2 * step / LN10, step / LN10
            for sign in (-1, 1):
                other_q, other_t, inside = query_cells(logq, logt, center_q + sign * shift_q, center_t + sign * shift_t)
                if not np.all(inside & (other_q == center_cell_q) & (other_t == center_cell_t)):
                    raise AssertionError("Entropy finite difference crosses a native face")
            plus = evaluator(center_q + shift_q, center_t + shift_t)["S"]
            minus = evaluator(center_q - shift_q, center_t - shift_t)["S"]
            numerical = (plus - minus) / (2 * step)
            scale = np.maximum.reduce([np.ones(len(indices)), abs(center["Cv"]), abs(center["T_s_T"])])
            error = (numerical - center["T_s_T"]) / scale
            report[name][str(step)] = summary(error)
            if report[name][str(step)]["abs_max"] > 1e-3:
                raise AssertionError(f"{name}: returned-entropy finite difference failed")
    return {"queries": int(len(indices)), "steps_natural_logT": [1e-4, 3e-5, 1e-5],
            "normalization": "max(1,abs(Cv),abs(T*s_T)); fixed density; within native cells",
            "scaled_limit": 1e-3, "results": report}


def plot(sample_q, sample_t, curves, face_q, traces):
    os.environ["MPLCONFIGDIR"] = str(OUTPUT / ".mplconfig")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(4, 1, figsize=(9.2, 9.1), sharex=True, constrained_layout=True)
    density = sample_q + 2 * sample_t - 12
    fields = ("Cv", "Cv_x", "Cv_t", "dse_Cv")
    cuts = np.r_[sample_q[0], face_q, sample_q[-1]]
    for name in NAMES:
        state = curves[name]
        for axis, field in zip(axes, fields):
            if name != "biquintic":
                values = np.maximum(abs(state[field]), FLOOR) if field == "dse_Cv" else state[field] / 1e8
                axis.plot(density, values, color=COLORS[name], linestyle="--" if name == "legacy" else "-",
                          linewidth=1.5, label=LABELS[name], zorder=4 if name == "legacy" else 2)
                continue
            for index in range(len(cuts) - 1):
                inside = (sample_q > cuts[index]) & (sample_q < cuts[index + 1])
                segment_q = np.r_[cuts[index], sample_q[inside], cuts[index + 1]]
                left = state[field][0] if index == 0 else traces["right"][field][index - 1]
                right = state[field][-1] if index == len(cuts) - 2 else traces["left"][field][index]
                values = np.r_[left, state[field][inside], right]
                displayed = np.maximum(abs(values), FLOOR) if field == "dse_Cv" else values / 1e8
                axis.plot(segment_q + 2 * sample_t[0] - 12, displayed, color=COLORS[name],
                          linewidth=1.5, label=LABELS[name] if index == 0 else None, zorder=3)
            if field != "Cv":
                for side, marker in (("left", "o"), ("right", "x")):
                    values = traces[side][field]
                    displayed = np.maximum(abs(values), FLOOR) if field == "dse_Cv" else values / 1e8
                    axis.plot(face_q + 2 * sample_t[0] - 12, displayed, linestyle="none", marker=marker,
                              markersize=3, color=COLORS[name], zorder=5)
    for coordinate in face_q:
        for axis in axes:
            axis.axvline(coordinate + 2 * sample_t[0] - 12, color="gray", alpha=.2, linewidth=.7)
    axes[0].set_ylabel("$C_V/10^8$\n[erg g$^{-1}$ K$^{-1}$]")
    axes[1].set_ylabel("$(\\partial C_V/\\partial\\ln\\rho)_T/10^8$")
    axes[2].set_ylabel("$(\\partial C_V/\\partial\\ln T)_\\rho/10^8$")
    axes[3].set_ylabel(r"$|\mathrm{dse}_{C_V}|=|T s_T/C_V-1|$")
    axes[3].set_yscale("log")
    axes[3].set_ylim(bottom=FLOOR / 3)
    axes[3].text(.01, .04, r"Exact zeros shown at $10^{-16}$ display floor; potential curves overlap",
                 transform=axes[3].transAxes, fontsize=8)
    axes[3].set_xlabel("$\\log_{10}\\rho$ [g cm$^{-3}$]")
    axes[0].legend(fontsize=8, ncol=2)
    figure.suptitle("FreeEOS-only diagnostic: audited default MESA selection is Skye-only\n"
                   "Original 036 fits and trajectory; blue circles/crosses = exact left/right face traces",
                   fontsize=10)
    for extension in ("pdf", "png"):
        figure.savefig(OUTPUT / f"figure2_consistency.{extension}", dpi=180)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot-only", action="store_true")
    arguments = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    if arguments.plot_only:
        result_path = OUTPUT / "results.json"
        result = json.loads(result_path.read_text())
        arrays_path = OUTPUT / "curves_and_traces.npz"
        if checksum(arrays_path) != result["output_sha256"]["curves_and_traces.npz"]:
            raise ValueError("Saved 068 arrays fail their checksum")
        saved_output = np.load(arrays_path)
        sample_q, sample_t, face_q = (saved_output[key] for key in
                                      ("sample_logQ", "sample_logT", "face_logQ"))
        fields = ("Cv", "Cv_x", "Cv_t", "S", "T_s_T", "dse_Cv")
        curves = {name: {field: saved_output[f"{name}_{field}"] for field in fields} for name in NAMES}
        traces = {side: {field: saved_output[f"{side}_{field}"] for field in fields}
                  for side in ("left", "right")}
        plot(sample_q, sample_t, curves, face_q, traces)
        result["provenance"]["script_sha256"] = checksum(__file__)
        result["output_sha256"] = {path.name: checksum(path) for path in
                                   (arrays_path, OUTPUT / "figure2_consistency.pdf",
                                    OUTPUT / "figure2_consistency.png")}
        result_path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        return
    archived_result_path = ARCHIVE / "results.json"
    archived_arrays_path = ARCHIVE / "curves_and_traces.npz"
    archived_result = json.loads(archived_result_path.read_text())
    for name, expected in archived_result["artifact_sha256"].items():
        if checksum(ARCHIVE / name) != expected:
            raise ValueError(f"036 archive hash mismatch: {name}")
    protected = [archived_result_path, *(ARCHIVE / name for name in archived_result["artifact_sha256"]),
                 CODE / "freeeos_potential.py",
                 CODE / "smooth_patch.py", WORKSPACE / "freeeos_report/audit_freeeos.py"]
    protected.extend(WORKSPACE / path for path in archived_result["input_sha256"])
    before = {str(path.relative_to(WORKSPACE)): checksum(path) for path in protected}
    for name, expected in archived_result["input_sha256"].items():
        if before[name] != expected:
            raise ValueError(f"036 input hash mismatch: {name}")
    saved = np.load(archived_arrays_path)
    sample_q, sample_t = saved["sample_logQ"], saved["sample_logT"]
    face_q = saved["face_logQ"]
    if len(sample_q) != 2001 or not np.array_equal(sample_q, np.linspace(2.25, 2.75, 2001)) or not np.all(sample_t == 6.5117106292763225):
        raise ValueError("036 trajectory mismatch")
    table = WORKSPACE / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
    cache_path = WORKSPACE / "mesa/data/eosDT_data/cache/mesa-FreeEOS_02z70x.bin"
    _, logq, logt, rows = read_table(table)
    low, valid, _ = source_constraints(rows, logq, logt)
    nodes, complete, _ = construct_twists(low, valid, np.diff(logq)[0] * LN10, np.diff(logt)[0] * LN10)
    _, cache = legacy.read_cache(cache_path, len(logq), len(logt))
    patch_q, patch_t = saved["patch_logQ"], saved["patch_logT"]
    evaluators = {"biquintic": lambda q, t: potential_state(evaluate(nodes, complete, logq, logt, q, t), q, t),
                  "legacy": lambda q, t: legacy_state(cache, logq, logt, q, t)}
    for name in ("native_1.0", "thermo_0.1"):
        coefficients = saved[f"{name}_coefficients"]
        evaluators[name] = lambda q, t, coefficients=coefficients: potential_state(
            evaluate_fit(coefficients, patch_q * LN10, patch_t * LN10, q * LN10, t * LN10), q, t)
    curves = {name: evaluator(sample_q, sample_t) for name, evaluator in evaluators.items()}
    face_indices = np.searchsorted(logq, face_q)
    if not np.array_equal(logq[face_indices], face_q):
        raise ValueError("036 exact faces mismatch")
    face_t = np.full_like(face_q, sample_t[0])
    temperature_indices = np.full(len(face_q), np.searchsorted(logt, sample_t[0]) - 1)
    traces = {side: potential_state(evaluate(nodes, complete, logq, logt, face_q, face_t,
              cell_indices=(face_indices + offset, temperature_indices)), face_q, face_t)
              for side, offset in (("left", -1), ("right", 0))}
    reproduction = {}
    for name in NAMES:
        reproduction[name] = {}
        for field in ("Cv", "Cv_x", "Cv_t"):
            difference = curves[name][field] - saved[f"{name}_{field}"]
            scale = np.maximum(1, abs(saved[f"{name}_{field}"]))
            reproduction[name][field] = summary(difference / scale)
            if reproduction[name][field]["abs_max"] > 1e-11:
                raise AssertionError(f"036 curve reproduction failed: {name} {field}")
    reproduction["exact_faces"] = {}
    for side in ("left", "right"):
        reproduction["exact_faces"][side] = {}
        for field in ("Cv", "Cv_x", "Cv_t"):
            difference = traces[side][field] - saved[f"{side}_{field}"]
            scale = np.maximum(1, abs(saved[f"{side}_{field}"]))
            reproduction["exact_faces"][side][field] = summary(difference / scale)
            if reproduction["exact_faces"][side][field]["abs_max"] > 1e-11:
                raise AssertionError(f"036 face reproduction failed: {side} {field}")
    finite_difference = finite_differences(evaluators, logq, logt, sample_q, sample_t, face_q)
    payload = {"sample_logQ": sample_q, "sample_logT": sample_t, "face_logQ": face_q}
    payload.update({f"{name}_{field}": values for name, state in curves.items() for field, values in state.items()})
    payload.update({f"{side}_{field}": values for side, state in traces.items() for field, values in state.items()})
    np.savez_compressed(OUTPUT / "curves_and_traces.npz", **payload)
    plot(sample_q, sample_t, curves, face_q, traces)
    after = {str(path.relative_to(WORKSPACE)): checksum(path) for path in protected}
    if before != after:
        raise AssertionError("Read-only input changed")
    result = {"task": "068 Figure 2 returned-entropy consistency", "archive_reproduction_relative": reproduction,
              "archive_guard": 1e-11, "finite_difference": finite_difference,
              "metric": "signed T*(partial s/partial T)_rho/Cv-1, returned entropy and independently returned Cv",
              "potential_convention": "036 source constraints with default radiation split and 036 saved fits",
              "legacy_chain_rule": "exp(lnS)*(d_lnS/d_logT-2*d_lnS/d_logQ)/ln(10)",
              "display_floor": FLOOR, "sampling": archived_result["sampling"],
              "dse_Cv": {name: {**summary(curves[name]["dse_Cv"]),
                        "exact_zeros": int(np.count_nonzero(curves[name]["dse_Cv"] == 0)),
                        "below_display_floor": int(np.count_nonzero(abs(curves[name]["dse_Cv"]) < FLOOR))}
                        for name in NAMES},
              "face_dse_Cv": {side: summary(traces[side]["dse_Cv"]) for side in traces},
              "provenance": {"input_sha256": before, "inputs_unchanged": True, "new_fit_calls": 0,
                             "script_sha256": checksum(__file__), "numpy": np.__version__},
              "output_sha256": {path.name: checksum(path) for path in
                                (OUTPUT / "curves_and_traces.npz", OUTPUT / "figure2_consistency.pdf",
                                 OUTPUT / "figure2_consistency.png")}}
    (OUTPUT / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"archive_reproduction": reproduction, "dse_Cv": result["dse_Cv"],
                      "finite_difference": finite_difference}, indent=2))


if __name__ == "__main__":
    main()
