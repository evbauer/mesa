"""Four-panel FreeEOS-only Figure 3 sidecar from authenticated 023 inputs."""

import argparse
import json
import os
from pathlib import Path
import shlex
import sys

import numpy as np

sys.dont_write_bytecode = True
CODE = Path(__file__).resolve().parents[1]
WORKSPACE = CODE.parents[3]
ARCHIVE = WORKSPACE / "freeeos_biquintic_progress/experiments/023_active_slice"
OUTPUT = WORKSPACE / "freeeos_biquintic_progress/experiments/066_extended_figure3"
sys.path.insert(0, str(CODE))
sys.path.insert(0, str(WORKSPACE / "freeeos_report"))
import audit_freeeos as baseline
from freeeos_potential import (EVAL_ORDERS, LN10, checksum, construct_twists, evaluate,
                               query_cells, read_table, source_constraints,
                               thermodynamic_partials, to_physical)

DISPLAY_FLOOR = 1e-16
FIELDS = ("Cv", "Cv_x", "Cv_t", "S", "T_s_T", "dse_Cv", "E", "E_T", "eT_Cv")
STEPS = (1e-4, 3e-5, 1e-5)
FD_GUARD = 1e-3
ARCHIVE_GUARD = 1e-8


def legacy_state(cache, logq, logt, sample_q, sample_t):
    values = baseline.evaluate_spline(cache, logq, logt, sample_q, sample_t)
    derivative_q = baseline.evaluate_spline(cache, logq, logt, sample_q, sample_t, (1, 0))
    derivative_t = baseline.evaluate_spline(cache, logq, logt, sample_q, sample_t, (0, 1))
    cv_field, entropy_field = baseline.FIELDS["Cv"], baseline.FIELDS["lnS"]
    entropy = np.exp(values[:, entropy_field])
    temperature_entropy_partial = entropy * (derivative_t[:, entropy_field] -
                                               2 * derivative_q[:, entropy_field]) / LN10
    cv = values[:, cv_field]
    energy_field = baseline.FIELDS["lnE"]
    energy = np.exp(values[:, energy_field])
    energy_partial = energy / (10.0**sample_t) * (derivative_t[:, energy_field] -
                                                  2 * derivative_q[:, energy_field]) / LN10
    return {"Cv": cv, "Cv_x": derivative_q[:, cv_field] / LN10,
            "Cv_t": (derivative_t[:, cv_field] - 2 * derivative_q[:, cv_field]) / LN10,
            "S": entropy, "T_s_T": temperature_entropy_partial,
            "dse_Cv": temperature_entropy_partial / cv - 1,
            "E": energy, "E_T": energy_partial, "eT_Cv": energy_partial / cv - 1}


def potential_state(nodes, complete, logq, logt, sample_q, sample_t, cell_indices=None):
    native = evaluate(nodes, complete, logq, logt, sample_q, sample_t, cell_indices=cell_indices)
    state = thermodynamic_partials(native, sample_q, sample_t)
    physical = to_physical(native, EVAL_ORDERS)
    temperature_entropy_partial = (physical[0, 1] - physical[0, 2]) / (10.0**sample_t)
    energy_partial = (physical[0, 1] - physical[0, 2]) / (10.0**sample_t)
    return {"Cv": state["Cv"], "Cv_x": state["Cv_x"], "Cv_t": state["Cv_t"],
            "S": state["S"], "T_s_T": temperature_entropy_partial,
            "dse_Cv": temperature_entropy_partial / state["Cv"] - 1,
            "E": state["E"], "E_T": energy_partial, "eT_Cv": energy_partial / state["Cv"] - 1}


def summary(values):
    values = np.asarray(values)
    if not np.all(np.isfinite(values)) or not values.size:
        raise ValueError("Nonfinite or empty diagnostic")
    return {"count": int(values.size), "abs_median": float(np.median(abs(values))),
            "abs_p99": float(np.quantile(abs(values), .99)), "abs_max": float(np.max(abs(values))),
            "signed_min": float(np.min(values)), "signed_max": float(np.max(values))}


def finite_differences(evaluators, states, logq, logt, sample_q, sample_t):
    margin_q = np.min(abs(sample_q[:, None] - logq), axis=1)
    margin_t = np.min(abs(sample_t[:, None] - logt), axis=1)
    safe = (margin_q > 2 * max(STEPS) / LN10 + 1e-12) & (margin_t > max(STEPS) / LN10 + 1e-12)
    center_q, center_t = sample_q[safe], sample_t[safe]
    cell_q, cell_t, _ = query_cells(logq, logt, center_q, center_t)
    report, arrays = {}, {"fd_safe_mask": safe}
    for name, evaluator in evaluators.items():
        report[name] = {}
        for step in STEPS:
            report[name][str(step)] = {}
            for direction in ("x", "t"):
                shift_q = step / LN10 if direction == "x" else -2 * step / LN10
                shift_t = 0 if direction == "x" else step / LN10
                for sign in (-1, 1):
                    other_q, other_t, inside = query_cells(logq, logt, center_q + sign * shift_q, center_t + sign * shift_t)
                    if not np.all(inside & (other_q == cell_q) & (other_t == cell_t)):
                        raise ValueError("Finite difference crosses a native cell face")
                plus = evaluator(center_q + shift_q, center_t + shift_t)
                minus = evaluator(center_q - shift_q, center_t - shift_t)
                field = "Cv_" + direction
                numerical = (plus["Cv"] - minus["Cv"]) / (2 * step)
                scale = np.maximum.reduce([np.ones(safe.sum()), abs(states[name]["Cv"][safe]), abs(states[name][field][safe])])
                error = (numerical - states[name][field][safe]) / scale
                report[name][str(step)][field] = summary(error)
                arrays[f"fd__{name}__{step}__{field}_scaled_error"] = error
                if direction == "t":
                    numerical_entropy = (plus["S"] - minus["S"]) / (2 * step)
                    entropy_scale = np.maximum.reduce([np.ones(safe.sum()), abs(states[name]["Cv"][safe]),
                                                       abs(states[name]["T_s_T"][safe])])
                    error_entropy = (numerical_entropy - states[name]["T_s_T"][safe]) / entropy_scale
                    report[name][str(step)]["T_s_T"] = summary(error_entropy)
                    arrays[f"fd__{name}__{step}__T_s_T_scaled_error"] = error_entropy
                    numerical_energy = (plus["E"] - minus["E"]) / (2 * step * 10.0**center_t)
                    energy_scale = np.maximum.reduce([np.ones(safe.sum()), abs(states[name]["Cv"][safe]),
                                                      abs(states[name]["E_T"][safe])])
                    error_energy = (numerical_energy - states[name]["E_T"][safe]) / energy_scale
                    report[name][str(step)]["E_T"] = summary(error_energy)
                    arrays[f"fd__{name}__{step}__E_T_scaled_error"] = error_energy
    passed = all(data["abs_max"] <= FD_GUARD for model in report.values()
                 for step in model.values() for data in step.values())
    return {"safe_count": int(safe.sum()), "excluded_count": int((~safe).sum()),
            "scaled_limit": FD_GUARD, "passed": passed, "steps": report}, arrays


def plot(output):
    os.environ["MPLCONFIGDIR"] = str(output / ".mplconfig")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    data = np.load(output / "curves.npz")
    result = json.loads((output / "results.json").read_text())
    figure, axes = plt.subplots(4, 1, figsize=(11.5, 10), sharex=True, constrained_layout=True)
    shading = ("#b7dfba", "#ffe0a3", "#d6d6d6")
    for axis in axes:
        for span in result["measured_region_spans"]:
            axis.axvspan(*span["shading_logQ_bounds_midpoint_estimate"], color=shading[span["code"]], alpha=.28, zorder=0)
        for knot in data["native_logQ_guides"]:
            axis.axvline(knot, color=".6", linewidth=.5, alpha=.45, zorder=0)
        for knot in (.77, .80):
            axis.axvline(knot, color=".25", linestyle=":", linewidth=.8, zorder=0)
        axis.grid(axis="y", alpha=.2)
    labels = {"BQ": "Biquintic free-energy reconstruction", "legacy": "Legacy interpolated functions"}
    for name, color, style in (("BQ", "tab:blue", "-"), ("legacy", "tab:orange", "--")):
        for segment in np.unique(data["trace_segment"]):
            selected = data["trace_segment"] == segment
            for axis, field in zip(axes, ("Cv", "Cv_x", "Cv_t", "dse_Cv")):
                values = data[f"trace__{name}__{field}"][selected]
                displayed = np.maximum(abs(values), DISPLAY_FLOOR) if field == "dse_Cv" else values / 1e8
                axis.plot(data["trace_logQ"][selected], displayed, color=color, linestyle=style, linewidth=1.4,
                          label=labels[name] if segment == 0 else None)
                if field in ("Cv_x", "Cv_t"):
                    axis.plot(data["trace_logQ"][selected][[0, -1]], displayed[[0, -1]], linestyle="none",
                              marker="o", markersize=2.5, markerfacecolor="white", markeredgecolor=color, markeredgewidth=.65)
    axes[0].legend(loc="best", fontsize=9, ncol=2)
    axes[0].set_title(f"FreeEOS-only slice: log T = {result['logT']:.15f}; X=.70, Z=.02\n"
                      "No thermodynamic blending; shading reuses measured default-MESA routing", fontsize=11)
    axes[0].set_ylabel(r"$C_V$ [$10^8$ erg g$^{-1}$ K$^{-1}$]")
    axes[1].set_ylabel(r"$(\partial C_V/\partial\ln\rho)_T$" + "\n" + r"[$10^8$ erg g$^{-1}$ K$^{-1}$]")
    axes[2].set_ylabel(r"$(\partial C_V/\partial\ln T)_\rho$" + "\n" + r"[$10^8$ erg g$^{-1}$ K$^{-1}$]")
    axes[3].set_ylabel(r"$|\mathrm{dse}_{C_V}|$" + "\n" + r"$=|T s_T/C_V-1|$")
    axes[3].set_yscale("log")
    axes[3].set_ylim(bottom=DISPLAY_FLOOR / 3)
    axes[3].text(.015, .06, r"Blue exact zeros shown at $10^{-16}$ display floor", transform=axes[3].transAxes, fontsize=8)
    axes[3].legend(handles=[Patch(facecolor=shading[span["code"]], alpha=.4, label=span["label"])
                            for span in result["measured_region_spans"]], loc="center right", fontsize=8)
    axes[3].set_xlabel(r"$\log_{10}Q=\log_{10}\rho-2\log_{10}T+12$; guides: native Q knots; circles: cell traces")
    axes[3].set_xlim(float(data["logQ"].min()), float(data["logQ"].max()))
    for extension in ("pdf", "png"):
        figure.savefig(output / f"extended_freeeos_slice.{extension}", dpi=180)
    plt.close(figure)
    comparison, comparison_axes = plt.subplots(2, 1, figsize=(10, 5.5), sharex=True, sharey=True, constrained_layout=True)
    for axis, field, label in zip(comparison_axes, ("eT_Cv", "dse_Cv"),
                                 (r"$|e_T/C_V-1|$", r"$|\mathrm{dse}_{C_V}|=|T s_T/C_V-1|$")):
        for span in result["measured_region_spans"]:
            axis.axvspan(*span["shading_logQ_bounds_midpoint_estimate"], color=shading[span["code"]], alpha=.28, zorder=0)
        for knot in data["native_logQ_guides"]:
            axis.axvline(knot, color=".6", linewidth=.5, alpha=.45, zorder=0)
        for name, color, style in (("BQ", "tab:blue", "-"), ("legacy", "tab:orange", "--")):
            for segment in np.unique(data["trace_segment"]):
                selected = data["trace_segment"] == segment
                axis.plot(data["trace_logQ"][selected], np.maximum(abs(data[f"trace__{name}__{field}"][selected]), DISPLAY_FLOOR),
                          color=color, linestyle=style, linewidth=1.4, label=labels[name] if segment == 0 else None)
        axis.set_yscale("log")
        axis.set_ylabel(label)
        axis.grid(axis="y", alpha=.2)
        axis.set_ylim(bottom=DISPLAY_FLOOR / 3)
    comparison_axes[0].legend(loc="center right", fontsize=8)
    comparison_axes[0].set_title("FreeEOS-only returned-function consistency on the same 023 slice")
    comparison_axes[1].set_xlabel(r"$\log_{10}Q$; blue exact zeros use the $10^{-16}$ display floor")
    for extension in ("pdf", "png"):
        comparison.savefig(output / f"consistency_comparison.{extension}", dpi=180)
    plt.close(comparison)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot-only", action="store_true")
    arguments = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    if arguments.plot_only:
        replay = json.loads((OUTPUT / "results.json").read_text())
        if checksum(OUTPUT / "curves.npz") != replay["output_sha256"]["curves.npz"]:
            raise ValueError("Saved066 plot arrays fail their checksum")
        plot(OUTPUT)
        replay["output_sha256"] = {name: checksum(OUTPUT / name) for name in replay["output_sha256"]}
        replay["plot_script_sha256"] = checksum(__file__)
        (OUTPUT / "results.json").write_text(json.dumps(replay, indent=2, allow_nan=False) + "\n")
        return
    source = WORKSPACE / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
    cache_path = WORKSPACE / "mesa/data/eosDT_data/cache/mesa-FreeEOS_02z70x.bin"
    metadata_path, archive_path = ARCHIVE / "results.json", ARCHIVE / "queries.npz"
    metadata = json.loads(metadata_path.read_text())
    protected = [source, cache_path, metadata_path, archive_path, ARCHIVE / "probe.tsv", ARCHIVE / "probe_metadata.json",
                 CODE / "freeeos_potential.py", WORKSPACE / "freeeos_report/audit_freeeos.py", Path(__file__),
                 WORKSPACE / "mesa/eos/private/eosdt_load_tables.f90", WORKSPACE / "mesa/eos/private/eosdt_eval.f90"]
    before = {str(path.relative_to(WORKSPACE)): checksum(path) for path in protected}
    if checksum(source) != metadata["source_sha256"] or checksum(cache_path) != metadata["legacy_cache_sha256"]:
        raise ValueError("Source/cache differ from original023")
    if checksum(archive_path) != metadata["output_artifacts_sha256"]["queries.npz"]:
        raise ValueError("023 query archive checksum mismatch")
    if checksum(ARCHIVE / "probe.tsv") != metadata["output_artifacts_sha256"]["probe.tsv"]:
        raise ValueError("023 measured probe checksum mismatch")
    saved = np.load(archive_path)
    sample_q, sample_t = saved["logQ"], saved["logT"]
    if len(sample_q) != 1001 or not np.all(sample_t == metadata["logT"]):
        raise ValueError("Unexpected 023 trajectory")
    _, logq, logt, rows = read_table(source)
    low, valid, _ = source_constraints(rows, logq, logt, split_radiation=False, minimum_logT=3.2)
    nodes, complete, _ = construct_twists(low, valid, np.diff(logq)[0] * LN10, np.diff(logt)[0] * LN10, radius=2, route="mean")
    _, cache = baseline.read_cache(cache_path, len(logq), len(logt))
    evaluators = {"BQ": lambda density, temperature: potential_state(nodes, complete, logq, logt, density, temperature),
                  "legacy": lambda density, temperature: legacy_state(cache, logq, logt, density, temperature)}
    states = {name: evaluator(sample_q, sample_t) for name, evaluator in evaluators.items()}
    archive_check = {}
    for name, state in states.items():
        old_name = "BQ_exact" if name == "BQ" else name
        if np.any(state["Cv"] <= 0) or any(not np.all(np.isfinite(values)) for values in state.values()):
            raise ValueError("Invalid Cv or nonfinite state")
        archive_check[name] = {}
        for field in ("Cv", "Cv_t"):
            scale = np.maximum.reduce([np.ones(len(sample_q)), abs(saved[f"{old_name}__Cv"]), abs(saved[f"{old_name}__{field}"])])
            archive_check[name][field] = summary((state[field] - saved[f"{old_name}__{field}"]) / scale)
            if archive_check[name][field]["abs_max"] > ARCHIVE_GUARD:
                raise ValueError("Original023 state reproduction guard failed")
    density_index, temperature_index, _ = query_cells(logq, logt, sample_q, sample_t)
    convention = {}
    for label, source_field, cache_field in (("lnS", 3, baseline.FIELDS["lnS"]), ("lnE", 2, baseline.FIELDS["lnE"]),
                                            ("Cv", 7, baseline.FIELDS["Cv"])):
        differences = []
        for density_end in (0, 1):
            for temperature_end in (0, 1):
                source_value = rows[density_index + density_end, temperature_index + temperature_end, source_field]
                if label in ("lnS", "lnE"):
                    source_value = source_value * LN10
                cache_value = cache[0, cache_field, density_index + density_end, temperature_index + temperature_end]
                differences.append((cache_value - source_value) / np.maximum(1, abs(source_value)))
        convention[label] = summary(np.concatenate(differences))
        if convention[label]["abs_max"] > 1e-12:
            raise ValueError("Source/cache nodal convention guard failed")
    fd, arrays = finite_differences(evaluators, states, logq, logt, sample_q, sample_t)
    arrays.update({key: saved[key] for key in saved.files if key in ("logQ", "logT", "logRho", "region_code", "native_logQ_guides") or key.startswith("MESA_frac_")})
    arrays.update({f"{name}__{field}": values for name, state in states.items() for field, values in state.items()})
    segments_q, segments_id = [], []
    segments = {name: {field: [] for field in FIELDS} for name in states}
    segment_count = 0
    for cell_index in range(len(logq) - 1):
        lower, upper = max(sample_q.min(), logq[cell_index]), min(sample_q.max(), logq[cell_index + 1])
        if upper - lower <= 1e-12:
            continue
        internal = sample_q[(sample_q > lower + 1e-12) & (sample_q < upper - 1e-12)]
        trace_q = np.r_[lower, internal, upper]
        trace_t = np.full(len(trace_q), sample_t[0])
        forced = np.full(len(trace_q), cell_index), np.full(len(trace_q), temperature_index[0])
        traces = {"BQ": potential_state(nodes, complete, logq, logt, trace_q, trace_t, forced),
                  "legacy": legacy_state(cache, logq, logt, trace_q, trace_t)}
        segments_q.append(trace_q)
        segments_id.append(np.full(len(trace_q), segment_count))
        for name, state in traces.items():
            for field in FIELDS:
                segments[name][field].append(state[field])
        segment_count += 1
    arrays.update(trace_logQ=np.concatenate(segments_q), trace_segment=np.concatenate(segments_id))
    arrays.update({f"trace__{name}__{field}": np.concatenate(values) for name, fields in segments.items() for field, values in fields.items()})
    np.savez_compressed(OUTPUT / "curves.npz", **arrays)
    result = {"task": "066 four-panel Figure3 sidecar", "command": shlex.join([sys.executable, *sys.argv]),
              "logT": float(sample_t[0]), "query_count": len(sample_q), "cell_trace_segments": segment_count,
              "metric": "signed dse_Cv=T*(partial s/partial T)_rho/Cv-1; returned lnS derivative and returned Cv",
              "legacy_chain_rule": "S=exp(lnS); T*s_T=S*(d_lnS/d_logT-2*d_lnS/d_logQ)/ln(10)",
              "legacy_energy_chain_rule": "E=exp(lnE); E_T=E/T*(d_lnE/d_logT-2*d_lnE/d_logQ)/ln(10)",
              "potential_convention": "original023 total-potential, source floor3.2, radius2/mean; no analytic radiation split",
              "display_floor": DISPLAY_FLOOR, "archive_guard": ARCHIVE_GUARD, "archive_reproduction": archive_check,
              "source_cache_nodal_conventions": convention, "finite_difference": fd,
              "measured_region_spans": metadata["measured_region_spans"],
              "dse_Cv": {name: {**summary(state["dse_Cv"]), "exact_zeros": int(np.count_nonzero(state["dse_Cv"] == 0)),
                                "below_display_floor": int(np.count_nonzero(abs(state["dse_Cv"]) < DISPLAY_FLOOR))} for name, state in states.items()},
              "eT_Cv": {name: {**summary(state["eT_Cv"]), "exact_zeros": int(np.count_nonzero(state["eT_Cv"] == 0)),
                               "below_display_floor": int(np.count_nonzero(abs(state["eT_Cv"]) < DISPLAY_FLOOR))} for name, state in states.items()},
              "trace_metric_display": {name: {field: {"count": int(len(arrays[f"trace__{name}__{field}"])),
                                                       "exact_zeros": int(np.count_nonzero(arrays[f"trace__{name}__{field}"] == 0)),
                                                       "below_display_floor": int(np.count_nonzero(abs(arrays[f"trace__{name}__{field}"]) < DISPLAY_FLOOR))}
                                              for field in ("dse_Cv", "eT_Cv")} for name in states},
              "legacy_metric_comparison": {
                  "signed_correlation": float(np.corrcoef(states["legacy"]["eT_Cv"], states["legacy"]["dse_Cv"])[0, 1]),
                  "regions": {span["label"]: {field: summary(states["legacy"][field][saved["region_code"] == span["code"]])
                                              for field in ("eT_Cv", "dse_Cv")} for span in metadata["measured_region_spans"]},
                  "near_feature_0.74_to_0.82": {field: summary(states["legacy"][field][(sample_q >= .74) & (sample_q <= .82)]) for field in ("eT_Cv", "dse_Cv")},
                  "marked_faces": {str(face): {field: float(states["legacy"][field][np.argmin(abs(sample_q - face))])
                                              for field in ("eT_Cv", "dse_Cv")} for face in (.77, .80)},
                  "maxima": {field: {"logQ": float(sample_q[np.argmax(abs(states["legacy"][field]))]),
                                      "signed_value": float(states["legacy"][field][np.argmax(abs(states["legacy"][field]))])}
                             for field in ("eT_Cv", "dse_Cv")}},
              "provenance": {"input_sha256": before, "original023_array_sha256": checksum(archive_path),
                             "new_MESA_probe_calls": 0, "python": sys.version, "numpy": np.__version__}}
    after = {str(path.relative_to(WORKSPACE)): checksum(path) for path in protected}
    result["provenance"]["inputs_unchanged"] = before == after
    (OUTPUT / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    if before != after or not fd["passed"]:
        raise ValueError("Input preservation or within-cell finite-difference guard failed; diagnostics saved")
    plot(OUTPUT)
    result["output_sha256"] = {name: checksum(OUTPUT / name) for name in ("curves.npz", "extended_freeeos_slice.pdf", "extended_freeeos_slice.png",
                                                                        "consistency_comparison.pdf", "consistency_comparison.png")}
    (OUTPUT / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"archive_reproduction": archive_check, "finite_difference": fd, "dse_Cv": result["dse_Cv"]}, indent=2))


if __name__ == "__main__":
    main()
