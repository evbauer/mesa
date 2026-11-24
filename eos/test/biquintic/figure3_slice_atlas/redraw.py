"""Authenticated, plot-only density-axis replay of experiment 066."""

import hashlib
import json
import os
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[5]
SOURCE = ROOT / "freeeos_biquintic_progress/experiments/066_extended_figure3"
OUTPUT = ROOT / "freeeos_biquintic_progress/experiments/071_figure3_slice_atlas"
FLOOR = 1e-16


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    archived = json.loads((SOURCE / "results.json").read_text())
    source_arrays = SOURCE / "curves.npz"
    if digest(source_arrays) != archived["output_sha256"]["curves.npz"]:
        raise ValueError("066 saved arrays differ from archived checksum")
    with np.load(source_arrays) as data:
        arrays = {key: data[key] for key in data.files}
    logt = float(archived["logT"])
    offset = 2 * logt - 12
    if not np.allclose(arrays["logRho"], arrays["logQ"] + offset, rtol=0, atol=3e-14):
        raise ValueError("Saved decimal density does not match exact fixed-T mapping")
    if not np.all(arrays["logT"] == logt):
        raise ValueError("Unexpected temperature variation")
    for name in ("BQ", "legacy"):
        residual = arrays[f"{name}__T_s_T"] / arrays[f"{name}__Cv"] - 1
        if not np.array_equal(residual, arrays[f"{name}__dse_Cv"]):
            raise ValueError(f"Saved {name} entropy/Cv metric differs from returned functions")
    os.environ["MPLCONFIGDIR"] = str(OUTPUT / ".mplconfig")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    styles = (("BQ", "tab:blue", "-", "Total-potential biquintic"),
              ("legacy", "tab:orange", "--", "Legacy interpolated functions"))
    guides = arrays["native_logQ_guides"] + offset
    segments = np.unique(arrays["trace_segment"])
    density = arrays["trace_logQ"] + offset
    figure, axes = plt.subplots(4, 1, figsize=(11.5, 10), sharex=True, constrained_layout=True)
    for axis in axes:
        for guide in guides:
            axis.axvline(guide, color=".65", linewidth=.5, alpha=.45, zorder=0)
        for guide in (.77, .80):
            axis.axvline(guide + offset, color=".3", linewidth=.8, linestyle=":", zorder=0)
        axis.grid(axis="y", alpha=.2)
    for name, color, linestyle, label in styles:
        for segment in segments:
            selected = arrays["trace_segment"] == segment
            for axis, field in zip(axes, ("Cv", "Cv_x", "Cv_t", "dse_Cv")):
                values = arrays[f"trace__{name}__{field}"][selected]
                shown = np.maximum(abs(values), FLOOR) if field == "dse_Cv" else values / 1e8
                axis.plot(density[selected], shown, color=color, linestyle=linestyle, linewidth=1.4,
                          label=label if segment == 0 else None)
                if field in ("Cv_x", "Cv_t"):
                    axis.plot(density[selected][[0, -1]], shown[[0, -1]], linestyle="none",
                              marker="o", markersize=2.5, markerfacecolor="white",
                              markeredgecolor=color, markeredgewidth=.65)
    axes[0].legend(loc="best", fontsize=9, ncol=2)
    axes[0].set_title(f"FreeEOS-only slice: $\\log_{{10}}T={logt:.15f}$; $X=0.70$, $Z=0.02$\n"
                      "Unblended returned functions; dotted lines: selected native density faces", fontsize=11)
    axes[0].set_ylabel(r"$C_V$ [$10^8$ erg g$^{-1}$ K$^{-1}$]")
    axes[1].set_ylabel(r"$(\partial C_V/\partial\ln\rho)_T$" + "\n" + r"[$10^8$ erg g$^{-1}$ K$^{-1}$]")
    axes[2].set_ylabel(r"$(\partial C_V/\partial\ln T)_\rho$" + "\n" + r"[$10^8$ erg g$^{-1}$ K$^{-1}$]")
    axes[3].set_ylabel(r"$|\mathrm{dse}_{C_V}|=|T s_T/C_V-1|$")
    axes[3].set_yscale("log")
    axes[3].set_ylim(bottom=FLOOR / 3)
    axes[3].text(.015, .06, r"Exact blue arithmetic zeros displayed at $10^{-16}$", transform=axes[3].transAxes,
                 fontsize=8)
    axes[3].set_xlabel(r"$\log_{10}\rho$ [g cm$^{-3}$]; guides: native $Q$ faces; circles: exact cell endpoints")
    axes[3].set_xlim(float(arrays["logQ"].min() + offset), float(arrays["logQ"].max() + offset))
    secondary = axes[0].secondary_xaxis("top", functions=(lambda rho: rho - offset, lambda logq: logq + offset))
    secondary.set_xlabel(r"$\log_{10}Q=\log_{10}\rho-2\log_{10}T+12$")
    for extension in ("pdf", "png"):
        figure.savefig(OUTPUT / f"figure3_density_redraw.{extension}", dpi=180)
    plt.close(figure)

    signed_figure, signed_axes = plt.subplots(2, 1, figsize=(10.5, 6), constrained_layout=True)
    for axis, bounds in zip(signed_axes, ((.50, 1.10), (.74, .82))):
        for guide in guides:
            if bounds[0] <= guide - offset <= bounds[1]:
                axis.axvline(guide, color=".65", linewidth=.5, alpha=.45)
        axis.axhline(0, color=".3", linewidth=.7)
        for name, color, linestyle, label in styles:
            for segment in segments:
                selected = arrays["trace_segment"] == segment
                axis.plot(density[selected], arrays[f"trace__{name}__dse_Cv"][selected],
                          color=color, linestyle=linestyle, linewidth=1.1,
                          label=label if segment == 0 else None)
        axis.set_xlim(bounds[0] + offset, bounds[1] + offset)
        axis.set_ylabel(r"Signed $T s_T/C_V-1$")
        axis.grid(axis="y", alpha=.2)
    signed_axes[0].legend(loc="best", fontsize=8)
    signed_axes[0].set_title("Signed returned-entropy/Cv residual: whole slice and selected-face detail")
    signed_axes[1].set_xlabel(r"$\log_{10}\rho$ [g cm$^{-3}$]")
    for extension in ("pdf", "png"):
        signed_figure.savefig(OUTPUT / f"figure3_signed_residual.{extension}", dpi=180)
    plt.close(signed_figure)

    roots = []
    trace_q = arrays["trace_logQ"]
    trace_value = arrays["trace__legacy__dse_Cv"]
    for segment in segments:
        indices = np.flatnonzero(arrays["trace_segment"] == segment)
        for first, second in zip(indices[:-1], indices[1:]):
            low, high = trace_value[first], trace_value[second]
            if low * high < 0:
                estimate = (trace_q[first] * abs(high) + trace_q[second] * abs(low)) / (abs(low) + abs(high))
                roots.append({"bracket_logQ": [float(trace_q[first]), float(trace_q[second])],
                              "signed_endpoints": [float(low), float(high)],
                              "linear_estimate_logQ": float(estimate), "linear_estimate_logRho": float(estimate + offset),
                              "nearest_guide_distance_logQ": float(np.min(abs(arrays["native_logQ_guides"] - estimate)))})
    minima = []
    sample_q = arrays["logQ"]
    sample_value = arrays["legacy__dse_Cv"]
    for index in range(1, len(sample_q) - 1):
        if abs(sample_value[index]) < abs(sample_value[index - 1]) and abs(sample_value[index]) <= abs(sample_value[index + 1]):
            minima.append({"logQ": float(sample_q[index]), "logRho": float(sample_q[index] + offset),
                           "signed_value": float(sample_value[index]),
                           "neighbor_sign_change": bool(sample_value[index - 1] * sample_value[index + 1] < 0),
                           "nearest_guide_distance_logQ": float(np.min(abs(arrays["native_logQ_guides"] - sample_q[index])))})
    report = {"source_curves_sha256": digest(source_arrays), "source_results_sha256": digest(SOURCE / "results.json"),
              "logT": logt, "logRho_minus_logQ": offset, "query_count": len(sample_q),
              "cell_segments": len(segments), "native_guides": len(guides),
              "signed_metric_identity_checked": True, "zero_brackets_legacy": roots,
              "sampled_local_abs_minima_legacy": minima,
              "exact_BQ_arithmetic_zeros": int(np.count_nonzero(arrays["BQ__dse_Cv"] == 0)),
              "redraw_sha256": {name: digest(OUTPUT / name) for name in
                                ("figure3_density_redraw.pdf", "figure3_density_redraw.png",
                                 "figure3_signed_residual.pdf", "figure3_signed_residual.png")}}
    (OUTPUT / "redraw_results.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"zero_brackets": len(roots), "sampled_local_minima": len(minima),
                      "largest_minimum_guide_offset": max((row["nearest_guide_distance_logQ"] for row in minima), default=0),
                      "files": list(report["redraw_sha256"])}, indent=2))


if __name__ == "__main__":
    main()
