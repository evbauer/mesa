"""Plot measured-routing slice atlas without routing fills or EOS blending."""

import json
import os
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
OLD = ROOT / "freeeos_biquintic_progress/experiments/060_pure_freeeos_faces"
NEW = ROOT / "freeeos_biquintic_progress/experiments/071_figure3_slice_atlas"


def routing_spans(logq, routing):
    changes = np.r_[0, np.flatnonzero(routing[1:] != routing[:-1]) + 1, len(routing)]
    return [{"code": int(routing[start]), "count": int(stop - start),
             "logQ_min": float(logq[start]), "logQ_max": float(logq[stop - 1])}
            for start, stop in zip(changes[:-1], changes[1:])]


def main():
    os.environ["MPLCONFIGDIR"] = str(NEW / ".mplconfig")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    selection = json.loads((NEW / "selection.json").read_text())
    metadata = json.loads((NEW / "curve_metadata.json").read_text())
    probe = np.loadtxt(NEW / "routing_probe.tsv", skiprows=1, ndmin=2)
    cases = selection["cases"][1:]
    if probe.shape != (401 * len(cases), 14):
        raise ValueError(f"Unexpected routing probe shape {probe.shape}")
    summaries = {}
    for case_number, case in enumerate(cases):
        name = case["id"]
        data = np.load(NEW / f"{name}_curves.npz")
        measured = probe[401 * case_number:401 * (case_number + 1)]
        if not np.allclose(measured[:, :2], np.column_stack((data["logRho"], data["logT"])), rtol=0, atol=1e-12):
            raise ValueError(f"Probe order or coordinates differ from {name} slice")
        fractions = measured[:, 3:11]
        good = (measured[:, 2] == 0) & np.all(np.isfinite(measured), axis=1)
        good &= np.all(fractions >= -1e-12, axis=1) & np.all(fractions <= 1 + 1e-12, axis=1)
        good &= abs(fractions.sum(axis=1) - 1) <= 1e-12
        pure = good & (abs(fractions[:, 2] - 1) <= 1e-12)
        pure &= (np.sum(abs(np.delete(fractions, 2, axis=1)), axis=1) <= 1e-12)
        inactive = good & (abs(fractions[:, 2]) <= 1e-12)
        routing = np.where(~good, 3, np.where(pure, 0, np.where(inactive, 2, 1)))
        anchor_sample = int(np.argmin(abs(data["logQ"] - case["logQ"])))
        summaries[name] = {"successful_valid_count": int(good.sum()), "pure_FreeEOS_count": int(pure.sum()),
                           "blend_count": int((routing == 1).sum()), "inactive_count": int(inactive.sum()),
                           "failed_or_invalid_count": int((~good).sum()),
                           "supported_count": int(data["supported"].sum()),
                           "supported_and_pure_count": int((data["supported"] & pure).sum()),
                           "spans": routing_spans(data["logQ"], routing),
                           "nearest_anchor_sample_logQ": float(data["logQ"][anchor_sample]),
                           "nearest_anchor_sample_FreeEOS_fraction": float(fractions[anchor_sample, 2])}

        figure, axes = plt.subplots(4, 1, figsize=(10, 9.5), sharex=True, constrained_layout=True)
        for axis in axes:
            axis.axvline(case["logRho"], color=".35", linestyle=":", linewidth=.85)
            axis.grid(axis="y", alpha=.2)
        styles = [("legacy", "tab:orange", "--", "Legacy columns")]
        if case["orientation"] == "tau":
            styles = [("lower", "tab:blue", "-", "Biquintic: lower T cell"),
                      ("upper", "tab:purple", "-.", "Biquintic: upper T cell"), *styles]
        else:
            styles = [("single", "tab:blue", "-", "Total-potential biquintic"), *styles]
        segments = np.unique(data["trace_segment"])
        trace_density = data["trace_logQ"] + 2 * case["logT"] - 12
        for label, color, line, display_name in styles:
            for segment in segments:
                selected = data["trace_segment"] == segment
                for axis, field in zip(axes, ("Cv", "Cv_x", "Cv_t", "dse_Cv")):
                    values = data[f"trace_{label}_{field}"][selected]
                    shown = np.maximum(abs(values), 1e-16) if field == "dse_Cv" else values / 1e8
                    axis.plot(trace_density[selected], shown, color=color, linestyle=line, linewidth=1.25,
                              label=display_name if segment == 0 and field == "Cv" else None)
        axes[0].legend(loc="best", fontsize=8)
        orientation = "$u$ (density)" if case["orientation"] == "u" else r"$\tau$ (temperature)"
        axes[0].set_title(f"{name}: exact {orientation} face; $\\log_{{10}}T={case['logT']:.6g}$, "
                          f"$\\log_{{10}}Q={case['logQ']:.6g}$\n"
                          f"$|\\Delta C_{{V,t}}|/D={case['jump_abs']:.4g}$; "
                          f"{summaries[name]['pure_FreeEOS_count']}/401 wider samples route fully FreeEOS", fontsize=10)
        axes[0].set_ylabel(r"$C_V$ [$10^8$ erg g$^{-1}$ K$^{-1}$]")
        axes[1].set_ylabel(r"$C_{V,x}$ [$10^8$ erg g$^{-1}$ K$^{-1}$]")
        axes[2].set_ylabel(r"$C_{V,t}$ [$10^8$ erg g$^{-1}$ K$^{-1}$]")
        axes[3].set_ylabel(r"$|T s_T/C_V-1|$")
        axes[3].set_yscale("log")
        axes[3].set_ylim(bottom=3e-17)
        axes[3].text(.015, .05, r"Exact biquintic arithmetic zeros use $10^{-16}$ display floor",
                     transform=axes[3].transAxes, fontsize=7)
        axes[3].set_xlabel(r"$\log_{10}\rho$ [g cm$^{-3}$]; vertical dotted line: anchor face or midpoint")
        figure.savefig(NEW / f"{name}.pdf")
        figure.savefig(NEW / f"{name}.png", dpi=180)
        plt.close(figure)

    manifest = np.load(OLD / "face_manifest.npz")
    old = np.load(OLD / "face_results.npz")
    figure, axis = plt.subplots(figsize=(9, 6), constrained_layout=True)
    for direction, marker in (("u", "|"), ("tau", "_")):
        mask = old[f"{direction}_pure"] & (old[f"{direction}_support_code"] == 0)
        q, t = manifest[f"{direction}_logQ"][mask], manifest[f"{direction}_logT"][mask]
        axis.scatter(q + 2 * t - 12, t, c=".8", s=1, marker=marker, rasterized=True)
    for number, case in enumerate(selection["cases"]):
        axis.scatter([case["logRho"]], [case["logT"]], s=65,
                     marker="o" if case["orientation"] == "u" else "s", facecolors="white",
                     edgecolors="tab:blue" if case["orientation"] == "u" else "tab:purple", linewidths=1.5)
        axis.annotate(str(number), (case["logRho"], case["logT"]), xytext=(4, 4), textcoords="offset points", fontsize=8)
    axis.set(xlabel=r"$\log_{10}\rho$ [g cm$^{-3}$]", ylabel=r"$\log_{10}T$ [K]",
             title="Atlas anchors (0 original, 1--7 alternatives); pale: supported 060 pure-FreeEOS faces")
    axis.grid(alpha=.15)
    figure.savefig(NEW / "anchor_locator.pdf")
    figure.savefig(NEW / "anchor_locator.png", dpi=180)
    plt.close(figure)
    (NEW / "routing_summary.json").write_text(json.dumps(summaries, indent=2) + "\n")
    print(json.dumps({name: {key: row[key] for key in ("pure_FreeEOS_count", "blend_count", "inactive_count",
                                                      "failed_or_invalid_count", "supported_count")}
                      for name, row in summaries.items()}, indent=2))


if __name__ == "__main__":
    main()
