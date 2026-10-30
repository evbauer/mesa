"""Plot saved matched-state identity residuals; does not evaluate any EOS."""

import argparse
import json
import os
from pathlib import Path

import numpy as np


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    parser.add_argument("--model", type=int, default=100)
    args = parser.parse_args()
    output = args.summary.resolve().parent
    os.environ.setdefault("MPLCONFIGDIR", str(output / ".mplconfig"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    summary = json.loads(args.summary.read_text())
    selected = [record for record in summary["profiles"] if int(record["metadata"]["model_number"]) == args.model]
    if len(selected) != 1:
        raise ValueError("Expected exactly one named profile")
    records = sorted(selected[0]["records"], key=lambda record: record["logT"])
    temperatures = [record["logT"] for record in records]
    figure, axes = plt.subplots(2, 1, figsize=(8.0, 5.8), sharex=True,
                               gridspec_kw={"height_ratios": [3, 1]}, layout="constrained")
    for key, label, color, style in (
            ("off_component", "Legacy FreeEOS component", "tab:orange", "--"),
            ("off_final", "Legacy final EOS", "tab:orange", "-"),
            ("candidate_component", "Potential FreeEOS component", "tab:blue", "--"),
            ("candidate_final", "Potential final EOS", "tab:blue", "-")):
        residual = [max(abs(record[key]["ratios"][name]) for name in ("Dpe", "Dse", "Dsp")) for record in records]
        axes[0].semilogy(temperatures, np.maximum(residual, 1e-16), color=color, linestyle=style,
                         linewidth=1.4, label=label)
    axes[0].axhline(1e-3, color=".4", linestyle=":", linewidth=.8)
    axes[0].set_ylabel(r"$\max(|D_{pe}|,|D_{se}|,|D_{sp}|)$")
    axes[0].set_ylim(1e-16, 2e-2)
    axes[0].legend(loc="center left", fontsize=9, framealpha=.95)
    axes[0].set_title(f"Matched EOS calls at all {len(records)} zones of candidate model {args.model}")
    axes[1].plot(temperatures, [record["fraction_FreeEOS"] for record in records], color="black", linewidth=1.5)
    axes[1].set_ylabel("Final FreeEOS\nweight")
    axes[1].set_xlabel(r"$\log_{10}(T/\mathrm{K})$")
    for axis in axes:
        axis.grid(alpha=.16)
    for extension in ("pdf", "png"):
        figure.savefig(output / f"matched_profile_{args.model}.{extension}", dpi=180)
    plt.close(figure)


if __name__ == "__main__":
    run()
