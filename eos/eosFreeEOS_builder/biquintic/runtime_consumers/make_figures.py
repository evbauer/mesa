"""Plot saved precision-control evidence without reevaluating an EOS."""

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from run_consumers import OUTPUT


def main():
    failures = json.loads((OUTPUT / "precision_control/frozen_failures.json").read_text())
    failure = next(record for record in failures if record["route"] == "component")
    actual = [json.loads(line) for line in (OUTPUT / "precision_control/representable_states.jsonl").read_text().splitlines()]
    actual = [record for record in actual if record["route"] == "component"]
    matched = json.loads((OUTPUT / "precision_matched/samples.json").read_text())
    traces = json.loads((OUTPUT / "precision_matched/traces.json").read_text())
    center = failure["newton"]["history"][-1]["y"]
    offsets = (np.array(matched["y"]) - center) / np.spacing(center)
    figure, axes = plt.subplots(1, 2, figsize=(10.5, 3.5), constrained_layout=True)
    axes[0].plot(offsets, np.array([record["residual"] for record in actual]) * 1e13,
                 ".", markersize=2.5, color="firebrick", alpha=.7, label="Actual API")
    axes[0].plot(offsets, np.array(matched["residuals"]["extended"]) * 1e13,
                 color="navy", linewidth=1.5, label="Extended cell control")
    axes[0].axhspan(-1, 1, color="gray", alpha=.15, label="Frozen tolerance")
    axes[0].set(xlabel="Representable-state offset from stalled iterate (ULPs)",
                ylabel="Residual / 1e-13", title="T0980: one native cell, pure FreeEOS")
    axes[0].legend(fontsize=8)
    histories = [(failure["newton"]["history"], "Actual API", "firebrick", "-"),
                 (traces["binary64"]["history"], "Polynomial binary64 control", "darkorange", "--"),
                 (traces["extended"]["history"], "Same polynomial, extended", "navy", "-")]
    for history, label, color, linestyle in histories:
        axes[1].semilogy(range(len(history)), [abs(record["residual"]) for record in history],
                         marker=".", color=color, linestyle=linestyle, label=label)
    axes[1].axhline(1e-13, color="gray", linestyle=":")
    axes[1].set(xlabel="Newton iteration", ylabel="Absolute residual", title="Original forcing and stopping threshold")
    axes[1].legend(fontsize=8)
    for axis in axes:
        axis.grid(alpha=.2)
    figure.savefig(OUTPUT / "precision_matched/precision_diagnostic.pdf")
    figure.savefig(OUTPUT / "precision_matched/precision_diagnostic.png", dpi=150)
    plt.close(figure)


if __name__ == "__main__":
    main()
