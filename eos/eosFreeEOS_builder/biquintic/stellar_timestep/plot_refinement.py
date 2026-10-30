"""Plot saved matched-endpoint contrasts, with no EOS or stellar calculation."""

import json
import math
import time

import timestep as task
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    task.frozen.placement()
    if time.time() >= task.NUMERICAL_STOP:
        raise RuntimeError("TASK057numerical analysis cutoff20:30 reached")
    data = json.loads((task.ANALYSIS / "endpoint_vectors.json").read_text())
    comparison = json.loads((task.ANALYSIS / "comparison.json").read_text())
    if not comparison["exactly_identical_endpoint_ages"]:
        raise ValueError("figure requires exactly identical ages")
    figure, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True, constrained_layout=True)
    labels = (("baseline_refinement", "Baseline: half - coarse", "--"),
              ("candidate_refinement", "Candidate: half - coarse", ":"),
              ("EOS_coarse", "EOS: candidate - baseline, coarse", "-"),
              ("EOS_half", "EOS: candidate - baseline, half", "-"))
    for key, label, style in labels:
        axes[0].plot(data["coarse_baseline_logT"], [1e6 * math.expm1(math.log(10) * value)
                     for value in data["differences"]["logRho"][key]], label=label, linestyle=style, linewidth=1.1)
        axes[1].plot(data["coarse_baseline_logT"], [1e6 * value / reference for value, reference in
                     zip(data["differences"]["cv"][key], data["coarse_baseline_cv"])], linestyle=style, linewidth=1.1)
    axes[0].set_ylabel("Density contrast (ppm)")
    axes[1].set_ylabel(r"$\Delta C_V/C_{V,\mathrm{baseline,coarse}}$ (ppm)")
    axes[1].set_xlabel(r"Coarse baseline $\log_{10}T$ (K)")
    axes[0].legend(fontsize=8, loc="best")
    for axis in axes:
        axis.axhline(0, color="black", linewidth=.5)
        axis.grid(alpha=.2)
    figure.suptitle("Frozen-composition temporal refinement at the same endpoint\n100/1000yr versus50/500yr initial/max dt; no tolerance changes", fontsize=11)
    figure.savefig(task.ANALYSIS / "endpoint_contrasts.pdf")
    figure.savefig(task.ANALYSIS / "endpoint_contrasts.png", dpi=150)
    plt.close(figure)


if __name__ == "__main__":
    main()
