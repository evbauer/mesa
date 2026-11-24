"""Plot saved model8 A/B profile differences; no stellar or EOS evaluation."""

import math
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analyze import mesa_profile
from harness import RUNS, SUMMARY, placement


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extension", action="store_true")
    arguments = parser.parse_args()
    placement()
    suffix, model = ("100", 100) if arguments.extension else ("00", 8)
    baseline_header, baseline = mesa_profile(RUNS / f"baseline_{suffix}/LOGS/profile{model}.data")
    candidate_header, candidate = mesa_profile(RUNS / f"candidate_{suffix}/LOGS/profile{model}.data")
    if baseline_header["star_age"] != candidate_header["star_age"]:
        raise ValueError("comparison profiles have different ages")
    figure, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True, constrained_layout=True)
    temperature = baseline["logT"]
    for quantity, label in (("logRho", r"$\Delta\rho/\rho$"), ("logT", r"$\Delta T/T$"), ("logR", r"$\Delta r/r$")):
        fractional = [1e6 * math.expm1(math.log(10) * (candidate_value - baseline_value))
                      for baseline_value, candidate_value in zip(baseline[quantity], candidate[quantity])]
        axes[0].plot(temperature, fractional, label=label, linewidth=1.1)
    axes[0].set_ylabel("Structural difference (ppm)")
    axes[0].legend(loc="best", ncol=3)
    for quantity, label in (("cp", r"$C_P$"), ("cv", r"$C_V$")):
        fractional = [1e6 * (candidate_value - baseline_value) / baseline_value
                      for baseline_value, candidate_value in zip(baseline[quantity], candidate[quantity])]
        axes[1].plot(temperature, fractional, label=label, linewidth=1.1)
    axes[1].set_ylabel("Heat-capacity difference (ppm)")
    axes[1].set_xlabel(r"Baseline $\log_{10} T$ (K)")
    axes[1].legend(loc="best")
    for axis in axes:
        axis.axhline(0, color="black", linewidth=.5)
        axis.grid(alpha=.2)
    figure.suptitle(f"Exploratory frozen-composition relaxation, model{model}\nCandidate minus baseline on the same mass grid and at the same age", fontsize=11)
    output = "extension_profile_comparison" if arguments.extension else "profile_comparison"
    figure.savefig(SUMMARY / (output + ".pdf"))
    figure.savefig(SUMMARY / (output + ".png"), dpi=150)
    plt.close(figure)


if __name__ == "__main__":
    main()
