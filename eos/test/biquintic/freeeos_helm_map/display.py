"""Recover old 060 plotted axes and summarize the corrected 070 display crop."""

import json

import numpy as np

from common import OLD, OUTPUT, digest, save_json


def old_axes():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 2, sharex=True, sharey=True)
    with np.load(OLD / "face_manifest.npz") as manifest, np.load(OLD / "face_results.npz") as results:
        for axis, direction in zip(axes, ("u", "tau")):
            face_q = manifest[f"{direction}_logQ"]
            face_t = manifest[f"{direction}_logT"]
            density = face_q + 2 * face_t - 12
            pure = results[f"{direction}_pure"]
            support = results[f"{direction}_support_code"] == 0
            jump = results[f"{direction}_Cv_t_jump"]
            gap = pure & ~support
            shown = pure & support & np.isfinite(jump)
            axis.scatter(density[gap], face_t[gap], s=2.2)
            axis.scatter(density[shown], face_t[shown], s=2.2)
        old_density = list(map(float, axes[0].get_xlim()))
        old_temperature = list(map(float, axes[0].get_ylim()))
    plt.close(figure)
    return old_density, old_temperature


def in_display(grid, bounds, ye):
    density = np.log10(grid["din"] / ye)
    temperature = np.log10(grid["temperature"])
    density_ok = (density >= bounds["new_logRho"][0]) & (density <= bounds["new_logRho"][1])
    temperature_ok = (temperature >= bounds["new_logT"][0]) & (temperature <= bounds["new_logT"][1])
    return density_ok[:, None] & temperature_ok[None, :]


def main():
    broad = json.loads((OUTPUT / "bounds.json").read_text())
    old_density, old_temperature = old_axes()
    corrected = {"old_logRho": old_density, "old_logT": old_temperature,
                 "new_logRho": [old_density[0], old_density[1] + 2],
                 "new_logT": [old_temperature[0], old_temperature[1] + 2],
                 "Ye": broad["Ye"],
                 "method": "Replayed 060's pure-supported and pure-gap scatter coordinates with shared axes and Matplotlib default 5% margins; autoscaled limits from get_xlim/get_ylim; +2.00 dex on upper limits",
                 "broad_survey_manifest_bounds": broad,
                 "old_manifest_sha256": digest(OLD / "face_manifest.npz"),
                 "old_results_sha256": digest(OLD / "face_results.npz")}
    if (OUTPUT / "display_bounds_corrected.json").exists():
        raise FileExistsError("corrected display bounds already saved")
    save_json(OUTPUT / "display_bounds_corrected.json", corrected)
    stats = {"display_bounds": corrected, "HELM": {}, "FreeEOS": {}}
    for direction in ("T", "d"):
        filename = "helm_T_traces_corrected.npz" if direction == "T" else "helm_d_traces.npz"
        with np.load(OUTPUT / f"helm_{direction}_grid.npz") as grid, np.load(OUTPUT / filename) as traces:
            crop = in_display(grid, corrected, broad["Ye"])
            pure = traces["pure"]
            route_failure = ~traces["route_success"] | ~traces["fraction_good"]
            external = ~pure & ~route_failure
            kernel_failure = external & ~traces["kernel_good"]
            shown = crop & traces["shown"]
            direction_stats = {"native_faces": int(crop.sum()),
                "direct_fully_FreeEOS": int((crop & pure).sum()),
                "routing_failure": int((crop & route_failure).sum()),
                "external_HELM_evaluator_failure": int((crop & kernel_failure).sum()),
                "shown": int(shown.sum()),
                "shown_exact_zero": int((shown & (traces["jump"] == 0)).sum()),
                "shown_below_display": int((shown & (np.abs(traces["jump"]) < 1e-5)).sum()),
                "shown_above_display": int((shown & (np.abs(traces["jump"]) > 10)).sum()),
                "minimum_denominator": float(np.min(traces["denominator"][shown])) if shown.any() else None}
            if direction == "T":
                direction_stats["precision_unresolved"] = int((crop & external & traces["kernel_good"] & ~traces["resolved"]).sum())
                selected = np.abs(traces["jump"][shown])
                direction_stats["jump"] = {"median": float(np.median(selected)),
                                            "p99": float(np.quantile(selected, .99)),
                                            "maximum": float(np.max(selected))}
            stats["HELM"][direction] = direction_stats
    with np.load(OUTPUT / "freeeos_traces.npz") as traces:
        for direction in ("u", "tau"):
            density = traces[f"{direction}_logQ"] + 2 * traces[f"{direction}_logT"] - 12
            temperature = traces[f"{direction}_logT"]
            crop = (density >= corrected["new_logRho"][0]) & (density <= corrected["new_logRho"][1]) & (
                temperature >= corrected["new_logT"][0]) & (temperature <= corrected["new_logT"][1])
            pure = traces[f"{direction}_pure"]
            support = traces[f"{direction}_support_code"] == 0
            stats["FreeEOS"][direction] = {"pure": int((crop & pure).sum()),
                "supported": int((crop & pure & support).sum()),
                "support_gap": int((crop & pure & ~support).sum())}
    save_json(OUTPUT / "in_display_stats_corrected.json", stats)
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
