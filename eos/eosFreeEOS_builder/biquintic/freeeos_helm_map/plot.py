"""Plot only stored 070 arrays; never evaluate an EOS during a redraw."""

import argparse
import json
import os
from pathlib import Path

import numpy as np

from common import OLD, OUTPUT, digest, save_json
from display import in_display


def extent(density, temperature, ye):
    log_density = np.log10(density / ye)
    log_temperature = np.log10(temperature)
    density_step = float(np.median(np.diff(log_density)))
    temperature_step = float(np.median(np.diff(log_temperature)))
    return (float(log_density[0] - density_step / 2), float(log_density[-1] + density_step / 2),
            float(log_temperature[0] - temperature_step / 2),
            float(log_temperature[-1] + temperature_step / 2))


def sampled_boundary(axis, manifest, pure, direction):
    density_index = manifest[f"{direction}_density_index"]
    temperature_index = manifest[f"{direction}_temperature_index"]
    shape = (int(density_index.max()) + 1, int(temperature_index.max()) + 1)
    face_q = manifest[f"{direction}_logQ"].reshape(shape)
    face_t = manifest[f"{direction}_logT"].reshape(shape)
    q_axis = face_q[:, 0]
    t_axis = face_t[0]
    padded_q = np.concatenate(([q_axis[0] - (q_axis[1] - q_axis[0])], q_axis,
                               [q_axis[-1] + (q_axis[-1] - q_axis[-2])]))
    padded_t = np.concatenate(([t_axis[0] - (t_axis[1] - t_axis[0])], t_axis,
                               [t_axis[-1] + (t_axis[-1] - t_axis[-2])]))
    q_mesh, t_mesh = np.meshgrid(padded_q, padded_t, indexing="ij")
    density = q_mesh + 2 * t_mesh - 12
    selection = np.pad(pure.reshape(shape), 1, constant_values=False).astype(float)
    axis.contourf(density, t_mesh, selection, levels=(.5, 1.5), colors=("#fff7e6",), alpha=.8)
    axis.contour(density, t_mesh, selection, levels=(.5,), colors="black", linewidths=4.0,
                 zorder=10)
    return {"native_logQ_step": float(np.median(np.diff(q_axis))),
            "native_logT_step": float(np.median(np.diff(t_axis))),
            "interpolation": "0.5 contour of zero-padded binary midpoint selection mask; not exact MESA routing geometry"}


def categorical_overlay(axis, category, bounds, color, opacity=1):
    import matplotlib.colors as colors

    masked = np.ma.masked_where(~category.T, category.T.astype(float))
    palette = colors.ListedColormap([color])
    palette.set_bad((0, 0, 0, 0))
    axis.imshow(masked, origin="lower", interpolation="nearest", extent=bounds,
                aspect="auto", cmap=palette, vmin=0, vmax=1, alpha=opacity, rasterized=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stem", default="freeeos_helm_native_faces_corrected")
    args = parser.parse_args()
    if Path(args.stem).name != args.stem or not args.stem:
        raise ValueError("stem must be a simple filename")
    figure_paths = [OUTPUT / f"{args.stem}.{extension}" for extension in ("png", "pdf")]
    presentation_path = OUTPUT / f"{args.stem}_presentation.json"
    if any(path.exists() for path in (*figure_paths, presentation_path)):
        raise FileExistsError("choose a new plot stem; previous artifacts are immutable")
    os.environ.setdefault("MPLCONFIGDIR", str(OUTPUT / "mpl_cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    bounds = json.loads((OUTPUT / "display_bounds_corrected.json").read_text())
    ye = bounds["Ye"]
    norm = TwoSlopeNorm(vmin=-5, vcenter=-2, vmax=1)
    cmap = plt.get_cmap("coolwarm").copy()
    cmap.set_bad((0, 0, 0, 0))
    figure, axes = plt.subplots(1, 2, figsize=(14.8, 6.0), sharex=True, sharey=True,
                                constrained_layout=True)
    presentation = {"cmap": "coolwarm", "vmin": -5, "vcenter": -2, "vmax": 1,
                    "display_zero_floor": 1e-17, "bounds": bounds,
                    "in_display_stats_sha256": digest(OUTPUT / "in_display_stats_corrected.json"),
                    "inputs_sha256": {name: digest(OUTPUT / name) for name in
                                      ("freeeos_traces.npz", "helm_T_traces_corrected.npz",
                                       "helm_d_traces.npz", "helm_T_grid.npz", "helm_d_grid.npz")},
                    "panels": {}}
    with np.load(OLD / "face_manifest.npz") as manifest, \
         np.load(OUTPUT / "freeeos_traces.npz") as freeeos:
        for axis, free_direction, helm_direction, label in zip(
                axes, ("u", "tau"), ("d", "T"),
                ("FreeEOS $u$ faces / HELM $d$ faces", "FreeEOS $\\tau$ faces / HELM $T$ faces")):
            axis.set_facecolor("#f5f5f5")
            old_pure = freeeos[f"{free_direction}_pure"]
            boundary = sampled_boundary(axis, manifest, old_pure, free_direction)
            with np.load(OUTPUT / f"helm_{helm_direction}_grid.npz") as grid, \
                 np.load(OUTPUT / ("helm_T_traces_corrected.npz" if helm_direction == "T" else
                                   "helm_d_traces.npz")) as helm:
                shown = helm["shown"]
                metric = np.abs(helm["jump"])
                helm_bounds = extent(grid["din"], grid["temperature"], ye)
                crop = in_display(grid, bounds, ye)
                log_metric = np.ma.masked_where(~shown, np.log10(np.maximum(metric, 1e-17)))
                image = axis.imshow(log_metric.T, origin="lower", interpolation="nearest",
                                    extent=helm_bounds, aspect="auto", cmap=cmap, norm=norm,
                                    rasterized=True)
                categorical_overlay(axis, helm["pure"], helm_bounds, "#fff0cb")
                failed = ~helm["route_success"] | ~helm["fraction_good"]
                failed |= ~helm["kernel_good"] & ~helm["pure"]
                if helm_direction == "T":
                    categorical_overlay(axis, ~helm["resolved"] & ~helm["pure"] &
                                        helm["kernel_good"] & helm["route_success"], helm_bounds, "#8c70a2")
                categorical_overlay(axis, helm["branch_near"] & ~helm["pure"], helm_bounds, "#b485a9")
                categorical_overlay(axis, failed, helm_bounds, "#b8a5ae")
                helm_counts = {"shown": int((shown & crop).sum()), "failed": int((failed & crop).sum()),
                               "zero_shown": int(np.count_nonzero(shown & crop & (metric == 0)))}
            face_q = freeeos[f"{free_direction}_logQ"]
            face_t = freeeos[f"{free_direction}_logT"]
            pure = freeeos[f"{free_direction}_pure"]
            supported = freeeos[f"{free_direction}_support_code"] == 0
            metric = np.abs(freeeos[f"{free_direction}_jump"])
            selected = pure & supported & np.isfinite(metric)
            gap = pure & ~supported
            axis.scatter((face_q + 2 * face_t - 12)[selected], face_t[selected],
                         c=np.log10(np.maximum(metric[selected], 1e-17)), s=2.8, marker="s",
                         linewidths=0, cmap=cmap, norm=norm, rasterized=True, zorder=3)
            axis.scatter((face_q + 2 * face_t - 12)[gap], face_t[gap],
                         c="#777777", s=3.2, marker="s", linewidths=0, rasterized=True, zorder=4)
            axis.set_title(label, fontsize=11)
            axis.set_xlim(bounds["new_logRho"])
            axis.set_ylim(bounds["new_logT"])
            axis.set_xlabel(r"$\log_{10}\rho$ [g cm$^{-3}$]")
            axis.grid(alpha=.1, linewidth=.5)
            presentation["panels"][free_direction] = {"boundary": boundary,
                "FreeEOS_shown": int(selected.sum()), "FreeEOS_support_gap": int(gap.sum()),
                "HELM": helm_counts}
    axes[0].set_ylabel(r"$\log_{10}T$ [K]")
    figure.colorbar(image, ax=axes, label=r"$\log_{10}J((C_V)_{\ln T})$",
                    extend="both", shrink=.81)
    handles = [Line2D([0], [0], color="black", lw=4.0, label="Inside: FreeEOS native faces"),
               Patch(facecolor="#3b4db4", label="HELM $d$-face exact zero"),
               Patch(facecolor="#fff0cb", label="Fully FreeEOS; native HELM hidden"),
               Patch(facecolor="#777777", label="FreeEOS reconstruction unavailable"),
               Patch(facecolor="#b8a5ae", label="HELM / routing failure"),
               Patch(facecolor="#8c70a2", label="HELM precision unresolved"),
               Patch(facecolor="#f5f5f5", label="No native HELM support")]
    figure.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=8,
                  bbox_to_anchor=(.5, -.095))
    figure.suptitle("Unblended native-face $C_V$ temperature-partial jumps; solar mixture\n"
                    "HELM where direct MESA face probes are not fully FreeEOS", fontsize=11)
    figure.savefig(figure_paths[0], dpi=190, bbox_inches="tight")
    figure.savefig(figure_paths[1], bbox_inches="tight")
    save_json(presentation_path, presentation)


if __name__ == "__main__":
    main()
