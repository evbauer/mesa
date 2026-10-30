"""Held-out-node validation of the WDM_AV Hermite tables.

For a source table on the fine grid (20 points/decade), build processed tables
on the *coarse* grid of every other node (10 points/decade) with the builder,
evaluate them with the independent :mod:`hermite_eval` evaluator at the omitted
fine nodes, and compare against the source's own values there.

Three build variants are compared:

* ``a_o5`` -- biquintic, second derivatives from the source dPdRho_T/dPdT_Rho/Cv;
* ``b_o5`` -- biquintic, second derivatives by finite differences
  (``--ignore-source-derivs``);
* ``a_o3`` -- bicubic (second derivatives come from the cubic interpolant).

Physical quantities and their source truths (eos_mixing_notes.tex Sec. 2)::

    P       = rho f_xi                 <->  source P
    s       = -f_tau / T               <->  source S
    c_V     = (f_tau - f_tautau) / T   <->  source Cv
    chi_rho = 1 + f_xixi / f_xi        <->  rho dPdRho_T / P
    chi_T   = f_xitau / f_xi           <->  T   dPdT_Rho / P
    e       = f - f_tau                <->  source E

A fine node is counted only if it is omitted from the coarse grid, its coarse
cell is fully valid (all four coarse corners valid), and the source truth is
finite and non-zero.  Statistics are reported over four regions:

* ``skye_valid``          -- logT >= 6.5 OR logRho >= 3.7 (Decision D20);
* ``skye_valid_no_pair``  -- the above, excluding the high-T / low-rho
  pair-radiation corner (NOT (logT > 8.0 AND logRho < 0));
* ``target``              -- the WDM_AV coverage hole, logRho in [-2, 4] and
  logT in [3, 7];
* ``whole_valid``         -- the whole valid plane.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

# make the builder package and this directory importable when run as a script
_HERE = Path(__file__).resolve().parent
for _p in (str(_HERE), str(_HERE.parent)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import hermite_eval as he
import numpy as np
from wdm_av_builder import constants as C
from wdm_av_builder import io as wio
from wdm_av_builder.derivs import build_nodal_data

VARIANTS: dict[str, dict] = {
    "a_o5": {"order": 5, "ignore_source_derivs": False},
    "b_o5": {"order": 5, "ignore_source_derivs": True},
    "a_o3": {"order": 3, "ignore_source_derivs": False},
}
VARIANT_LABELS = {
    "a_o5": "(a) biquintic, source 2nd derivs",
    "b_o5": "(b) biquintic, FD 2nd derivs",
    "a_o3": "(a) bicubic",
}
QUANTITIES = ("P", "s", "c_V", "chi_rho", "chi_T", "e")

# statistics regions (evaluated on the fine (logRho, logT) grid).  "skye_valid"
# and "whole_valid" are the original two; "target" is the WDM_AV coverage hole
# (Decision D20); "skye_valid_no_pair" drops the high-T / low-rho pair-radiation
# corner where P is nearly flat in rho and chi_rho cancels.
REGION_NAMES = ("skye_valid", "skye_valid_no_pair", "target", "whole_valid")


def _regions(logrho2d: np.ndarray, logt2d: np.ndarray) -> dict[str, np.ndarray]:
    skye = (logt2d >= 6.5) | (logrho2d >= 3.7)
    pair_corner = (logt2d > 8.0) & (logrho2d < 0.0)
    all_true = np.ones_like(skye, dtype=bool)
    return {
        "skye_valid": skye,
        "skye_valid_no_pair": skye & ~pair_corner,
        "target": (logrho2d >= -2.0)
        & (logrho2d <= 4.0)
        & (logt2d >= 3.0)
        & (logt2d <= 7.0),
        "whole_valid": all_true,
    }


DEFAULT_FIGDIR = (
    Path(os.environ.get("WDM_AV_MESA_DIR") or _HERE.parents[2]).parent
    / "claude_tex_notes/figures/heldout"
)
DEFAULT_RESULTSDIR = Path(__file__).resolve().parent / "results"


def coarsen(source: wio.SourceTable) -> wio.SourceTable:
    """Return the every-other-node coarsening of a source table."""
    return wio.SourceTable(
        header=dict(source.header),
        logrho=source.logrho[::2].copy(),
        logt=source.logt[::2].copy(),
        columns={k: v[::2, ::2].copy() for k, v in source.columns.items()},
    )


def _coarse_table(source_coarse: wio.SourceTable, variant: str) -> he.HermiteTable:
    opts = VARIANTS[variant]
    nodal = build_nodal_data(
        source_coarse,
        order=opts["order"],
        ignore_source_derivs=opts["ignore_source_derivs"],
    )
    return he.table_from_arrays(nodal.attrs, nodal.datasets)


def _truth(source: wio.SourceTable):
    """Source-exact physical quantities on the fine grid (shape (ntau, nxi))."""
    rho = np.exp(C.ln10 * source.logrho)[None, :]
    T = np.exp(C.ln10 * source.logt)[:, None]
    cols = source.columns
    P = cols["P"]
    return {
        "P": P,
        "s": cols["S"],
        "c_V": cols["Cv"],
        "chi_rho": rho * cols["dPdRho_T"] / P,
        "chi_T": T * cols["dPdT_Rho"] / P,
        "e": cols["E"],
    }


def _computed(table: he.HermiteTable, source: wio.SourceTable):
    """Evaluator-derived physical quantities at every fine node, plus the coarse
    cell indices used (for validity gating)."""
    xi = C.ln10 * source.logrho
    tau = C.ln10 * source.logt
    XX, TT = np.meshgrid(xi, tau)  # (ntau, nxi)
    out = he.evaluate(table, XX, TT)
    rho = np.exp(XX)
    T = np.exp(TT)
    comp = {
        "P": rho * out["f_xi"],
        "s": -out["f_tau"] / T,
        "c_V": (out["f_tau"] - out["f_tautau"]) / T,
        "chi_rho": 1.0 + out["f_xixi"] / out["f_xi"],
        "chi_T": out["f_xitau"] / out["f_xi"],
        "e": out["f"] - out["f_tau"],
    }
    return comp, out["ixi"], out["itau"]


def _omitted_mask(ntau: int, nxi: int) -> np.ndarray:
    """Fine nodes not present on the coarse (every-other) grid."""
    j = np.arange(ntau)[:, None]
    i = np.arange(nxi)[None, :]
    in_coarse = (i % 2 == 0) & (j % 2 == 0)
    return ~in_coarse


def _cellcenter_mask(ntau: int, nxi: int) -> np.ndarray:
    """Cell-center fine nodes: odd index in BOTH directions.

    These sit at the geometric centre of a coarse cell, the farthest from every
    coarse node, so they are a uniform, maximally-interpolated subset of the
    held-out nodes.  Plotting only these removes the alternating stripe pattern
    that the edge-centred held-out nodes (odd in one direction only, lying on a
    coarse grid line) produce on the error maps.
    """
    j = np.arange(ntau)[:, None]
    i = np.arange(nxi)[None, :]
    return (i % 2 == 1) & (j % 2 == 1)


def _cell_valid(valid_coarse: np.ndarray, ixi: np.ndarray, itau: np.ndarray):
    """True where all four coarse corners of the enclosing cell are valid."""
    v = valid_coarse.astype(bool)
    return v[itau, ixi] & v[itau, ixi + 1] & v[itau + 1, ixi] & v[itau + 1, ixi + 1]


def _stats(relerr: np.ndarray, mask: np.ndarray) -> dict:
    a = relerr[mask]
    a = a[np.isfinite(a)]
    if a.size == 0:
        return {"count": 0, "median": None, "p90": None, "p99": None, "max": None}
    return {
        "count": int(a.size),
        "median": float(np.median(a)),
        "p90": float(np.percentile(a, 90)),
        "p99": float(np.percentile(a, 99)),
        "max": float(a.max()),
    }


def run_one(
    source_path: Path,
    label: str,
    *,
    figdir: Path = DEFAULT_FIGDIR,
    variants=tuple(VARIANTS),
) -> dict:
    """Run the held-out validation for one source over the requested variants.

    Returns a nested summary dict and writes one error-map PNG per variant.
    """
    source = wio.read_source(source_path)
    required = ("P", "S", "Cv", "dPdRho_T", "dPdT_Rho", "E")
    missing = [c for c in required if c not in source.columns]
    if missing:
        raise ValueError(f"{source_path}: source missing columns {missing}")

    ntau, nxi = source.nlogt, source.nlogrho
    logrho, logt = source.logrho, source.logt
    logrho2d, logt2d = np.meshgrid(logrho, logt)
    omitted = _omitted_mask(ntau, nxi)
    cellcenter = _cellcenter_mask(ntau, nxi)
    regions = _regions(logrho2d, logt2d)

    truth = _truth(source)
    figdir = Path(figdir)

    summary: dict = {
        "label": label,
        "source_file": str(source_path),
        "isotope": source.name,
        "nfine": [ntau, nxi],
        "variants": {},
    }

    for variant in variants:
        coarse_src = coarsen(source)
        table = _coarse_table(coarse_src, variant)
        valid_coarse = table.datasets["valid"]
        comp, ixi, itau = _computed(table, source)
        cellok = _cell_valid(valid_coarse, ixi, itau)
        base = omitted & cellok

        vsum: dict = {"label": VARIANT_LABELS[variant], "quantities": {}}
        maps: dict[str, np.ndarray] = {}
        for q in QUANTITIES:
            tr = truth[q]
            cp = comp[q]
            good = base & np.isfinite(tr) & (tr != 0.0) & np.isfinite(cp)
            relerr = np.full((ntau, nxi), np.nan)
            relerr[good] = np.abs(cp[good] - tr[good]) / np.abs(tr[good])
            maps[q] = relerr
            qstat = {
                name: _stats(relerr, good & regions[name]) for name in REGION_NAMES
            }
            # cell-center-only statistics (odd index in both directions): the
            # uniform, maximally-interpolated held-out subset drawn on the maps.
            qstat["cell_center"] = {
                name: _stats(relerr, good & cellcenter & regions[name])
                for name in REGION_NAMES
            }
            # location of the worst Skye-region error, with the values there
            region = good & regions["skye_valid"]
            if region.any():
                r = np.where(region, relerr, -np.inf)
                jmax, imax = np.unravel_index(np.argmax(r), r.shape)
                qstat["worst_skye_at"] = {
                    "logRho": float(logrho[imax]),
                    "logT": float(logt[jmax]),
                    "relerr": float(relerr[jmax, imax]),
                    "truth": float(tr[jmax, imax]),
                    "interp": float(cp[jmax, imax]),
                }
            vsum["quantities"][q] = qstat
        summary["variants"][variant] = vsum

        _plot_variant(label, variant, logrho, logt, maps, figdir)

    return summary


def _plot_variant(label, variant, logrho, logt, maps, figdir: Path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    # Uniform error-map convention (all WDM_AV error maps): log10 relative error on
    # a diverging RdBu_r scale centred at -3, blue (better) to -7, red (worse) to
    # +1, saturating beyond both ends; no-data cells masked and drawn light grey.
    norm = TwoSlopeNorm(vmin=-7.0, vcenter=-3.0, vmax=1.0)
    cmap = plt.get_cmap("coolwarm").copy()
    cmap.set_under(cmap(0.0))
    cmap.set_over(cmap(1.0))
    cmap.set_bad("white")

    # cell-center subgrid (odd index in both directions): a regular grid, so the
    # map is continuous with no white gaps between the plotted held-out nodes.
    lr_cc = logrho[1::2]
    lt_cc = logt[1::2]

    figdir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5))
    axes = axes.ravel()
    im = None
    for ax, q in zip(axes, QUANTITIES):
        with np.errstate(divide="ignore"):
            z = np.log10(maps[q][1::2, 1::2])  # NaN stays NaN; exact 0 -> -inf
        z = np.ma.masked_invalid(z)
        ax.set_facecolor("white")
        im = ax.pcolormesh(lr_cc, lt_cc, z, shading="auto", cmap=cmap, norm=norm)
        # outline the no-data regions so white cannot be read as the gray midpoint
        _m = np.ma.getmaskarray(z).astype(float)
        if _m.any() and not _m.all():
            ax.contour(lr_cc, lt_cc, _m, levels=[0.5], colors="0.3", linewidths=0.5)
        ax.axhline(6.5, color="k", lw=0.6, ls="--", alpha=0.5)
        ax.axvline(3.7, color="k", lw=0.6, ls="--", alpha=0.5)
        ax.set_title(f"{q}  log10|rel err|", fontsize=9)
        ax.set_xlabel("logRho")
        ax.set_ylabel("logT")
    fig.suptitle(f"held-out error: {label}  {VARIANT_LABELS[variant]}")
    fig.tight_layout(rect=(0, 0, 0.92, 1))
    cbar = fig.colorbar(im, ax=axes.tolist(), extend="both", fraction=0.046, pad=0.02)
    cbar.set_label("log10 |relative error|")
    path = figdir / f"{label}_{variant}.png"
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------
def _label_for(source_path: Path) -> str:
    """Auto label: <skye|helm|...>_<isotope>, from the parent dir and file stem."""
    parent = source_path.parent.name
    tag = parent.replace("source_", "")
    iso = source_path.stem.replace("eosWDM_AV_", "")
    return f"{tag}_{iso}"


def write_reports(summaries: list[dict], resultsdir: Path) -> dict[str, Path]:
    resultsdir = Path(resultsdir)
    resultsdir.mkdir(parents=True, exist_ok=True)
    json_path = resultsdir / "heldout_summary.json"
    with json_path.open("w") as fh:
        json.dump(summaries, fh, indent=2)

    # flat CSV: one row per (label, variant, quantity, region)
    csv_path = resultsdir / "heldout_summary.csv"
    fields = [
        "label",
        "isotope",
        "variant",
        "quantity",
        "node_set",
        "region",
        "count",
        "median",
        "p90",
        "p99",
        "max",
    ]
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for s in summaries:
            for variant, vsum in s["variants"].items():
                for q, qstat in vsum["quantities"].items():
                    for region in REGION_NAMES:
                        w.writerow(
                            {
                                "label": s["label"],
                                "isotope": s["isotope"],
                                "variant": variant,
                                "quantity": q,
                                "node_set": "all",
                                "region": region,
                                **qstat[region],
                            }
                        )
                    for region in REGION_NAMES:
                        w.writerow(
                            {
                                "label": s["label"],
                                "isotope": s["isotope"],
                                "variant": variant,
                                "quantity": q,
                                "node_set": "cellcenter",
                                "region": region,
                                **qstat["cell_center"][region],
                            }
                        )

    md_path = resultsdir / "heldout_summary.md"
    md_path.write_text(_markdown_table(summaries))
    return {"json": json_path, "csv": csv_path, "md": md_path}


def _fmt(x) -> str:
    if x is None:
        return "  -  "
    return f"{x:.1e}"


_REGION_DESC = {
    "skye_valid": "Skye-valid region: logT >= 6.5 OR logRho >= 3.7",
    "skye_valid_no_pair": (
        "Skye-valid region excluding the pair-radiation corner "
        "(NOT (logT > 8.0 AND logRho < 0))"
    ),
    "target": "WDM_AV target region: logRho in [-2, 4] AND logT in [3, 7]",
    "whole_valid": "whole valid plane",
}


def _markdown_table(summaries: list[dict]) -> str:
    """median/p99/max for the five primary quantities, one block per region."""
    qs = ("P", "s", "c_V", "chi_rho", "chi_T")
    blurb = (
        "Relative error of each quantity at fine nodes omitted from the "
        "10-points/decade coarse grid, vs the source values. "
        "Variants: a_o5 = biquintic/source 2nd derivs, b_o5 = biquintic/FD "
        "2nd derivs, a_o3 = bicubic. Full per-region data (incl. p90 and the "
        "whole valid plane) is in heldout_summary.csv / .json."
    )
    lines = ["# Held-out-node interpolation error", "", blurb, ""]
    for region in ("skye_valid", "skye_valid_no_pair", "target"):
        lines.append(f"## {_REGION_DESC[region]}")
        lines.append("")
        lines.append("| table | variant | stat | " + " | ".join(qs) + " |")
        lines.append("|" + "---|" * (3 + len(qs)))
        for s in summaries:
            for variant in ("a_o5", "b_o5", "a_o3"):
                if variant not in s["variants"]:
                    continue
                vsum = s["variants"][variant]
                for stat in ("median", "p99", "max"):
                    row = [s["label"], variant, stat]
                    for q in qs:
                        row.append(_fmt(vsum["quantities"][q][region][stat]))
                    lines.append("| " + " | ".join(row) + " |")
        lines.append("")
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m heldout",
        description="held-out-node validation of WDM_AV Hermite tables",
    )
    p.add_argument("sources", nargs="+", help="source table paths")
    p.add_argument("--figdir", default=str(DEFAULT_FIGDIR))
    p.add_argument("--resultsdir", default=str(DEFAULT_RESULTSDIR))
    p.add_argument(
        "--variants",
        nargs="+",
        default=list(VARIANTS),
        choices=list(VARIANTS),
    )
    p.add_argument("--label", default=None, help="label (single source only)")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    summaries = []
    for src in args.sources:
        src_path = Path(src)
        label = (
            args.label
            if (args.label and len(args.sources) == 1)
            else _label_for(src_path)
        )
        print(f"heldout: {label}  ({src_path})")
        summary = run_one(
            src_path, label, figdir=Path(args.figdir), variants=tuple(args.variants)
        )
        summaries.append(summary)
        for variant in args.variants:
            sk = summary["variants"][variant]["quantities"]
            msg = "  ".join(
                f"{q}={_fmt(sk[q]['skye_valid']['median'])}" for q in QUANTITIES
            )
            print(f"    {variant}: median(skye) {msg}")
    paths = write_reports(summaries, Path(args.resultsdir))
    print("wrote", paths["json"], paths["csv"], paths["md"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
