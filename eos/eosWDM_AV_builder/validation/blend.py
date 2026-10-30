"""Blend-zone analysis of the integrated WDM_AV eosDT component.

Reads the 1-D blend-edge sweeps and the 2-D fraction map written by
``test_fortran/blend_sweep.f90`` (the FULL eosDT cascade, ``eosDT_get``), and
for each cut produces:

* a **value figure** -- lnPgas, lnE, lnS, chiT, Cv, grad_ad versus the sweep
  coordinate, with the component blend fractions shown as a strip underneath and
  the ramp shaded on every panel;
* a **derivative-consistency figure** -- the returned partial (d_dlnd for a
  logRho sweep, d_dlnT for a logT sweep) against a central finite difference of
  the returned value, for lnPgas, lnE, lnS, Cv, grad_ad, plus their relative
  difference.

and two quantitative summaries:

* **(a)** the maximum relative mismatch between returned partials and the finite
  differences, inside the ramp versus outside it.  ``Do_Blend`` propagates the
  blend-fraction derivatives, so truncation-level agreement is expected
  everywhere (central-difference error ~ ``step^2``).
* **(b)** the value step at each edge: a smooth trend is fit on each side of the
  ramp (pure WDM_AV on one side, the pure lower-cascade component on the other)
  and extrapolated to the ramp center; the jump in lnE, lnS, lnPgas, Cv, grad_ad
  there is reported, for the mixture and for pure c12.

Internal-consistency only: no reference EOS is read, and radiation is left in
(both WDM_AV and Skye add it via pack_for_export).  The analysis therefore makes
no radiation correction.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

# make the builder package and this directory importable when run as a script
_HERE = Path(__file__).resolve().parent
for _p in (str(_HERE), str(_HERE.parent)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from mixtures import delta_s_ideal
from wdm_av_builder import constants as C

LN10 = C.ln10

# --------------------------------------------------------------------------
# paths
# --------------------------------------------------------------------------
DATA_ROOT = Path(os.environ.get("WDM_AV_MESA_DIR") or _HERE.parents[2]) / "data/eosWDM_AV_data"
BLEND_ROOT = DATA_ROOT / "blend"
DEFAULT_FIGDIR = DATA_ROOT.parents[1].parent / "claude_tex_notes/figures/blend"
DEFAULT_RESULTSDIR = _HERE / "results"

SETS = ("skye", "helm", "skye_ne", "helm_ne")
MIXTURES = ("c12_o16", "c12")
COMP = {"c12_o16": {"c12": 0.5, "o16": 0.5}, "c12": {"c12": 1.0}}
CUTS = (
    "rho_right_T6.0",
    "rho_right_T6.8",
    "T_top_R2.0",
    "T_top_R4.0",
    "rho_left_T6.8",
)

# column layout written by blend_sweep.f90 (0-based)
COLS = {
    "coord": 0,
    "logRho": 1,
    "logT": 2,
    "flag": 3,
    "lnPgas": 4,
    "lnE": 5,
    "lnS": 6,
    "chiRho": 7,
    "chiT": 8,
    "Cv": 9,
    "Cp": 10,
    "grad_ad": 11,
    "gamma1": 12,
    "dlnPgas_dlnRho": 13,
    "dlnPgas_dlnT": 14,
    "dlnE_dlnRho": 15,
    "dlnE_dlnT": 16,
    "dlnS_dlnRho": 17,
    "dlnS_dlnT": 18,
    "dCv_dlnRho": 19,
    "dCv_dlnT": 20,
    "dgrad_ad_dlnRho": 21,
    "dgrad_ad_dlnT": 22,
    "frac_HELM": 23,
    "frac_OPAL_SCVH": 24,
    "frac_FreeEOS": 25,
    "frac_PC": 26,
    "frac_Skye": 27,
    "frac_CMS": 28,
    "frac_ideal": 29,
    "frac_WDM_AV": 30,
}
NCOL = len(COLS)

FRAC_NAMES = (
    "HELM",
    "OPAL_SCVH",
    "FreeEOS",
    "PC",
    "Skye",
    "CMS",
    "ideal",
    "WDM_AV",
)

# (value, d_dlnd, d_dlnT) column names for the five derivative-checked quantities
DERIV_QUANTS = {
    "lnPgas": ("lnPgas", "dlnPgas_dlnRho", "dlnPgas_dlnT"),
    "lnE": ("lnE", "dlnE_dlnRho", "dlnE_dlnT"),
    "lnS": ("lnS", "dlnS_dlnRho", "dlnS_dlnT"),
    "Cv": ("Cv", "dCv_dlnRho", "dCv_dlnT"),
    "grad_ad": ("grad_ad", "dgrad_ad_dlnRho", "dgrad_ad_dlnT"),
}
VALUE_PANELS = ("lnPgas", "lnE", "lnS", "chiT", "Cv", "grad_ad")
STEP_QUANTS = ("lnPgas", "lnE", "lnS", "Cv", "grad_ad")

RAMP_LO, RAMP_HI = 1.0e-3, 1.0 - 1.0e-3  # frac_WDM_AV window defining the ramp


# --------------------------------------------------------------------------
# readers
# --------------------------------------------------------------------------
def read_cut(path: Path) -> tuple[np.ndarray, dict]:
    """Return (data[nrow, NCOL], header_dict) for a cut file."""
    header: dict = {}
    rows = []
    with open(path) as fh:
        for ln in fh:
            s = ln.rstrip("\n")
            st = s.strip()
            if not st or st.startswith("#"):
                continue
            tok = st.split()
            # data rows start with a signed number; header rows start with a key
            if st[0].isdigit() or st[0] in "+-.":
                rows.append([float(x) for x in tok[:NCOL]])
            else:
                header[tok[0]] = " ".join(tok[1:])
    return np.array(rows, dtype=float), header


# --------------------------------------------------------------------------
# geometry of a cut
# --------------------------------------------------------------------------
def sweep_axis(cut: str) -> str:
    """'logRho' if the cut varies density, else 'logT'."""
    return "logRho" if cut.startswith("rho_") else "logT"


def partial_name(cut: str, value_name: str) -> str:
    """The returned-partial column that matches the swept coordinate."""
    _, dd, dt = DERIV_QUANTS[value_name]
    return dd if sweep_axis(cut) == "logRho" else dt


def ramp_center(coord: np.ndarray, frac_wdm: np.ndarray) -> float | None:
    """Coordinate where frac_WDM_AV crosses 0.5 (linear interpolation)."""
    f = frac_wdm - 0.5
    sign = np.sign(f)
    idx = np.where(np.diff(sign) != 0)[0]
    if idx.size == 0:
        return None
    i = idx[0]
    x0, x1 = coord[i], coord[i + 1]
    y0, y1 = frac_wdm[i], frac_wdm[i + 1]
    if y1 == y0:
        return float(0.5 * (x0 + x1))
    return float(x0 + (0.5 - y0) * (x1 - x0) / (y1 - y0))


# --------------------------------------------------------------------------
# (a) finite-difference consistency
# --------------------------------------------------------------------------
def central_fd(value: np.ndarray, coord: np.ndarray) -> np.ndarray:
    """Central finite difference d(value)/d(ln x), x = 10**coord.

    Interior points only; endpoints are NaN.  coord is uniform (step dex), so
    d ln x = step * ln10.
    """
    fd = np.full_like(value, np.nan)
    step = coord[1] - coord[0]
    denom = 2.0 * step * LN10
    fd[1:-1] = (value[2:] - value[:-2]) / denom
    return fd


def fd_consistency(data: np.ndarray, cut: str) -> dict[str, dict]:
    """Returned partial vs central FD of the returned value, per quantity.

    The mismatch is normalised by a single robust scale per quantity (the max
    |returned partial| over the cut), not pointwise by |FD|: the latter blows up
    wherever the partial crosses zero and would mask the real question, which is
    whether the ramp degrades agreement.  Reports the max/median normalised
    mismatch inside the ramp versus outside, plus the arrays for plotting.
    """
    coord = data[:, COLS["coord"]]
    frac = data[:, COLS["frac_WDM_AV"]]
    in_ramp = (frac > RAMP_LO) & (frac < RAMP_HI)

    out: dict[str, dict] = {}
    for q in STEP_QUANTS:
        val = data[:, COLS[DERIV_QUANTS[q][0]]]
        ret = data[:, COLS[partial_name(cut, q)]]
        fd = central_fd(val, coord)
        interior = np.isfinite(fd) & np.isfinite(ret)
        # robust per-quantity scale: the span of the returned partial on the cut
        scale = np.nanmax(np.abs(ret[interior])) if interior.any() else 1.0
        scale = max(scale, 1.0e-300)
        rel = np.full_like(val, np.nan)
        rel[interior] = np.abs(ret[interior] - fd[interior]) / scale

        m_in = interior & in_ramp
        m_out = interior & ~in_ramp
        out[q] = {
            "coord": coord,
            "returned": ret,
            "fd": fd,
            "rel": rel,
            "scale": float(scale),
            "max_rel_in_ramp": _nanmax(rel[m_in]),
            "median_rel_in_ramp": _nanmedian(rel[m_in]),
            "max_rel_out_ramp": _nanmax(rel[m_out]),
            "median_rel_out_ramp": _nanmedian(rel[m_out]),
        }
    return out


def _nanmax(a: np.ndarray) -> float | None:
    a = a[np.isfinite(a)]
    return float(a.max()) if a.size else None


def _nanmedian(a: np.ndarray) -> float | None:
    a = a[np.isfinite(a)]
    return float(np.median(a)) if a.size else None


# --------------------------------------------------------------------------
# (b) value step across the ramp
# --------------------------------------------------------------------------
def _side_fit(coord, val, mask, center, nfit=40, deg=2):
    """Quadratic fit of val vs coord over the nfit points of mask nearest the
    ramp center; returns the value extrapolated to center (or None)."""
    idx = np.where(mask & np.isfinite(val))[0]
    if idx.size < deg + 1:
        return None
    order = idx[np.argsort(np.abs(coord[idx] - center))]
    use = order[: min(nfit, order.size)]
    deg_use = min(deg, use.size - 1)
    c = np.polyfit(coord[use], val[use], deg_use)
    return float(np.polyval(c, center))


def value_steps(data: np.ndarray, cut: str, mixture: str) -> dict:
    """Jump at the ramp center = (WDM_AV-side trend) - (lower-component trend).

    Both sides are extrapolated from their pure-component points.  The sign is
    (value WDM_AV approaches) minus (value the neighbour approaches).
    """
    coord = data[:, COLS["coord"]]
    frac = data[:, COLS["frac_WDM_AV"]]
    center = ramp_center(coord, frac)
    result: dict = {"ramp_center": center, "neighbour": None, "steps": {}}
    if center is None:
        return result

    wdm_mask = frac >= RAMP_HI
    nbr_mask = frac <= RAMP_LO
    if not wdm_mask.any() or not nbr_mask.any():
        return result

    # dominant neighbour component: evaluated at the pure-neighbour point nearest
    # the ramp center
    nbr_idx = np.where(nbr_mask)[0]
    j = nbr_idx[np.argmin(np.abs(coord[nbr_idx] - center))]
    fr = np.array([data[j, COLS[f"frac_{n}"]] for n in FRAC_NAMES[:-1]])
    result["neighbour"] = FRAC_NAMES[int(np.argmax(fr))]
    result["neighbour_frac"] = float(fr.max())

    for q in STEP_QUANTS:
        val = data[:, COLS[DERIV_QUANTS[q][0]]]
        w = _side_fit(coord, val, wdm_mask, center)
        n = _side_fit(coord, val, nbr_mask, center)
        if w is None or n is None:
            result["steps"][q] = None
            continue
        jump = w - n
        rel = jump / n if (q in ("Cv", "grad_ad") and n != 0.0) else jump
        # lnX quantities: the jump IS already d(lnX) ~ relative step in X
        result["steps"][q] = {
            "wdm_side": w,
            "neighbour_side": n,
            "jump": jump,
            "rel_or_dln": rel,
        }

    # expected lnS step from the ideal mixing entropy (Findings F1/F4): WDM_AV
    # counts it once, Skye twice, HELM zero times.  So against a Skye neighbour
    # d(lnS) ~ -Delta_s_ideal / S, against HELM ~ +Delta_s_ideal / S; ~0 for a
    # pure species.  Only strictly predictive for the matching table set.
    result["expected_lnS_step"] = _expected_lnS(result, COMP.get(mixture))
    return result


def _expected_lnS(result: dict, comp: dict | None) -> float | None:
    if comp is None:
        return None
    nbr = result.get("neighbour")
    s_step = result["steps"].get("lnS")
    if s_step is None:
        return None
    # S at the ramp centre from the neighbour-side extrapolation of lnS
    s_center = np.exp(s_step["neighbour_side"])
    ds = delta_s_ideal(comp)
    if nbr == "Skye":
        return -ds / s_center
    if nbr == "HELM":
        return +ds / s_center
    return None


# --------------------------------------------------------------------------
# plotting
# --------------------------------------------------------------------------
def _bounds_for_axis(header: dict, axis: str) -> list[float]:
    key = (
        "box_logRho_any_all_all_any" if axis == "logRho" else "box_logT_any_all_all_any"
    )
    return [float(x) for x in header[key].split()]


def plot_values(data, header, table_set, mixture, cut, figdir: Path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    coord = data[:, COLS["coord"]]
    frac = data[:, COLS["frac_WDM_AV"]]
    axis = sweep_axis(cut)
    bounds = [
        b for b in _bounds_for_axis(header, axis) if coord.min() <= b <= coord.max()
    ]
    in_ramp = (frac > RAMP_LO) & (frac < RAMP_HI)
    ramp_span = (coord[in_ramp].min(), coord[in_ramp].max()) if in_ramp.any() else None

    fig = plt.figure(figsize=(12, 11))
    gs = fig.add_gridspec(4, 2, height_ratios=[1, 1, 1, 0.7])
    axes = [fig.add_subplot(gs[r, c]) for r in range(3) for c in range(2)]
    axf = fig.add_subplot(gs[3, :])

    for ax, q in zip(axes, VALUE_PANELS):
        ax.plot(coord, data[:, COLS[q]], lw=1.2, color="C0")
        if ramp_span:
            ax.axvspan(*ramp_span, color="0.85", zorder=0)
        for b in bounds:
            ax.axvline(b, color="0.6", lw=0.6, ls=":")
        ax.set_title(q, fontsize=10)
        ax.set_xlabel(axis)
        ax.grid(alpha=0.25)

    # fractions strip: every component that is ever non-negligible
    for n in FRAC_NAMES:
        fr = data[:, COLS[f"frac_{n}"]]
        if np.nanmax(fr) > 1.0e-3:
            axf.plot(coord, fr, lw=1.3, label=n)
    if ramp_span:
        axf.axvspan(*ramp_span, color="0.85", zorder=0)
    axf.set_ylim(-0.05, 1.05)
    axf.set_ylabel("blend fraction")
    axf.set_xlabel(axis)
    axf.legend(ncol=4, fontsize=8, loc="center right")
    axf.grid(alpha=0.25)

    fig.suptitle(f"WDM_AV({table_set}) full cascade: {mixture}  cut {cut}", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    figdir.mkdir(parents=True, exist_ok=True)
    path = figdir / f"cut_{table_set}_{mixture}_{cut}.png"
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path


def plot_derivs(data, header, fdc, table_set, mixture, cut, figdir: Path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    coord = data[:, COLS["coord"]]
    frac = data[:, COLS["frac_WDM_AV"]]
    axis = sweep_axis(cut)
    dln = "dlnRho" if axis == "logRho" else "dlnT"
    in_ramp = (frac > RAMP_LO) & (frac < RAMP_HI)
    ramp_span = (coord[in_ramp].min(), coord[in_ramp].max()) if in_ramp.any() else None

    fig, axes = plt.subplots(len(STEP_QUANTS), 2, figsize=(12, 13))
    for row, q in enumerate(STEP_QUANTS):
        d = fdc[q]
        axL, axR = axes[row]
        axL.plot(coord, d["returned"], lw=1.3, color="C0", label="returned")
        axL.plot(coord, d["fd"], lw=1.0, color="C1", ls="--", label="central FD")
        axL.set_title(f"d({q})/{dln}", fontsize=9)
        axL.set_ylabel("partial")
        axL.legend(fontsize=8)
        with np.errstate(divide="ignore"):
            axR.semilogy(coord, d["rel"], lw=1.0, color="C3")
        axR.set_title(f"{q}: |returned - FD| / max|partial|", fontsize=9)
        axR.set_ylabel("rel. mismatch")
        for ax in (axL, axR):
            if ramp_span:
                ax.axvspan(*ramp_span, color="0.85", zorder=0)
            ax.set_xlabel(axis)
            ax.grid(alpha=0.25)

    fig.suptitle(
        f"WDM_AV({table_set}) partials vs central FD: {mixture}  cut {cut}",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    figdir.mkdir(parents=True, exist_ok=True)
    path = figdir / f"fd_{table_set}_{mixture}_{cut}.png"
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------
def run_one(table_set, mixture, cut, *, figdir: Path) -> dict:
    path = BLEND_ROOT / table_set / f"{mixture}_{cut}.txt"
    data, header = read_cut(path)
    fdc = fd_consistency(data, cut)
    steps = value_steps(data, cut, mixture)
    plot_values(data, header, table_set, mixture, cut, figdir)
    plot_derivs(data, header, fdc, table_set, mixture, cut, figdir)

    fd_summary = {
        q: {
            k: fdc[q][k]
            for k in (
                "max_rel_in_ramp",
                "median_rel_in_ramp",
                "max_rel_out_ramp",
                "median_rel_out_ramp",
            )
        }
        for q in STEP_QUANTS
    }
    return {
        "table_set": table_set,
        "mixture": mixture,
        "cut": cut,
        "n_nodes": int(data.shape[0]),
        "n_bad": int((data[:, COLS["flag"]] != 0).sum()),
        "fd_consistency": fd_summary,
        "value_steps": steps,
    }


def _fmt(x) -> str:
    if x is None:
        return "  -  "
    return f"{x:.2e}"


def markdown_report(summaries: list[dict]) -> str:
    lines = [
        "# WDM_AV blend-zone analysis (full eosDT cascade)",
        "",
        (
            "Internal consistency of `eosDT_get` across the WDM_AV blend box: "
            "returned partials vs central finite differences of the returned "
            "values, and the value step at each edge (pure-WDM_AV trend minus "
            "pure-neighbour trend, both extrapolated to the ramp centre)."
        ),
        "",
        "## (a) FD consistency: max relative mismatch (returned vs central FD)",
        "",
        "| set | mixture | cut | quantity | max in-ramp | max out-ramp |",
        "|---|---|---|---|---|---|",
    ]
    for s in summaries:
        for q in STEP_QUANTS:
            fc = s["fd_consistency"][q]
            lines.append(
                f"| {s['table_set']} | {s['mixture']} | {s['cut']} | {q} | "
                f"{_fmt(fc['max_rel_in_ramp'])} | {_fmt(fc['max_rel_out_ramp'])} |"
            )
    lines += [
        "",
        "## (b) value steps at the ramp centre",
        "",
        (
            "Jump = (WDM_AV-side trend) - (neighbour-side trend) at the ramp "
            "centre. For lnPgas/lnE/lnS the jump is d(lnX), i.e. the relative "
            "step in X; for Cv/grad_ad the relative step is jump/neighbour."
        ),
        "",
        "| set | mixture | cut | neighbour | lnPgas | lnE | lnS | lnS(exp) | Cv(rel) | grad_ad(rel) |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for s in summaries:
        st = s["value_steps"]
        steps = st.get("steps", {})

        def g(q, steps=steps):
            v = steps.get(q)
            return _fmt(v["rel_or_dln"]) if v else "  -  "

        lines.append(
            f"| {s['table_set']} | {s['mixture']} | {s['cut']} | "
            f"{st.get('neighbour')} | {g('lnPgas')} | {g('lnE')} | {g('lnS')} | "
            f"{_fmt(st.get('expected_lnS_step'))} | {g('Cv')} | {g('grad_ad')} |"
        )
    lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="python blend.py")
    p.add_argument("--sets", nargs="+", default=list(SETS), choices=list(SETS))
    p.add_argument("--mixtures", nargs="+", default=list(MIXTURES))
    p.add_argument("--cuts", nargs="+", default=list(CUTS))
    p.add_argument("--figdir", default=str(DEFAULT_FIGDIR))
    p.add_argument("--resultsdir", default=str(DEFAULT_RESULTSDIR))
    args = p.parse_args(argv)

    figdir = Path(args.figdir)
    summaries = []
    for table_set in args.sets:
        for mixture in args.mixtures:
            for cut in args.cuts:
                print(f"blend: {table_set}  {mixture}  {cut}")
                summaries.append(run_one(table_set, mixture, cut, figdir=figdir))

    resultsdir = Path(args.resultsdir)
    resultsdir.mkdir(parents=True, exist_ok=True)
    (resultsdir / "blend_summary.json").write_text(json.dumps(summaries, indent=2))
    (resultsdir / "blend_summary.md").write_text(markdown_report(summaries))
    print("wrote", resultsdir / "blend_summary.json", resultsdir / "blend_summary.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
