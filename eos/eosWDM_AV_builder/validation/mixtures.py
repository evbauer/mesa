"""Mixture-level validation of the integrated WDM_AV eosDT component.

Compares the additive-volume (AV) WDM_AV mixture, evaluated through the
integrated eosDT path (``eosDT_get_component(handle, i_eos_WDM_AV, ...)`` in
``test_fortran/mix_vs_reference.f90``), against MESA's Skye and HELM run
directly on the same five reference mixtures, over the (logRho, logT) plane.

Two independent physics questions (Decision D20, brief WP4b):

* **(a) AV-of-Skye-tables vs Skye-on-mixture** -- the differences are
  interpolation error plus the additive-volume approximation of Skye's Coulomb
  mixing and electron sharing.  Both sides include radiation (Skye ignores
  ``include_radiation``; WDM_AV's ``pack_for_export`` always adds it), so e, s,
  Cv, chiRho, chiT are directly comparable.  The Skye reference entropy is first
  corrected for the double-counted ideal mixing entropy (Finding F1): the
  reference ``S`` is high by ``Delta s_ideal = -(k_B/m_u) sum_i (X_i/A_i) ln y_i``
  (number fractions ``y_i``), which this module subtracts before comparing.  The
  uncorrected comparison is reported alongside as evidence.

* **(b) AV-of-HELM-tables vs HELM-on-mixture** -- HELM has no Coulomb term, so
  the differences are interpolation error plus the pure electron-fragmentation
  error of additive volume (which should vanish in the non-degenerate and fully
  degenerate limits and peak between).  The HELM reference is gas-only
  (``include_radiation = .false.``), so radiation is stripped analytically from
  WDM_AV's e, s, Cv; gas-only chiRho/chiT come from WDM_AV's own
  ``dlnPgas/dln{rho,T}`` columns.  The pure-species held-out interpolation level
  (from ``heldout.py``) is reported as the interpolation floor so the AV model
  error can be isolated.

Statistics (median, p90, p99, max of the relative error, absolute error for the
lnS partials) are reported over four regions:

* ``whole_valid``        -- every node where both sides are finite;
* ``skye_valid``         -- ``logT >= 6.5 OR logRho >= 3.7`` (Decision D20);
* ``skye_valid_no_pair`` -- the above minus the pair-radiation corner
  (``NOT (logT > 8 AND logRho < 0)``);
* ``target``             -- the WDM_AV coverage hole, ``logRho in [-2, 4]`` and
  ``logT in [3, 7]``.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

import numpy as np

# make the builder package and this directory importable when run as a script
_HERE = Path(__file__).resolve().parent
for _p in (str(_HERE), str(_HERE.parent)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from wdm_av_builder import constants as C

# --------------------------------------------------------------------------
# paths and constants
# --------------------------------------------------------------------------
DATA_ROOT = Path("/Users/bauer39/work/mesas/evbauer/mesa/data/eosWDM_AV_data")
REF_DIRS = {
    "skye": DATA_ROOT / "reference_skye_mixtures",
    "helm": DATA_ROOT / "reference_helm_mixtures",
}
OUT_DIRS = {
    "skye": DATA_ROOT / "validation_mixtures" / "skye",
    "helm": DATA_ROOT / "validation_mixtures" / "helm",
}
SOURCE_DIR = DATA_ROOT / "source_skye"  # atomic weights (identical for helm)

DEFAULT_FIGDIR = Path(
    "/Users/bauer39/work/mesas/evbauer/claude_tex_notes/figures/mixtures"
)
DEFAULT_RESULTSDIR = _HERE / "results"
DEFAULT_HELDOUT = DEFAULT_RESULTSDIR / "heldout_summary.json"

MIXTURES = (
    "mix_c12_o16",
    "mix_o16_ne20_mg24",
    "mix_he4_c12",
    "mix_h1_he4",
    "mix_c12_o16_ne22",
)

KBM = C.kerg / C.amu  # k_B / m_u  (erg g^-1 K^-1)

# the five relative-error quantities drawn on the maps, in panel order
MAP_QUANTITIES = (
    "P",
    "E",
    "S",
    "chiRho",
    "chiT",
    "Cv",
    "grad_ad",
    "gamma1",
)
# quantities carried in the summary (S is the F1-corrected entropy for skye)
REL_QUANTITIES = MAP_QUANTITIES
ABS_QUANTITIES = ("dlnS_dlnRho", "dlnS_dlnT")

REGION_NAMES = ("whole_valid", "skye_valid", "skye_valid_no_pair", "target")

# map our quantity names to the held-out study's names (interpolation floor)
HELDOUT_QNAME = {
    "P": "P",
    "E": "e",
    "S": "s",
    "Cv": "c_V",
    "chiRho": "chi_rho",
    "chiT": "chi_T",
}

_REF_NAMES = [
    "logRho",
    "logT",
    "lnPgas",
    "lnE",
    "lnS",
    "chiRho",
    "chiT",
    "Cv",
    "Cp",
    "grad_ad",
    "gamma1",
    "gamma3",
    "eta",
    "mu",
    "lnfree_e",
    "dlnPgas_dlnRho",
    "dlnPgas_dlnT",
    "dlnE_dlnRho",
    "dlnE_dlnT",
    "dlnS_dlnRho",
    "dlnS_dlnT",
    "dCv_dlnRho",
    "dCv_dlnT",
    "dgrad_ad_dlnRho",
    "dgrad_ad_dlnT",
]
_WDM_NAMES = [
    "logRho",
    "logT",
    "flag",
    "lnPgas",
    "lnE",
    "lnS",
    "chiRho",
    "chiT",
    "Cv",
    "Cp",
    "grad_ad",
    "gamma1",
    "gamma3",
    "eta",
    "mu",
    "lnfree_e",
    "dlnPgas_dlnRho",
    "dlnPgas_dlnT",
    "dlnE_dlnRho",
    "dlnE_dlnT",
    "dlnS_dlnRho",
    "dlnS_dlnT",
    "dCv_dlnRho",
    "dCv_dlnT",
    "dgrad_ad_dlnRho",
    "dgrad_ad_dlnT",
]
REF_COLS = {n: i for i, n in enumerate(_REF_NAMES)}
WDM_COLS = {n: i for i, n in enumerate(_WDM_NAMES)}


# --------------------------------------------------------------------------
# readers
# --------------------------------------------------------------------------
def _is_data_line(s: str) -> bool:
    return bool(s) and s[0] in "0123456789+-."


def read_table(path: Path, ncols: int) -> np.ndarray:
    """Read a whitespace table, keeping exactly the first ``ncols`` columns.

    The reference and output writers emit one extra column on declined (all-NaN)
    rows (an off-by-one in the Fortran NaN fill), so rows are ragged; trimming to
    ``ncols`` makes the array rectangular.  'NaN' tokens parse to ``np.nan``.
    """
    rows = []
    with open(path) as fh:
        for ln in fh:
            s = ln.strip()
            if not _is_data_line(s):
                continue
            toks = s.split()
            if len(toks) < ncols:
                toks = toks + ["NaN"] * (ncols - len(toks))
            rows.append(toks[:ncols])
    return np.array(rows, dtype=float)


def read_header(path: Path) -> dict:
    """Parse composition and grid definition from a reference/output header."""
    comp_str = ""
    grid = {}
    with open(path) as fh:
        for ln in fh:
            s = ln.strip()
            if not s:
                continue
            if _is_data_line(s):
                break
            if s.startswith("#"):
                continue
            key = s.split()[0]
            if key == "composition":
                m = re.search(r'"([^"]*)"', s)
                if m:
                    comp_str = m.group(1)
            elif key in ("nlogRho", "nlogT"):
                p = s.split()
                grid[key] = int(p[1])
                grid[key + "_min"] = float(p[3])
                grid[key + "_max"] = float(p[5])
    return {"composition": comp_str, **grid}


def parse_composition(comp_str: str) -> dict[str, float]:
    """'c12=0.50 o16=0.50' -> {'c12': 0.5, 'o16': 0.5} (normalised)."""
    comp = {}
    for tok in comp_str.split():
        if "=" not in tok:
            continue
        name, frac = tok.split("=", 1)
        comp[name.strip()] = float(frac)
    tot = sum(comp.values())
    if tot > 0:
        comp = {k: v / tot for k, v in comp.items()}
    return comp


_A_CACHE: dict[str, float] = {}


def atomic_weight(name: str) -> float:
    """Atomic weight A for an isotope, read from the source-table header."""
    if name not in _A_CACHE:
        src = SOURCE_DIR / f"eosWDM_AV_{name}.txt"
        a = None
        with open(src) as fh:
            for ln in fh:
                if ln.startswith("A "):
                    a = float(ln.split()[1])
                    break
        if a is None:
            raise ValueError(f"no A header in {src}")
        _A_CACHE[name] = a
    return _A_CACHE[name]


# --------------------------------------------------------------------------
# physics helpers
# --------------------------------------------------------------------------
def delta_s_ideal(comp: dict[str, float]) -> float:
    """Ideal entropy of mixing per gram (erg g^-1 K^-1), Finding F1.

    ``Delta s_ideal = -(k_B/m_u) sum_i (X_i/A_i) ln y_i`` with number fractions
    ``y_i``.  Skye's mixture S is high by this amount; subtract it from the Skye
    reference before comparing to WDM_AV (which counts mixing once).
    """
    n = {k: comp[k] / atomic_weight(k) for k in comp}
    ntot = sum(n.values())
    if ntot <= 0:
        return 0.0
    y = {k: n[k] / ntot for k in comp}
    return -KBM * sum(n[k] * np.log(y[k]) for k in comp if y[k] > 0)


def _col(arr: np.ndarray, cols: dict, name: str) -> np.ndarray:
    return arr[:, cols[name]]


def comparison_quantities(
    table_set: str,
    ref: np.ndarray,
    wdm: np.ndarray,
    comp: dict[str, float],
    logrho: np.ndarray,
    logt: np.ndarray,
) -> dict:
    """Build (wdm, ref) arrays for each compared quantity for one mixture/set.

    Returns dict ``name -> {"w": ..., "r": ..., "kind": "rel"|"abs"}`` plus the
    auxiliary ``_extra`` dict (S_uncorrected, radiation-clean mask, delta_s).
    """
    rho = np.power(10.0, logrho)
    T = np.power(10.0, logt)

    def rc(name):  # reference column
        return _col(ref, REF_COLS, name)

    def wc(name):  # wdm output column
        return _col(wdm, WDM_COLS, name)

    out: dict = {}
    extra: dict = {}

    # gas pressure: gas-only on both sides, directly comparable
    out["P"] = {"w": np.exp(wc("lnPgas")), "r": np.exp(rc("lnPgas")), "kind": "rel"}

    if table_set == "skye":
        ds = delta_s_ideal(comp)
        extra["delta_s_ideal"] = ds
        S_ref = np.exp(rc("lnS"))
        S_ref_corr = S_ref - ds
        S_wdm = np.exp(wc("lnS"))
        out["E"] = {"w": np.exp(wc("lnE")), "r": np.exp(rc("lnE")), "kind": "rel"}
        out["S"] = {"w": S_wdm, "r": S_ref_corr, "kind": "rel"}
        # evidence: entropy WITHOUT the F1 correction
        extra["S_uncorrected"] = {"w": S_wdm, "r": S_ref, "kind": "rel"}
        out["chiRho"] = {"w": wc("chiRho"), "r": rc("chiRho"), "kind": "rel"}
        out["chiT"] = {"w": wc("chiT"), "r": rc("chiT"), "kind": "rel"}
        out["Cv"] = {"w": wc("Cv"), "r": rc("Cv"), "kind": "rel"}
        out["grad_ad"] = {"w": wc("grad_ad"), "r": rc("grad_ad"), "kind": "rel"}
        out["gamma1"] = {"w": wc("gamma1"), "r": rc("gamma1"), "kind": "rel"}
        # lnS partials: dS/dlnX is unchanged by the constant F1 shift, so the
        # corrected lnS partial is ref_partial * S_ref / S_ref_corr.
        scale = S_ref / S_ref_corr
        out["dlnS_dlnRho"] = {
            "w": wc("dlnS_dlnRho"),
            "r": rc("dlnS_dlnRho") * scale,
            "kind": "abs",
        }
        out["dlnS_dlnT"] = {
            "w": wc("dlnS_dlnT"),
            "r": rc("dlnS_dlnT") * scale,
            "kind": "abs",
        }
    else:  # helm: strip radiation from WDM_AV (reference is gas-only)
        e_rad = C.crad * np.power(T, 4) / rho
        s_rad = 4.0 * C.crad * np.power(T, 3) / (3.0 * rho)
        cv_rad = 4.0 * C.crad * np.power(T, 3) / rho
        S_wdm_tot = np.exp(wc("lnS"))
        S_wdm_gas = S_wdm_tot - s_rad
        out["E"] = {
            "w": np.exp(wc("lnE")) - e_rad,
            "r": np.exp(rc("lnE")),
            "kind": "rel",
        }
        out["S"] = {"w": S_wdm_gas, "r": np.exp(rc("lnS")), "kind": "rel"}
        out["Cv"] = {"w": wc("Cv") - cv_rad, "r": rc("Cv"), "kind": "rel"}
        # gas-only chiRho/chiT: WDM_AV's chi refer to TOTAL pressure (like Skye),
        # but its dlnPgas/dln{rho,T} columns are exactly the gas-only chi.
        out["chiRho"] = {"w": wc("dlnPgas_dlnRho"), "r": rc("chiRho"), "kind": "rel"}
        out["chiT"] = {"w": wc("dlnPgas_dlnT"), "r": rc("chiT"), "kind": "rel"}
        # grad_ad / gamma1: WDM_AV includes radiation, the HELM reference does
        # not.  Kept for completeness but only meaningful where radiation is
        # negligible; a radiation-clean mask (Prad/Pgas < 1e-2) is recorded.
        out["grad_ad"] = {"w": wc("grad_ad"), "r": rc("grad_ad"), "kind": "rel"}
        out["gamma1"] = {"w": wc("gamma1"), "r": rc("gamma1"), "kind": "rel"}
        Pgas = np.exp(rc("lnPgas"))
        Prad = C.crad * np.power(T, 4) / 3.0
        extra["rad_clean"] = (Prad / Pgas) < 1.0e-2
        # gas-only lnS partials of WDM_AV, built from the total-S partials:
        #   dS_gas/dlnRho = dS_tot/dlnRho + s_rad   (s_rad ~ 1/rho)
        #   dS_gas/dlnT   = dS_tot/dlnT   - 3 s_rad (s_rad ~ T^3)
        dStot_dlnRho = wc("dlnS_dlnRho") * S_wdm_tot
        dStot_dlnT = wc("dlnS_dlnT") * S_wdm_tot
        with np.errstate(divide="ignore", invalid="ignore"):
            w_dlnS_dlnRho = (dStot_dlnRho + s_rad) / S_wdm_gas
            w_dlnS_dlnT = (dStot_dlnT - 3.0 * s_rad) / S_wdm_gas
        out["dlnS_dlnRho"] = {
            "w": w_dlnS_dlnRho,
            "r": rc("dlnS_dlnRho"),
            "kind": "abs",
        }
        out["dlnS_dlnT"] = {"w": w_dlnS_dlnT, "r": rc("dlnS_dlnT"), "kind": "abs"}

    return {"quantities": out, "extra": extra}


# --------------------------------------------------------------------------
# regions and statistics
# --------------------------------------------------------------------------
def regions(logrho: np.ndarray, logt: np.ndarray) -> dict[str, np.ndarray]:
    skye = (logt >= 6.5) | (logrho >= 3.7)
    pair = (logt > 8.0) & (logrho < 0.0)
    return {
        "whole_valid": np.ones_like(skye, dtype=bool),
        "skye_valid": skye,
        "skye_valid_no_pair": skye & ~pair,
        "target": (logrho >= -2.0) & (logrho <= 4.0) & (logt >= 3.0) & (logt <= 7.0),
    }


def error_array(w: np.ndarray, r: np.ndarray, kind: str, good: np.ndarray):
    """Elementwise error (|w-r|/|r| for rel, |w-r| for abs) and its valid mask."""
    err = np.full(w.shape, np.nan)
    if kind == "rel":
        ok = good & np.isfinite(w) & np.isfinite(r) & (r != 0.0)
        err[ok] = np.abs(w[ok] - r[ok]) / np.abs(r[ok])
    else:
        ok = good & np.isfinite(w) & np.isfinite(r)
        err[ok] = np.abs(w[ok] - r[ok])
    return err, ok


def stats(err: np.ndarray, mask: np.ndarray) -> dict:
    a = err[mask]
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


# --------------------------------------------------------------------------
# held-out interpolation floor
# --------------------------------------------------------------------------
def load_heldout_floor(path: Path) -> dict:
    if not path.exists():
        return {}
    with open(path) as fh:
        data = json.load(fh)
    floor = {}  # (tag, isotope) -> {region -> {hq -> {median,p99}}}
    for s in data:
        iso = s["isotope"]
        tag = s["label"].split("_")[0]  # skye / helm
        qd = s["variants"].get("a_o5", {}).get("quantities", {})
        floor[(tag, iso)] = qd
    return floor


def mixture_floor(
    floor: dict, table_set: str, comp: dict[str, float], region: str
) -> dict:
    """Mass-weighted pure-species interpolation floor for the mixture.

    Uses the matching table set's held-out numbers where available, else the
    skye held-out as a proxy (noted in the markdown).
    """
    res = {}
    for our_q, hq in HELDOUT_QNAME.items():
        acc, wsum, proxy = 0.0, 0.0, False
        for name, x in comp.items():
            qd = floor.get((table_set, name))
            if qd is None:
                qd = floor.get(("skye", name))
                proxy = True
            if qd is None or hq not in qd:
                continue
            med = qd[hq][region]["median"]
            if med is None:
                continue
            acc += x * med
            wsum += x
        if wsum > 0:
            res[our_q] = {"median": acc / wsum, "proxy": proxy}
    return res


# --------------------------------------------------------------------------
# per mixture driver
# --------------------------------------------------------------------------
def run_one(table_set: str, mixture: str, *, figdir: Path, floor: dict) -> dict:
    ref_path = REF_DIRS[table_set] / f"{mixture}.txt"
    out_path = OUT_DIRS[table_set] / f"{mixture}.txt"
    hdr = read_header(out_path)
    comp = parse_composition(hdr["composition"])

    ref = read_table(ref_path, len(_REF_NAMES))
    wdm = read_table(out_path, len(_WDM_NAMES))
    if ref.shape[0] != wdm.shape[0]:
        raise ValueError(f"{mixture}/{table_set}: row count mismatch")

    logrho = _col(wdm, WDM_COLS, "logRho")
    logt = _col(wdm, WDM_COLS, "logT")
    flag = _col(wdm, WDM_COLS, "flag")
    if not np.allclose(logrho, _col(ref, REF_COLS, "logRho")):
        raise ValueError(f"{mixture}/{table_set}: grid mismatch")

    wdm_valid = flag == 0.0
    ref_valid = np.isfinite(_col(ref, REF_COLS, "lnPgas"))
    both = wdm_valid & ref_valid

    cq = comparison_quantities(table_set, ref, wdm, comp, logrho, logt)
    reg = regions(logrho, logt)

    summary = {
        "table_set": table_set,
        "mixture": mixture,
        "composition": comp,
        "delta_s_ideal": cq["extra"].get("delta_s_ideal"),
        "n_nodes": int(logrho.size),
        "n_wdm_valid": int(wdm_valid.sum()),
        "n_ref_valid": int(ref_valid.sum()),
        "n_both_valid": int(both.sum()),
        "n_ref_valid_wdm_declined": int((ref_valid & ~wdm_valid).sum()),
        "coverage": {},
        "quantities": {},
        "S_uncorrected": {},
        "floor": {},
    }

    # coverage per region
    for rn in REGION_NAMES:
        m = reg[rn]
        summary["coverage"][rn] = {
            "n_region": int(m.sum()),
            "n_ref_valid": int((m & ref_valid).sum()),
            "n_both_valid": int((m & both).sum()),
        }

    maps: dict[str, np.ndarray] = {}
    for q in REL_QUANTITIES + ABS_QUANTITIES:
        spec = cq["quantities"][q]
        err, ok = error_array(spec["w"], spec["r"], spec["kind"], both)
        if q in MAP_QUANTITIES:
            maps[q] = err
        qstat = {"kind": spec["kind"]}
        for rn in REGION_NAMES:
            m = reg[rn] & ok
            if table_set == "helm" and q in ("grad_ad", "gamma1"):
                m = m & cq["extra"]["rad_clean"]
            qstat[rn] = stats(err, m)
        summary["quantities"][q] = qstat

    # uncorrected-entropy evidence (skye only)
    if "S_uncorrected" in cq["extra"]:
        spec = cq["extra"]["S_uncorrected"]
        err, ok = error_array(spec["w"], spec["r"], spec["kind"], both)
        for rn in REGION_NAMES:
            summary["S_uncorrected"][rn] = stats(err, reg[rn] & ok)

    # interpolation floor
    for rn in ("skye_valid", "target"):
        summary["floor"][rn] = mixture_floor(floor, table_set, comp, rn)

    plot_maps(table_set, mixture, logrho, logt, maps, figdir)
    return summary


# --------------------------------------------------------------------------
# plotting
# --------------------------------------------------------------------------
def plot_maps(table_set, mixture, logrho, logt, maps, figdir: Path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # infer the grid shape (logT outer, logRho inner)
    nR = int(np.sum(logt == logt[0]))
    nT = logt.size // nR
    LR = logrho.reshape(nT, nR)
    LT = logt.reshape(nT, nR)

    figdir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 4, figsize=(17, 8))
    axes = axes.ravel()
    for ax, q in zip(axes, MAP_QUANTITIES):
        z = maps[q].reshape(nT, nR)
        with np.errstate(divide="ignore"):
            lz = np.log10(z)
        finite = np.isfinite(lz)
        if finite.any():
            vmax = np.nanpercentile(lz[finite], 99.5)
            vmin = max(np.nanmin(lz[finite]), vmax - 8.0)
        else:
            vmin, vmax = -16.0, 0.0
        im = ax.pcolormesh(
            LR, LT, lz, shading="auto", cmap="viridis", vmin=vmin, vmax=vmax
        )
        ax.axhline(6.5, color="w", lw=0.6, ls="--", alpha=0.6)
        ax.axvline(3.7, color="w", lw=0.6, ls="--", alpha=0.6)
        # target box
        ax.plot([-2, 4, 4, -2, -2], [3, 3, 7, 7, 3], color="r", lw=0.7, alpha=0.7)
        ax.set_title(f"{q}  log10|rel err|", fontsize=9)
        ax.set_xlabel("logRho")
        ax.set_ylabel("logT")
        plt.colorbar(im, ax=ax)
    fig.suptitle(f"WDM_AV({table_set}) vs {table_set}-on-mixture: {mixture}")
    fig.tight_layout()
    path = figdir / f"{table_set}_{mixture}.png"
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------
def write_reports(summaries: list[dict], resultsdir: Path) -> dict[str, Path]:
    resultsdir = Path(resultsdir)
    resultsdir.mkdir(parents=True, exist_ok=True)

    json_path = resultsdir / "mixtures_summary.json"
    with json_path.open("w") as fh:
        json.dump(summaries, fh, indent=2)

    csv_path = resultsdir / "mixtures_summary.csv"
    fields = [
        "table_set",
        "mixture",
        "quantity",
        "kind",
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
            for q, qstat in s["quantities"].items():
                for region in REGION_NAMES:
                    st = qstat[region]
                    w.writerow(
                        {
                            "table_set": s["table_set"],
                            "mixture": s["mixture"],
                            "quantity": q,
                            "kind": qstat["kind"],
                            "region": region,
                            **st,
                        }
                    )
            for region in REGION_NAMES:
                st = s["S_uncorrected"].get(region)
                if st:
                    w.writerow(
                        {
                            "table_set": s["table_set"],
                            "mixture": s["mixture"],
                            "quantity": "S_uncorrected",
                            "kind": "rel",
                            "region": region,
                            **st,
                        }
                    )

    md_path = resultsdir / "mixtures_summary.md"
    md_path.write_text(markdown_report(summaries))
    return {"json": json_path, "csv": csv_path, "md": md_path}


def _fmt(x) -> str:
    if x is None:
        return "  -  "
    return f"{x:.1e}"


_REGION_DESC = {
    "whole_valid": "whole common-valid region",
    "skye_valid": "Skye-valid region (logT >= 6.5 OR logRho >= 3.7)",
    "skye_valid_no_pair": "Skye-valid minus pair corner (NOT (logT>8 AND logRho<0))",
    "target": "target region (logRho in [-2,4], logT in [3,7])",
}


def _agg_median(summaries, table_set, q, region) -> float | None:
    vals = [
        s["quantities"][q][region]["median"]
        for s in summaries
        if s["table_set"] == table_set
        and s["quantities"][q][region]["median"] is not None
    ]
    return float(np.median(vals)) if vals else None


def markdown_report(summaries: list[dict]) -> str:
    qs = ("P", "E", "S", "chiT", "Cv", "grad_ad")
    lines = [
        "# Mixture-level WDM_AV validation vs Skye / HELM",
        "",
        "Relative error of the additive-volume WDM_AV mixture (integrated eosDT "
        "path, `eosDT_get_component(i_eos_WDM_AV)`) against Skye / HELM run "
        "directly on the same five reference mixtures, at the reference grid "
        "nodes where both sides are finite. `S` is the F1-corrected entropy for "
        "Skye (reference minus `Delta s_ideal`); for HELM, radiation is stripped "
        "from WDM_AV's e/s/Cv and gas-only chiRho/chiT come from WDM_AV's "
        "dlnPgas columns. HELM grad_ad/gamma1 mix a radiation-definition "
        "difference and are taken only where Prad/Pgas < 1e-2. Full per-mixture, "
        "per-region data (p90, p99, max, all quantities) is in the CSV/JSON.",
        "",
    ]

    # (A) aggregate median table, per set and region
    lines.append("## Aggregate median relative error (median over the 5 mixtures)")
    lines.append("")
    for table_set, tag in (("skye", "(a) Skye"), ("helm", "(b) HELM")):
        present = [s for s in summaries if s["table_set"] == table_set]
        if not present:
            continue
        lines.append(f"### {tag}")
        lines.append("")
        lines.append("| region | " + " | ".join(qs) + " |")
        lines.append("|" + "---|" * (1 + len(qs)))
        for region in REGION_NAMES:
            row = [region]
            for q in qs:
                row.append(_fmt(_agg_median(summaries, table_set, q, region)))
            lines.append("| " + " | ".join(row) + " |")
        lines.append("")

    # (B) per-mixture detail (median / p99 / max), target + skye_valid
    for region in ("target", "skye_valid"):
        lines.append(f"## Per-mixture, {_REGION_DESC[region]}")
        lines.append("")
        lines.append("| set | mixture | stat | " + " | ".join(qs) + " |")
        lines.append("|" + "---|" * (3 + len(qs)))
        for s in summaries:
            for stat in ("median", "p99", "max"):
                row = [s["table_set"], s["mixture"], stat]
                for q in qs:
                    row.append(_fmt(s["quantities"][q][region][stat]))
                lines.append("| " + " | ".join(row) + " |")
        lines.append("")

    # (C) F1 evidence: S uncorrected vs corrected (Skye)
    lines.append("## Entropy: F1 correction evidence (Skye)")
    lines.append("")
    lines.append(
        "Median relative error of S without and with the F1 correction "
        "(subtracting `Delta s_ideal` from the Skye reference). The correction "
        "removes the systematic entropy offset from Skye's double-counted ideal "
        "mixing entropy."
    )
    lines.append("")
    lines.append(
        "| mixture | Delta_s_ideal [erg/g/K] | S uncorr (skye_valid) | "
        "S corr (skye_valid) | S uncorr (target) | S corr (target) |"
    )
    lines.append("|" + "---|" * 6)
    for s in summaries:
        if s["table_set"] != "skye":
            continue
        unc_sv = s["S_uncorrected"].get("skye_valid", {}).get("median")
        cor_sv = s["quantities"]["S"]["skye_valid"]["median"]
        unc_tg = s["S_uncorrected"].get("target", {}).get("median")
        cor_tg = s["quantities"]["S"]["target"]["median"]
        lines.append(
            f"| {s['mixture']} | {s['delta_s_ideal']:.3e} | {_fmt(unc_sv)} | "
            f"{_fmt(cor_sv)} | {_fmt(unc_tg)} | {_fmt(cor_tg)} |"
        )
    lines.append("")

    # (D) HELM interpolation floor vs observed AV error
    lines.append("## HELM: interpolation floor vs observed error")
    lines.append("")
    lines.append(
        "Pure-species held-out interpolation median (mass-weighted over the "
        "mixture, a_o5 variant) vs the observed AV mixture median, so the AV "
        "model error can be isolated. A `*` marks a floor taken from the Skye "
        "held-out as a proxy (no HELM held-out table for that isotope)."
    )
    lines.append("")
    fq = ("P", "E", "S", "Cv", "chiRho", "chiT")
    for region in ("target", "skye_valid"):
        lines.append(f"### {_REGION_DESC[region]}")
        lines.append("")
        lines.append("| mixture | source | " + " | ".join(fq) + " |")
        lines.append("|" + "---|" * (2 + len(fq)))
        for s in summaries:
            if s["table_set"] != "helm":
                continue
            fl = s["floor"].get(region, {})
            frow = [s["mixture"], "interp floor"]
            for q in fq:
                v = fl.get(q)
                if v is None:
                    frow.append("  -  ")
                else:
                    frow.append(_fmt(v["median"]) + ("*" if v["proxy"] else ""))
            lines.append("| " + " | ".join(frow) + " |")
            orow = ["", "observed AV"]
            for q in fq:
                orow.append(_fmt(s["quantities"][q][region]["median"]))
            lines.append("| " + " | ".join(orow) + " |")
        lines.append("")

    # (E) coverage / declines
    lines.append("## Coverage (integrated-path declines)")
    lines.append("")
    lines.append(
        "`eosDT_get_component(i_eos_WDM_AV)` calls the component directly (no "
        "blend box); a node is declined when the AV solver fails (no-converge / "
        "bad-domain / out-of-table), in which case WDM_AV returns the all-zero "
        "sentinel that `mix_vs_reference` records as flag = -1."
    )
    lines.append("")
    lines.append(
        "| set | mixture | both-valid / total | both-valid in target | "
        "ref-valid but WDM declined |"
    )
    lines.append("|" + "---|" * 5)
    for s in summaries:
        cov_t = s["coverage"]["target"]
        lines.append(
            f"| {s['table_set']} | {s['mixture']} | "
            f"{s['n_both_valid']}/{s['n_nodes']} | "
            f"{cov_t['n_both_valid']}/{cov_t['n_region']} | "
            f"{s['n_ref_valid_wdm_declined']} |"
        )
    lines.append("")
    return "\n".join(lines)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python mixtures.py",
        description="mixture-level WDM_AV validation vs Skye/HELM",
    )
    p.add_argument(
        "--sets", nargs="+", default=["skye", "helm"], choices=["skye", "helm"]
    )
    p.add_argument("--mixtures", nargs="+", default=list(MIXTURES))
    p.add_argument("--figdir", default=str(DEFAULT_FIGDIR))
    p.add_argument("--resultsdir", default=str(DEFAULT_RESULTSDIR))
    p.add_argument("--heldout", default=str(DEFAULT_HELDOUT))
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    floor = load_heldout_floor(Path(args.heldout))
    figdir = Path(args.figdir)
    summaries = []
    for table_set in args.sets:
        for mixture in args.mixtures:
            print(f"mixtures: {table_set}  {mixture}")
            summary = run_one(table_set, mixture, figdir=figdir, floor=floor)
            summaries.append(summary)
            tgt = summary["quantities"]
            msg = "  ".join(
                f"{q}={_fmt(tgt[q]['target']['median'])}"
                for q in ("P", "E", "S", "chiT", "Cv")
            )
            print(
                f"    target median: {msg}   "
                f"(both-valid {summary['coverage']['target']['n_both_valid']}/"
                f"{summary['coverage']['target']['n_region']})"
            )
    paths = write_reports(summaries, Path(args.resultsdir))
    print("wrote", paths["json"], paths["csv"], paths["md"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
