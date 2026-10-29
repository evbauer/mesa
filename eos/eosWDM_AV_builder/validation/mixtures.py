"""Mixture-level validation of the integrated WDM_AV eosDT component.

Compares the additive-volume (AV) WDM_AV mixture, evaluated through the
integrated eosDT path (``eosDT_get_component(handle, i_eos_WDM_AV, ...)`` in
``test_fortran/mix_vs_reference.f90``), against MESA's Skye and HELM run
directly on the same eleven reference mixtures, over the (logRho, logT) plane.

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
    "skye_fd": DATA_ROOT
    / "reference_skye_mixtures",  # P,E,F-only tables, FD 2nd derivs
    # common-electron-density mixing rule (Decision D27); same references, the
    # WDM_AV output comes from the '<set>_ne' directories (mix_vs_reference arg2).
    "skye_ne": DATA_ROOT / "reference_skye_mixtures",
    "helm_ne": DATA_ROOT / "reference_helm_mixtures",
}
OUT_DIRS = {
    "skye": DATA_ROOT / "validation_mixtures" / "skye",
    "helm": DATA_ROOT / "validation_mixtures" / "helm",
    "skye_fd": DATA_ROOT / "validation_mixtures" / "skye_fd",
    "skye_ne": DATA_ROOT / "validation_mixtures" / "skye_ne",
    "helm_ne": DATA_ROOT / "validation_mixtures" / "helm_ne",
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
    # WP4c: asymmetric binaries and ternaries (compositions in the file headers)
    "mix_c12_o16_90_10",
    "mix_c12_o16_99_01",
    "mix_he4_c12_99_01",
    "mix_he4_c12_10_90",
    "mix_c12_o16_ne20",
    "mix_he4_c12_o16",
)

KBM = C.kerg / C.amu  # k_B / m_u  (erg g^-1 K^-1)

# eV -> erg (CODATA 2018, matches MESA const_def ev2erg)
EV2ERG = 1.602176634e-12

# Total electron binding (full-ionization) energies I_Z in eV, indexed by nuclear
# charge Z, exactly as in eos/private/ion_offset.f90 (compute_ion_offset, NIST
# H-Ni total binding energies).  Only the first 12 are needed for the reference
# mixtures (H..Mg).
ION_EV = {
    1: 13.598434599702,
    2: 79.0051545,
    3: 203.4861711,
    4: 399.14864,
    5: 670.9809,
    6: 1030.1085,
    7: 1486.058,
    8: 2043.8428,
    9: 2715.89,
    10: 3511.696,
    11: 4420.0,
    12: 5451.06,
}

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
_Z_CACHE: dict[str, int] = {}


def _source_header_value(name: str, key: str) -> str:
    src = SOURCE_DIR / f"eosWDM_AV_{name}.txt"
    with open(src) as fh:
        for ln in fh:
            toks = ln.split()
            if len(toks) >= 2 and toks[0] == key:
                return toks[1]
    raise ValueError(f"no {key} header in {src}")


def atomic_weight(name: str) -> float:
    """Atomic weight A for an isotope, read from the source-table header."""
    if name not in _A_CACHE:
        _A_CACHE[name] = float(_source_header_value(name, "A"))
    return _A_CACHE[name]


def charge(name: str) -> int:
    """Nuclear charge Z for an isotope, read from the source-table header."""
    if name not in _Z_CACHE:
        _Z_CACHE[name] = int(_source_header_value(name, "Z"))
    return _Z_CACHE[name]


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


def ion_offset_dE(comp: dict[str, float]) -> float:
    """Skye ion-offset non-additivity correction for energy (erg/g).

    Skye's ``compute_ion_offset`` (eos/private/ion_offset.f90) adds the mixture
    ion offset ``sum_k y_k I_{Z_k} / m_u`` (number fractions ``y_k``), whereas
    WDM_AV builds the mixture from per-species tables and so carries the
    mass-weighted pure-species offset ``sum_k X_k I_{Z_k} / m_u``.  The two differ
    by a composition-only constant.  Adding this to the Skye reference E puts it on
    the WDM_AV (mass-weighted) convention:

        dE = (sum_k X_k I_k - sum_k y_k I_k) * ev2erg / m_u .
    """
    n = {k: comp[k] / atomic_weight(k) for k in comp}
    ntot = sum(n.values())
    if ntot <= 0:
        return 0.0
    y = {k: n[k] / ntot for k in comp}
    sum_X_I = sum(comp[k] * ION_EV[charge(k)] for k in comp)
    sum_y_I = sum(y[k] * ION_EV[charge(k)] for k in comp)
    return (sum_X_I - sum_y_I) * EV2ERG / C.amu


def _col(arr: np.ndarray, cols: dict, name: str) -> np.ndarray:
    return arr[:, cols[name]]


def comparison_quantities(
    table_set: str,
    ref: np.ndarray,
    wdm: np.ndarray,
    comp: dict[str, float],
    logrho: np.ndarray,
    logt: np.ndarray,
    apply_skye_e_offset: bool = False,
) -> dict:
    """Build (wdm, ref) arrays for each compared quantity for one mixture/set.

    Returns dict ``name -> {"w": ..., "r": ..., "kind": "rel"|"abs"}`` plus the
    auxiliary ``extra`` dict (uncorrected-comparison evidence, radiation-clean
    mask, delta_s, ion-offset dE).

    ``apply_skye_e_offset`` (default False): Skye's ion-offset bug is now fixed
    upstream (PR #1097), so Skye's mixture E is additive apart from genuine
    Coulomb mixing and the comparison uses the RAW Skye E.  Set True to restore
    the old ``dE`` reference-energy correction (kept for provenance).
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

    if table_set.startswith("skye"):
        ds = delta_s_ideal(comp)
        dE = ion_offset_dE(comp)
        extra["delta_s_ideal"] = ds
        extra["ion_offset_dE"] = dE
        S_ref = np.exp(rc("lnS"))
        S_ref_corr = S_ref - ds
        S_wdm = np.exp(wc("lnS"))
        E_wdm = np.exp(wc("lnE"))
        E_ref = np.exp(rc("lnE"))
        # (c) Skye ion-offset non-additivity.  Fixed upstream (PR #1097): the
        # mixture E is now additive, so by default compare the RAW Skye E.  The
        # dE-corrected comparison is retained as evidence (and can be made primary
        # with apply_skye_e_offset=True).
        E_ref_corr = E_ref + dE
        if apply_skye_e_offset:
            out["E"] = {"w": E_wdm, "r": E_ref_corr, "kind": "rel"}
            extra["E_raw"] = {"w": E_wdm, "r": E_ref, "kind": "rel"}
        else:
            out["E"] = {"w": E_wdm, "r": E_ref, "kind": "rel"}
            extra["E_offset_applied"] = {"w": E_wdm, "r": E_ref_corr, "kind": "rel"}
        out["S"] = {"w": S_wdm, "r": S_ref_corr, "kind": "rel"}
        # evidence: entropy comparison WITHOUT the F1 correction
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
        # (b) HELM carries NO ideal mixing entropy (mean-ion treatment), while
        # WDM_AV adds it once.  Put the HELM reference S on the WDM_AV convention
        # by ADDING Delta s_ideal.  Report both raw and corrected.
        ds = delta_s_ideal(comp)
        extra["delta_s_ideal"] = ds
        S_ref_helm = np.exp(rc("lnS"))
        out["S"] = {"w": S_wdm_gas, "r": S_ref_helm + ds, "kind": "rel"}
        extra["S_uncorrected"] = {"w": S_wdm_gas, "r": S_ref_helm, "kind": "rel"}
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
def run_one(
    table_set: str,
    mixture: str,
    *,
    figdir: Path,
    floor: dict,
    apply_skye_e_offset: bool = False,
) -> dict:
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

    cq = comparison_quantities(
        table_set, ref, wdm, comp, logrho, logt, apply_skye_e_offset
    )
    reg = regions(logrho, logt)

    summary = {
        "table_set": table_set,
        "mixture": mixture,
        "composition": comp,
        "delta_s_ideal": cq["extra"].get("delta_s_ideal"),
        "ion_offset_dE": cq["extra"].get("ion_offset_dE"),
        "n_nodes": int(logrho.size),
        "n_wdm_valid": int(wdm_valid.sum()),
        "n_ref_valid": int(ref_valid.sum()),
        "n_both_valid": int(both.sum()),
        "n_ref_valid_wdm_declined": int((ref_valid & ~wdm_valid).sum()),
        "coverage": {},
        "quantities": {},
        "uncorrected": {},
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

    # reference-correction evidence: the alternative (uncorrected / offset-applied)
    # comparisons, so the raw-vs-corrected effect is on record.
    for key, spec in cq["extra"].items():
        if isinstance(spec, dict) and "kind" in spec:
            err, ok = error_array(spec["w"], spec["r"], spec["kind"], both)
            summary["uncorrected"][key] = {
                rn: stats(err, reg[rn] & ok) for rn in REGION_NAMES
            }

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
    from matplotlib.colors import TwoSlopeNorm

    # infer the grid shape (logT outer, logRho inner)
    nR = int(np.sum(logt == logt[0]))
    nT = logt.size // nR
    LR = logrho.reshape(nT, nR)
    LT = logt.reshape(nT, nR)

    # Uniform error-map convention (all WDM_AV error maps): log10 relative error on
    # a diverging RdBu_r scale centred at -3, blue (better) to -7, red (worse) to
    # +1, saturating beyond both ends; no-data cells (declined / NaN) masked and
    # drawn light grey so they are never confused with the scale's midpoint.
    norm = TwoSlopeNorm(vmin=-7.0, vcenter=-3.0, vmax=1.0)
    cmap = plt.get_cmap("coolwarm").copy()
    cmap.set_under(cmap(0.0))
    cmap.set_over(cmap(1.0))
    cmap.set_bad("white")

    figdir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 4, figsize=(17, 8))
    axes = axes.ravel()
    im = None
    for ax, q in zip(axes, MAP_QUANTITIES):
        with np.errstate(divide="ignore"):
            lz = np.log10(maps[q].reshape(nT, nR))
        lz = np.ma.masked_invalid(lz)
        ax.set_facecolor("white")
        im = ax.pcolormesh(LR, LT, lz, shading="auto", cmap=cmap, norm=norm)
        # outline the no-data regions so white cannot be read as the gray midpoint
        _m = np.ma.getmaskarray(lz).astype(float)
        if _m.any() and not _m.all():
            ax.contour(LR, LT, _m, levels=[0.5], colors="0.3", linewidths=0.5)
        ax.axhline(6.5, color="k", lw=0.6, ls="--", alpha=0.5)
        ax.axvline(3.7, color="k", lw=0.6, ls="--", alpha=0.5)
        # target box
        ax.plot([-2, 4, 4, -2, -2], [3, 3, 7, 7, 3], color="k", lw=0.8, alpha=0.7)
        ax.set_title(f"{q}  log10|rel err|", fontsize=9)
        ax.set_xlabel("logRho")
        ax.set_ylabel("logT")
    fig.suptitle(f"WDM_AV({table_set}) vs {table_set}-on-mixture: {mixture}")
    fig.tight_layout(rect=(0, 0, 0.93, 1))
    cbar = fig.colorbar(im, ax=axes.tolist(), extend="both", fraction=0.046, pad=0.02)
    cbar.set_label("log10 |relative error|")
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
            for key, regs in s.get("uncorrected", {}).items():
                for region in REGION_NAMES:
                    st = regs.get(region)
                    if st:
                        w.writerow(
                            {
                                "table_set": s["table_set"],
                                "mixture": s["mixture"],
                                "quantity": key,
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


def _mix_median(summaries, table_set, mixture, q, region):
    by = {(s["table_set"], s["mixture"]): s for s in summaries}
    s = by.get((table_set, mixture))
    if s is None:
        return None
    return s["quantities"][q][region]["median"]


def two_rule_section(summaries: list[dict]) -> list[str]:
    """Additive volume (pressure rule) vs common-electron-density rule (D27).

    Emits, for each (base, ne) set pair present, the aggregate median side by
    side and the six mixtures where additive volume is worst (ranked by the
    pressure-rule target-gap P median), both rules shown.
    """
    qs = ("P", "E", "S", "chiT", "Cv")
    pairs = [("skye", "skye_ne", "Skye"), ("helm", "helm_ne", "HELM")]
    sets_present = {s["table_set"] for s in summaries}
    lines = ["## Two mixing rules: additive volume vs common electron density (D27)", ""]
    lines.append(
        "Median relative error against the same reference, additive volume "
        "(`pressure`) vs the common-electron-density rule (`electron_density`), "
        "for the two regions that matter. Lower is better; the common-n_e rule is "
        "expected to help against Skye (linear mixing reference) and to leave the "
        "HELM mean-ion Coulomb difference essentially unchanged."
    )
    lines.append("")
    for base, ne, tag in pairs:
        if base not in sets_present or ne not in sets_present:
            continue
        for region in ("target", "skye_valid"):
            lines.append(f"### {tag}: {_REGION_DESC[region]}")
            lines.append("")
            lines.append("| rule | " + " | ".join(qs) + " |")
            lines.append("|" + "---|" * (1 + len(qs)))
            for ts, rlabel in ((base, "additive volume"), (ne, "common n_e")):
                row = [rlabel]
                for q in qs:
                    row.append(_fmt(_agg_median(summaries, ts, q, region)))
                lines.append("| " + " | ".join(row) + " |")
            lines.append("")
        # six worst mixtures (by pressure-rule target-gap P median)
        base_s = [s for s in summaries if s["table_set"] == base]
        ranked = sorted(
            base_s,
            key=lambda s: (s["quantities"]["P"]["target"]["median"] or 0.0),
            reverse=True,
        )[:6]
        lines.append(
            f"### {tag}: six worst mixtures (ranked by additive-volume target-gap "
            "P median), both rules"
        )
        lines.append("")
        lines.append(
            "| mixture | region | rule | " + " | ".join(qs) + " |"
        )
        lines.append("|" + "---|" * (3 + len(qs)))
        for s in ranked:
            mix = s["mixture"]
            for region in ("target", "skye_valid"):
                for ts, rlabel in ((base, "AV"), (ne, "n_e")):
                    row = [mix, region, rlabel]
                    for q in qs:
                        row.append(_fmt(_mix_median(summaries, ts, mix, q, region)))
                    lines.append("| " + " | ".join(row) + " |")
        lines.append("")

    # composition series for both rules (C/O and He/C), target-gap P and E
    lines.append("### Composition series, both rules (target-gap median P and E)")
    lines.append("")
    series = {
        "C/O (minority O)": ["mix_c12_o16_99_01", "mix_c12_o16_90_10", "mix_c12_o16"],
        "He/C": ["mix_he4_c12_99_01", "mix_he4_c12_10_90", "mix_he4_c12"],
    }
    by_key = {(s["table_set"], s["mixture"]): s for s in summaries}
    for sname, stems in series.items():
        lines.append(f"#### {sname}")
        lines.append("")
        lines.append(
            "| mixture | X_min | skye P (AV) | skye P (n_e) | skye E (AV) | "
            "skye E (n_e) | helm P (AV) | helm P (n_e) |"
        )
        lines.append("|" + "---|" * 8)
        for stem in stems:
            ref = by_key.get(("skye", stem)) or by_key.get(("helm", stem))
            if ref is None:
                continue
            xmin = min(ref["composition"].values())
            row = [stem, f"{xmin:.2f}"]
            row.append(_fmt(_mix_median(summaries, "skye", stem, "P", "target")))
            row.append(_fmt(_mix_median(summaries, "skye_ne", stem, "P", "target")))
            row.append(_fmt(_mix_median(summaries, "skye", stem, "E", "target")))
            row.append(_fmt(_mix_median(summaries, "skye_ne", stem, "E", "target")))
            row.append(_fmt(_mix_median(summaries, "helm", stem, "P", "target")))
            row.append(_fmt(_mix_median(summaries, "helm_ne", stem, "P", "target")))
            lines.append("| " + " | ".join(row) + " |")
        lines.append("")
    return lines


def markdown_report(summaries: list[dict]) -> str:
    qs = ("P", "E", "S", "chiT", "Cv", "grad_ad")
    lines = [
        "# Mixture-level WDM_AV validation vs Skye / HELM",
        "",
        (
            "Relative error of the additive-volume WDM_AV mixture (integrated eosDT "
            "path, `eosDT_get_component(i_eos_WDM_AV)`) against Skye / HELM run "
            "directly on the same eleven reference mixtures, at the reference grid "
            "nodes where both sides are finite. `S` is the F1-corrected entropy for "
            "Skye (reference minus `Delta s_ideal`); for HELM, radiation is stripped "
            "from WDM_AV's e/s/Cv and gas-only chiRho/chiT come from WDM_AV's "
            "dlnPgas columns. HELM grad_ad/gamma1 mix a radiation-definition "
            "difference and are taken only where Prad/Pgas < 1e-2. Full per-mixture, "
            "per-region data (p90, p99, max, all quantities) is in the CSV/JSON."
        ),
        "",
    ]

    # (A) aggregate median table, per set and region
    lines.append(
        "## Aggregate median relative error (median over all mixtures of the set)"
    )
    lines.append("")
    for table_set, tag in (
        ("skye", "(a) Skye"),
        ("helm", "(b) HELM"),
        (
            "skye_fd",
            "(c) Skye reference, tables built from P, E, F only (FD second derivatives)",
        ),
        ("skye_ne", "(d) Skye, common-electron-density rule (D27)"),
        ("helm_ne", "(e) HELM, common-electron-density rule (D27)"),
    ):
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

    # (A2) two-rule comparison: additive volume (pressure) vs common-n_e (D27)
    lines += two_rule_section(summaries)

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

    # (B2) composition trend: target-gap median P and E vs minority mass fraction
    lines.append("## Composition trend (target gap): error vs minority mass fraction")
    lines.append("")
    lines.append(
        "Target-gap median relative error in P and E for the C/O and He/C series, "
        "ordered by the minority (smaller) mass fraction X_min. Additive volume is "
        "exact for a pure species, and the electron-fragmentation error is quadratic "
        "in the electron-density split, so for small X_min the error is expected to "
        "fall roughly in proportion to X_min. Compare P/E against "
        "`validation/fragmentation_estimate.py`."
    )
    lines.append("")
    series = {
        "C/O (minority O)": [
            "mix_c12_o16_99_01",
            "mix_c12_o16_90_10",
            "mix_c12_o16",
        ],
        "He/C (minority = smaller of He,C)": [
            "mix_he4_c12_99_01",
            "mix_he4_c12_10_90",
            "mix_he4_c12",
        ],
    }
    by_key = {(s["table_set"], s["mixture"]): s for s in summaries}
    for sname, stems in series.items():
        lines.append(f"### {sname}")
        lines.append("")
        lines.append("| mixture | X_min | skye P | skye E | helm P | helm E |")
        lines.append("|" + "---|" * 6)
        for stem in stems:
            ssk = by_key.get(("skye", stem))
            she = by_key.get(("helm", stem))
            ref = ssk or she
            if ref is None:
                continue
            xmin = min(ref["composition"].values())
            row = [stem, f"{xmin:.2f}"]
            row.append(
                _fmt(ssk["quantities"]["P"]["target"]["median"]) if ssk else "  -  "
            )
            row.append(
                _fmt(ssk["quantities"]["E"]["target"]["median"]) if ssk else "  -  "
            )
            row.append(
                _fmt(she["quantities"]["P"]["target"]["median"]) if she else "  -  "
            )
            row.append(
                _fmt(she["quantities"]["E"]["target"]["median"]) if she else "  -  "
            )
            lines.append("| " + " | ".join(row) + " |")
        lines.append("")

    # (C) reference corrections: primary vs alternative evidence
    lines.append("## Reference corrections: primary vs alternative")
    lines.append("")
    lines.append(
        "Corrections before comparison. (1) Skye S: subtract `Delta s_ideal` "
        "(Finding F1, Skye double-counts the ideal mixing entropy). (2) HELM S: "
        "ADD `Delta s_ideal` (HELM's mean-ion treatment carries no ideal mixing "
        "entropy, while WDM_AV adds it once). (3) Skye E: the ion-offset bug is "
        "fixed upstream (PR #1097), so the Skye mixture E is now additive and the "
        "comparison uses the RAW Skye E by default; the old `dE` offset correction "
        "is retained but OFF. The `primary` column is the reported comparison; "
        "`alt` is the other convention (S: uncorrected; Skye E: dE-offset applied; "
        "HELM S: raw). Median relative error."
    )
    lines.append("")
    lines.append(
        "| set | mixture | quantity | Delta_s / dE | primary (skye_valid) | "
        "alt (skye_valid) | primary (target) | alt (target) |"
    )
    lines.append("|" + "---|" * 8)
    for s in summaries:
        unc = s.get("uncorrected", {})
        if s["table_set"].startswith("skye"):
            su = unc.get("S_uncorrected", {})
            lines.append(
                f"| skye | {s['mixture']} | S | {s['delta_s_ideal']:.3e} | "
                f"{_fmt(s['quantities']['S']['skye_valid']['median'])} | "
                f"{_fmt(su.get('skye_valid', {}).get('median'))} | "
                f"{_fmt(s['quantities']['S']['target']['median'])} | "
                f"{_fmt(su.get('target', {}).get('median'))} |"
            )
            # alt E = whichever of E_raw / E_offset_applied is present
            ea = unc.get("E_offset_applied") or unc.get("E_raw") or {}
            lines.append(
                f"| skye | {s['mixture']} | E | {s['ion_offset_dE']:.3e} | "
                f"{_fmt(s['quantities']['E']['skye_valid']['median'])} | "
                f"{_fmt(ea.get('skye_valid', {}).get('median'))} | "
                f"{_fmt(s['quantities']['E']['target']['median'])} | "
                f"{_fmt(ea.get('target', {}).get('median'))} |"
            )
        else:
            su = unc.get("S_uncorrected", {})
            lines.append(
                f"| helm | {s['mixture']} | S | {s['delta_s_ideal']:.3e} | "
                f"{_fmt(su.get('skye_valid', {}).get('median'))} | "
                f"{_fmt(s['quantities']['S']['skye_valid']['median'])} | "
                f"{_fmt(su.get('target', {}).get('median'))} | "
                f"{_fmt(s['quantities']['S']['target']['median'])} |"
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
        "--sets",
        nargs="+",
        default=["skye", "helm", "skye_fd"],
        choices=["skye", "helm", "skye_fd", "skye_ne", "helm_ne"],
    )
    p.add_argument("--mixtures", nargs="+", default=list(MIXTURES))
    p.add_argument("--figdir", default=str(DEFAULT_FIGDIR))
    p.add_argument("--resultsdir", default=str(DEFAULT_RESULTSDIR))
    p.add_argument("--heldout", default=str(DEFAULT_HELDOUT))
    p.add_argument(
        "--apply-skye-e-offset",
        action="store_true",
        help="apply the old Skye ion-offset dE energy correction (default off: "
        "the ion-offset bug is fixed upstream, so the raw Skye E is additive)",
    )
    return p


def plot_rule_diff(base: str, ne: str, mixture: str, figdir: Path):
    """Difference-of-rules map: log10|err(n_e)| - log10|err(AV)| per quantity.

    Blue (negative) => the common-electron-density rule is more accurate than
    additive volume at that node; red => worse. Same reference for both rules.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    ref = read_table(REF_DIRS[base] / f"{mixture}.txt", len(_REF_NAMES))
    wa = read_table(OUT_DIRS[base] / f"{mixture}.txt", len(_WDM_NAMES))
    wn = read_table(OUT_DIRS[ne] / f"{mixture}.txt", len(_WDM_NAMES))
    hdr = read_header(OUT_DIRS[ne] / f"{mixture}.txt")
    comp = parse_composition(hdr["composition"])
    logrho = _col(wn, WDM_COLS, "logRho")
    logt = _col(wn, WDM_COLS, "logT")

    both_a = (_col(wa, WDM_COLS, "flag") == 0.0) & np.isfinite(
        _col(ref, REF_COLS, "lnPgas")
    )
    both_n = (_col(wn, WDM_COLS, "flag") == 0.0) & np.isfinite(
        _col(ref, REF_COLS, "lnPgas")
    )
    both = both_a & both_n

    cqa = comparison_quantities(base, ref, wa, comp, logrho, logt)["quantities"]
    cqn = comparison_quantities(ne, ref, wn, comp, logrho, logt)["quantities"]

    nR = int(np.sum(logt == logt[0]))
    nT = logt.size // nR
    LR = logrho.reshape(nT, nR)
    LT = logt.reshape(nT, nR)
    norm = TwoSlopeNorm(vmin=-3.0, vcenter=0.0, vmax=3.0)
    cmap = plt.get_cmap("coolwarm_r").copy()  # blue = ne better (negative)
    cmap.set_bad("white")

    figdir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    axes = axes.ravel()
    im = None
    for ax, q in zip(axes, ("P", "E", "S", "chiT", "Cv", "grad_ad")):
        ea, oka = error_array(cqa[q]["w"], cqa[q]["r"], cqa[q]["kind"], both)
        en, okn = error_array(cqn[q]["w"], cqn[q]["r"], cqn[q]["kind"], both)
        with np.errstate(divide="ignore", invalid="ignore"):
            diff = np.log10(en) - np.log10(ea)
        diff = np.ma.masked_invalid(diff.reshape(nT, nR))
        ax.set_facecolor("white")
        im = ax.pcolormesh(LR, LT, diff, shading="auto", cmap=cmap, norm=norm)
        ax.axhline(6.5, color="k", lw=0.6, ls="--", alpha=0.5)
        ax.axvline(3.7, color="k", lw=0.6, ls="--", alpha=0.5)
        ax.plot([-2, 4, 4, -2, -2], [3, 3, 7, 7, 3], color="k", lw=0.8, alpha=0.7)
        ax.set_title(f"{q}: log10|err(n_e)| - log10|err(AV)|", fontsize=9)
        ax.set_xlabel("logRho")
        ax.set_ylabel("logT")
    fig.suptitle(f"{base}: rule difference (blue = common-n_e better): {mixture}")
    fig.tight_layout(rect=(0, 0, 0.92, 1))
    cbar = fig.colorbar(im, ax=axes.tolist(), extend="both", fraction=0.046, pad=0.02)
    cbar.set_label("log10 error ratio (n_e / AV); negative = n_e better")
    path = figdir / f"rulediff_{base}_{mixture}.png"
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    floor = load_heldout_floor(Path(args.heldout))
    figdir = Path(args.figdir)
    summaries = []
    for table_set in args.sets:
        for mixture in args.mixtures:
            print(f"mixtures: {table_set}  {mixture}")
            summary = run_one(
                table_set,
                mixture,
                figdir=figdir,
                floor=floor,
                apply_skye_e_offset=args.apply_skye_e_offset,
            )
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

    # difference-of-rules maps (skye set) for the key mixtures, when both the
    # additive-volume and common-n_e Skye outputs are present
    if "skye" in args.sets and "skye_ne" in args.sets:
        for mixture in ("mix_c12_o16", "mix_he4_c12", "mix_he4_c12_o16"):
            if mixture in args.mixtures:
                try:
                    p = plot_rule_diff("skye", "skye_ne", mixture, figdir)
                    print("wrote rule-diff map", p)
                except (FileNotFoundError, ValueError) as e:
                    print(f"  (skip rule-diff {mixture}: {e})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
