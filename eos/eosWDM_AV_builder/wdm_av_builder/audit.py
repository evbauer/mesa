"""Table-consistency audits (eos_mixing_notes.tex Sec. 9.3 and 10).

The checks operate on a :class:`~wdm_av_builder.io.SourceTable` together with
the nodal data built from it:

* **precision** -- ``F - E + T S`` (which should vanish if the table is
  thermodynamically consistent); reports the printed precision of the table.
  ``n/a`` when no S column is present.
* **Euler-Maclaurin** residuals of ``f`` against ``P/rho`` along xi and against
  ``f_tau`` along tau (trapezoid-rule residual, relative to the trapezoid term).
* **Maxwell** residual between route 1 (``d/dtau(P/rho)``) and route 2
  (``d/dxi(f - e)``) for ``f_xitau``.
* **FD stability** maps ``P_rho = f_xi + f_xixi`` and
  ``c_V = (f_tau - f_tautau)/T``, which must be positive.
* **Bicubic node-jump** maps for chi_rho (along xi) and c_V (along tau), from
  the cubic-Hermite second-derivative endpoint formulas of Sec. 8.2.  These are
  fixed by the tabulated ``f`` and first derivatives, independent of how
  ``f_xitau`` is estimated.
"""

from __future__ import annotations

import numpy as np

from . import constants as C
from .derivs import NodalData, build_nodal_data
from .io import SourceTable

_TINY = 1e-300

# audit maps written into the processed file's ``audit/`` group (FORMAT.md 2).
PROCESSED_AUDIT_KEYS: tuple[str, ...] = (
    "maxwell_residual",
    "euler_maclaurin_xi",
    "euler_maclaurin_tau",
    "chi_rho_jump",
    "cv_jump",
    "P_rho_fd",
    "cv_fd",
    "precision",
)


def _euler_maclaurin(f: np.ndarray, fp: np.ndarray, h: float, axis: int):
    """Trapezoid-rule residual of f against its derivative fp along ``axis``.

    Returns a map the same shape as f (cell residual stored at the left node of
    each cell; last index along the axis left zero), relative to the trapezoid
    term.
    """
    f = np.moveaxis(f, axis, -1)
    fp = np.moveaxis(fp, axis, -1)
    df = f[..., 1:] - f[..., :-1]
    trap = 0.5 * h * (fp[..., :-1] + fp[..., 1:])
    res = df - trap
    # robust relative measure: small (~O(h)) on a smooth consistent table,
    # large where the table is inconsistent/under-resolved; stays finite at
    # zero crossings of the derivative (unlike dividing by |trap| alone).
    rel = res / (np.abs(df) + np.abs(trap) + _TINY)
    out = np.zeros_like(f)
    out[..., :-1] = rel
    return np.moveaxis(out, -1, axis)


def _node_jump(f: np.ndarray, fp: np.ndarray, h: float, axis: int):
    """Bicubic second-derivative node-jump map along ``axis``.

    ``fp`` is the exact first derivative on grid lines (f_xi along xi, f_tau
    along tau).  Returns the jump in the Hermite second derivative at each
    interior node (same shape as f; boundaries zero).
    """
    f = np.moveaxis(f, axis, -1)
    fp = np.moveaxis(fp, axis, -1)
    diff = f[..., 1:] - f[..., :-1]  # per cell
    # cubic-Hermite d2f/dt2 at the two cell ends (Sec. 8.2)
    right_cell = 6.0 * diff - 4.0 * h * fp[..., :-1] - 2.0 * h * fp[..., 1:]
    left_cell = -6.0 * diff + 2.0 * h * fp[..., :-1] + 4.0 * h * fp[..., 1:]
    out = np.zeros_like(f)
    # node i (interior): right cell at t=0 minus left cell at t=1, scaled h^-2
    out[..., 1:-1] = (right_cell[..., 1:] - left_cell[..., :-1]) / h**2
    return np.moveaxis(out, -1, axis)


def _stats(arr: np.ndarray, mask: np.ndarray | None, thresh: float) -> dict:
    a = np.abs(np.asarray(arr, dtype=float))
    if mask is not None:
        a = a[mask]
    a = a[np.isfinite(a)]
    if a.size == 0:
        return {
            "max": float("nan"),
            "median": float("nan"),
            "frac_flagged": float("nan"),
        }
    return {
        "max": float(a.max()),
        "median": float(np.median(a)),
        "frac_flagged": float(np.mean(a > thresh)),
    }


def run_audit(
    source: SourceTable,
    nodal: NodalData | None = None,
    *,
    flag_threshold: float = 1e-3,
) -> tuple[dict[str, np.ndarray], dict[str, dict]]:
    """Run all audits.  Returns ``(maps, summary)``.

    ``maps`` holds 2-D arrays shaped like the grid; ``summary`` holds per-map
    max / median / fraction-flagged (and positivity fractions for the stability
    maps).
    """
    if nodal is None:
        nodal = build_nodal_data(source)
    ds = nodal.datasets
    xi = C.ln10 * source.logrho
    tau = C.ln10 * source.logt
    hxi = float(xi[1] - xi[0])
    htau = float(tau[1] - tau[0])
    T = np.exp(tau)[:, None]
    valid = ds["valid"].astype(bool)

    f = ds["f"]
    f_xi = ds["f_xi"]
    f_tau = ds["f_tau"]

    maps: dict[str, np.ndarray] = {}

    # --- precision (needs independent F, E, S) ---
    cols = source.columns
    if "S" in cols and "E" in cols and "F" in cols:
        resid = cols["F"] - cols["E"] + T * cols["S"]
        scale = np.abs(cols["F"]) + np.abs(cols["E"]) + _TINY
        maps["precision"] = resid / scale
        precision_available = True
    else:
        maps["precision"] = np.full(f.shape, np.nan)
        precision_available = False

    # --- Euler-Maclaurin ---
    maps["euler_maclaurin_xi"] = _euler_maclaurin(f, f_xi, hxi, axis=1)
    maps["euler_maclaurin_tau"] = _euler_maclaurin(f, f_tau, htau, axis=0)

    # --- Maxwell residual (route 1 vs route 2), relative ---
    r1 = nodal.diagnostics["fxitau_route1"]
    r2 = nodal.diagnostics["fxitau_route2"]
    maps["maxwell_residual"] = (r1 - r2) / (np.abs(r1) + np.abs(r2) + _TINY)

    # --- FD stability ---
    P_rho_fd = f_xi + ds["f_xixi"]
    cv_fd = (f_tau - ds["f_tautau"]) / T
    maps["P_rho_fd"] = P_rho_fd
    maps["cv_fd"] = cv_fd

    # --- bicubic node jumps (fractional: jump in the quantity / the quantity) ---
    chi_rho = 1.0 + ds["f_xixi"] / (f_xi + _TINY)
    fxixi_jump = _node_jump(f, f_xi, hxi, axis=1)
    maps["chi_rho_jump"] = (fxixi_jump / (f_xi + _TINY)) / (chi_rho + _TINY)
    ftautau_jump = _node_jump(f, f_tau, htau, axis=0)
    maps["cv_jump"] = (-ftautau_jump / T) / (cv_fd + _TINY)

    # --- summary ---
    summary: dict[str, dict] = {}
    summary["precision"] = (
        _stats(maps["precision"], valid, flag_threshold)
        if precision_available
        else {"status": "n/a (no S column)"}
    )
    for key in (
        "euler_maclaurin_xi",
        "euler_maclaurin_tau",
        "maxwell_residual",
        "chi_rho_jump",
        "cv_jump",
    ):
        summary[key] = _stats(maps[key], valid, flag_threshold)
    summary["P_rho_fd"] = {
        "min": float(np.nanmin(maps["P_rho_fd"][valid])),
        "frac_nonpositive": float(np.mean(maps["P_rho_fd"][valid] <= 0.0)),
    }
    summary["cv_fd"] = {
        "min": float(np.nanmin(maps["cv_fd"][valid])),
        "frac_nonpositive": float(np.mean(maps["cv_fd"][valid] <= 0.0)),
    }
    summary["valid"] = {"frac_valid": float(valid.mean())}
    return maps, summary


def format_summary(name: str, summary: dict[str, dict]) -> str:
    """One-screen text report of an audit summary."""
    lines = [f"audit summary for {name}:"]
    for key, st in summary.items():
        parts = ", ".join(
            f"{k}={v:.3e}" if isinstance(v, float) else f"{k}={v}"
            for k, v in st.items()
        )
        lines.append(f"  {key:20s} {parts}")
    return "\n".join(lines)
