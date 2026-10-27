"""Nodal-data builder: turn a :class:`~wdm_av_builder.io.SourceTable` into the
nine Hermite nodal fields the processed format needs, in log coordinates
(xi = ln rho, tau = ln T).

Exact relations (FORMAT.md Section 2; eos_mixing_notes.tex Sec. 2, 8)::

    f      = F                      (the tabulated Helmholtz free energy)
    f_xi   = P / rho               (exact from the source P)
    f_tau  = -T s = f - e          (prefer S; else F - E)

Second derivatives, when the source supplies dPdRho_T, dPdT_Rho, Cv::

    f_xixi   = (dP/drho)_T - f_xi          [see note below]
    f_xitau  = T (dP/dT)_rho / rho
    f_tautau = f_tau - T c_V

.. note::

   FORMAT.md line 105 and the WP1 brief both write ``f_xixi = rho (dP/drho)_T
   - f_xi``.  That ``rho`` factor is dimensionally inconsistent (``rho dP/drho``
   is erg/cm^3, not erg/g) and gives the wrong value: for an ideal gas
   ``(dP/drho)_T = f_xi`` so ``f_xixi`` must be 0, which only the factor-free
   form delivers.  Derivation: ``f_xi = P/rho`` and
   ``f_xixi = d/dxi(P/rho) = rho d/drho(P/rho) = (dP/drho)_T - P/rho``.  This
   builder uses the correct relation ``f_xixi = (dP/drho)_T - f_xi``.

Otherwise the second and the three higher cross derivatives come from
fourth-order (or better) finite differences on the uniform grid, using
Fornberg weights (centered in the interior, one-sided at edges).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import constants as C
from .io import SourceTable
from .synthetic import ELEMENTS

# ---------------------------------------------------------------------------
# generic finite-difference machinery (Fornberg 1988)
# ---------------------------------------------------------------------------


def fornberg_weights(z: float, x: np.ndarray, m: int) -> np.ndarray:
    """Fornberg (1988) finite-difference weights.

    Returns an array ``c`` of shape ``(len(x), m + 1)`` such that the k-th
    derivative at ``z`` is ``sum_i c[i, k] * f(x[i])`` for k = 0..m, exact for
    polynomials up to degree ``len(x) - 1``.

    Parameters
    ----------
    z : point at which the derivatives are approximated.
    x : 1-D array of (distinct) node locations.
    m : highest derivative order wanted.
    """
    x = np.asarray(x, dtype=float)
    n = x.size - 1
    c = np.zeros((n + 1, m + 1))
    c1 = 1.0
    c4 = x[0] - z
    c[0, 0] = 1.0
    for i in range(1, n + 1):
        mn = min(i, m)
        c2 = 1.0
        c5 = c4
        c4 = x[i] - z
        for j in range(i):
            c3 = x[i] - x[j]
            c2 *= c3
            if j == i - 1:
                for k in range(mn, 0, -1):
                    c[i, k] = c1 * (k * c[i - 1, k - 1] - c5 * c[i - 1, k]) / c2
                c[i, 0] = -c1 * c5 * c[i - 1, 0] / c2
            for k in range(mn, 0, -1):
                c[j, k] = (c4 * c[j, k] - k * c[j, k - 1]) / c3
            c[j, 0] = c4 * c[j, 0] / c3
        c1 = c2
    return c


def fd_weights(offsets: np.ndarray, deriv: int) -> np.ndarray:
    """Weights for ``deriv``-th derivative at 0 on integer ``offsets`` (unit
    spacing).  Divide the result by ``h**deriv`` for spacing ``h``."""
    offsets = np.asarray(offsets, dtype=float)
    return fornberg_weights(0.0, offsets, deriv)[:, deriv]


def _stencil_npts(deriv: int, acc: int) -> int:
    """Number of points for accuracy ``acc`` of derivative order ``deriv``,
    rounded up to odd so interior stencils can be centered."""
    npts = deriv + acc
    if npts % 2 == 0:
        npts += 1
    return npts


def apply_fd(
    arr: np.ndarray, axis: int, deriv: int, h: float, acc: int = 4
) -> np.ndarray:
    """Finite-difference derivative of ``arr`` along ``axis`` on a uniform grid.

    Centered Fornberg stencils in the interior, sliding to one-sided at the
    edges.  First derivatives use a 5-point stencil (4th order everywhere,
    matching eos_mixing_notes.tex Sec. 9.4); second derivatives use a 7-point
    stencil (>= 4th order everywhere).  Vectorised: one matmul per stencil
    group, no Python loop over the off-axis dimension.
    """
    arr = np.asarray(arr, dtype=float)
    n = arr.shape[axis]
    npts = _stencil_npts(deriv, acc)
    if n < npts:
        raise ValueError(f"axis length {n} < stencil size {npts}")
    work = np.moveaxis(arr, axis, 0)  # (n, ...)
    out = np.empty_like(work)
    # stencil start index for each point, clamped into range
    starts = np.clip(np.arange(n) - npts // 2, 0, n - npts)
    rel = starts - np.arange(n)  # stencil offset relative to the point
    # all interior points share one offset group; only edges differ, so this is
    # O(npts) vectorised tensordots, with no per-node Python loop.
    for r in np.unique(rel):
        rows = np.nonzero(rel == r)[0]
        offs = np.arange(npts) + r
        w = fd_weights(offs, deriv) / h**deriv
        idx = rows[:, None] + offs[None, :]  # (nrows, npts) in-range by construction
        gathered = work[idx]  # (nrows, npts, ...)
        out[rows] = np.tensordot(gathered, w, axes=([1], [0]))
    return np.moveaxis(out, 0, axis)


def _d_dxi(a, hxi):
    return apply_fd(a, axis=1, deriv=1, h=hxi)


def _d_dtau(a, htau):
    return apply_fd(a, axis=0, deriv=1, h=htau)


def _d2_dxi2(a, hxi):
    return apply_fd(a, axis=1, deriv=2, h=hxi)


def _d2_dtau2(a, htau):
    return apply_fd(a, axis=0, deriv=2, h=htau)


# ---------------------------------------------------------------------------
# nodal-data builder
# ---------------------------------------------------------------------------
@dataclass
class NodalData:
    """Result of :func:`build_nodal_data`."""

    attrs: dict = field(default_factory=dict)
    datasets: dict[str, np.ndarray] = field(default_factory=dict)
    diagnostics: dict[str, np.ndarray] = field(default_factory=dict)


def _check_uniform(logv: np.ndarray, name: str) -> float:
    d = np.diff(logv)
    h = float(d.mean())
    if not np.allclose(d, h, rtol=1e-9, atol=0):
        raise ValueError(f"{name} grid not uniform (spread {d.ptp():.3e})")
    return h


def build_nodal_data(
    source: SourceTable,
    *,
    order: int = 5,
    fxt_route: int = 1,
    zero_twist: bool = False,
    e_offset: float = 0.0,
    f_offset: float = 0.0,
    ignore_source_derivs: bool = False,
) -> NodalData:
    """Build the 9 nodal fields (and optional Zbar/eta) from a source table.

    Parameters
    ----------
    order : 3 (bicubic) or 5 (biquintic).  For order 3 the second and higher
        cross derivatives are set to zero (only f, f_xi, f_tau, f_xitau matter).
    fxt_route : which estimate of f_xitau to store when it is not supplied by
        the source: 1 = d/dtau(P/rho); 2 = d/dxi(f - e); 3 = cross FD of f.
        Route 2' (P/rho - d/dxi e) is always computed for the diagnostics.
    zero_twist : set the three higher cross terms
        (f_xixitau, f_xitautau, f_xixitautau) to zero (experiment), recorded in
        ``fd_scheme``.
    e_offset, f_offset : constants added to the (post-radiation-strip) E and F.
    ignore_source_derivs : when True, never use the supplied ``dPdRho_T``,
        ``dPdT_Rho``, ``Cv`` columns for the second derivatives; take the
        finite-difference path instead (so ``second_derivs_from_source = 0``).
        Used to compare the source-derivative and pure-FD builds.
    """
    if order not in (3, 5):
        raise ValueError("order must be 3 or 5")
    if fxt_route not in (1, 2, 3):
        raise ValueError("fxt_route must be 1, 2 or 3")

    logrho, logt = source.logrho, source.logt
    xi = C.ln10 * logrho
    tau = C.ln10 * logt
    hxi = _check_uniform(xi, "xi")
    htau = _check_uniform(tau, "tau")
    nxi, ntau = logrho.size, logt.size

    rho = np.exp(xi)[None, :]
    T = np.exp(tau)[:, None]

    P = np.array(source.columns["P"], dtype=float)
    E = np.array(source.columns["E"], dtype=float) if "E" in source.columns else None
    F = np.array(source.columns["F"], dtype=float) if "F" in source.columns else None
    S = np.array(source.columns["S"], dtype=float) if "S" in source.columns else None
    if F is None and (E is None or S is None):
        raise ValueError("source must contain F, or both E and S")

    # --- radiation stripping (analytic) ---
    if source.includes_radiation:
        a = C.crad
        P = P - a * T**4 / 3.0
        if E is not None:
            E = E - a * T**4 / rho
        if F is not None:
            F = F + a * T**4 / (3.0 * rho)
        if S is not None:
            S = S - (4.0 / 3.0) * a * T**3 / rho

    # --- energy zero-point offsets ---
    if E is not None:
        E = E + e_offset
    if F is not None:
        F = F + f_offset

    # --- exact first derivatives ---
    f_xi = P / rho
    if S is not None:
        f_tau = -T * S
    else:
        f_tau = F - E
    if F is not None:
        f = F
    else:  # reconstruct from e and f_tau = f - e
        f = E + f_tau
    e = f - f_tau if E is None else E

    # --- second derivatives ---
    have_src2 = (not ignore_source_derivs) and all(
        c in source.columns for c in ("dPdRho_T", "dPdT_Rho", "Cv")
    )
    if have_src2:
        dPdRho_T = np.array(source.columns["dPdRho_T"], dtype=float)
        dPdT_Rho = np.array(source.columns["dPdT_Rho"], dtype=float)
        Cv = np.array(source.columns["Cv"], dtype=float)
        f_xixi = dPdRho_T - f_xi  # corrected relation (see module docstring)
        f_xitau = T * dPdT_Rho / rho
        f_tautau = f_tau - T * Cv
        second_from_source = 1
    else:
        f_xixi = _d_dxi(f_xi, hxi)
        f_tautau = _d_dtau(f_tau, htau)
        f_xitau = None  # chosen below
        second_from_source = 0

    # --- f_xitau routes (always computed for diagnostics) ---
    r1 = _d_dtau(f_xi, htau)  # d/dtau (P/rho)
    r2 = _d_dxi(f_tau, hxi)  # d/dxi (f - e)
    r2p = f_xi - (_d_dxi(e, hxi) if e is not None else 0.0)  # P/rho - d_xi e
    r3 = _d_dxi(_d_dtau(f, htau), hxi)  # cross FD of f
    diagnostics = {
        "fxitau_route1": r1,
        "fxitau_route2": r2,
        "fxitau_route2p": r2p,
        "fxitau_route3": r3,
        "maxwell_residual": r1 - r2,
    }
    if not have_src2:
        f_xitau = {1: r1, 2: r2, 3: r3}[fxt_route]

    # --- higher cross terms (biquintic) from P/rho = f_xi ---
    f_xitautau = _d2_dtau2(f_xi, htau)  # d^2/dtau^2 (P/rho)
    d_tau_fxi = r1
    f_xixitau = _d_dxi(d_tau_fxi, hxi)  # d/dxi d/dtau (P/rho)
    f_xixitautau = _d_dxi(f_xitautau, hxi)  # d/dxi d^2/dtau^2 (P/rho)
    # cross-check routes via f_tau = f - e
    diagnostics["fxixitau_alt"] = _d2_dxi2(f_tau, hxi)  # d^2/dxi^2 (f - e)
    diagnostics["fxitautau_alt"] = _d_dxi(_d_dtau(f_tau, htau), hxi)
    diagnostics["fxixitautau_alt"] = _d_dxi(_d2_dxi2(f_tau, hxi), htau)
    diagnostics["fxixitau_disagree"] = f_xixitau - diagnostics["fxixitau_alt"]
    diagnostics["fxitautau_disagree"] = f_xitautau - diagnostics["fxitautau_alt"]
    diagnostics["fxixitautau_disagree"] = f_xixitautau - diagnostics["fxixitautau_alt"]

    # stability is judged on the *real* second derivatives (source or FD),
    # before any order-3 / zero-twist zeroing of the stored datasets.
    P_rho_stab = f_xi + f_xixi
    cv_stab = (f_tau - f_tautau) / T

    fd_notes = []
    if ignore_source_derivs and all(
        c in source.columns for c in ("dPdRho_T", "dPdT_Rho", "Cv")
    ):
        fd_notes.append("forced_fd_ignore_source_derivs")
    if zero_twist:
        f_xixitau = np.zeros_like(f_xi)
        f_xitautau = np.zeros_like(f_xi)
        f_xixitautau = np.zeros_like(f_xi)
        fd_notes.append("twist_terms_zeroed")

    if order == 3:
        zero = np.zeros_like(f_xi)
        f_xixi = zero.copy()
        f_tautau = zero.copy()
        f_xixitau = zero.copy()
        f_xitautau = zero.copy()
        f_xixitautau = zero.copy()
        fd_notes.append("bicubic_order3")

    # --- valid mask: P > 0, all finite, and mechanically/thermally stable
    #     (P_rho = f_xi + f_xixi > 0 and c_V = (f_tau - f_tautau)/T > 0) ---
    stack = [
        f,
        f_xi,
        f_tau,
        f_xixi,
        f_tautau,
        f_xitau,
        f_xixitau,
        f_xitautau,
        f_xixitautau,
    ]
    finite = np.ones_like(f, dtype=bool)
    for arr in stack:
        finite &= np.isfinite(arr)
    finite &= np.isfinite(P_rho_stab) & np.isfinite(cv_stab)
    valid = (P > 0.0) & finite & (P_rho_stab > 0.0) & (cv_stab > 0.0)

    datasets = {
        "f": f,
        "f_xi": f_xi,
        "f_tau": f_tau,
        "f_xixi": f_xixi,
        "f_tautau": f_tautau,
        "f_xitau": f_xitau,
        "f_xixitau": f_xixitau,
        "f_xitautau": f_xitautau,
        "f_xixitautau": f_xixitautau,
        "valid": valid.astype(np.int8),
    }

    has_Zbar = 0
    if "Zbar" in source.columns:
        zb = np.array(source.columns["Zbar"], dtype=float)
        datasets["Zbar"] = zb
        datasets["Zbar_xi"] = _d_dxi(zb, hxi)
        datasets["Zbar_tau"] = _d_dtau(zb, htau)
        datasets["Zbar_xitau"] = _d_dxi(_d_dtau(zb, htau), hxi)
        has_Zbar = 1
    has_eta = 0
    if "eta" in source.columns:
        et = np.array(source.columns["eta"], dtype=float)
        datasets["eta"] = et
        datasets["eta_xi"] = _d_dxi(et, hxi)
        datasets["eta_tau"] = _d_dtau(et, htau)
        datasets["eta_xitau"] = _d_dxi(_d_dtau(et, htau), hxi)
        has_eta = 1

    npts1 = _stencil_npts(1, 4)
    npts2 = _stencil_npts(2, 4)
    fd_scheme = (
        f"fornberg uniform; d1={npts1}pt(4th), d2={npts2}pt(>=4th); "
        f"f_xitau=route{fxt_route}; higher cross from d(P/rho)"
    )
    if fd_notes:
        fd_scheme += "; " + ",".join(fd_notes)

    element = source.header.get("element", ELEMENTS.get(source.Z, f"Z{source.Z}"))
    attrs = {
        "format_version": 1,
        "isotope": source.name,
        "element": element,
        "Z": source.Z,
        "A": source.A,
        "interp_order": order,
        "nxi": nxi,
        "ntau": ntau,
        "xi_min": float(xi[0]),
        "dxi": hxi,
        "tau_min": float(tau[0]),
        "dtau": htau,
        "includes_radiation": 0,
        "e_offset": float(e_offset),
        "f_offset": float(f_offset),
        "has_Zbar": has_Zbar,
        "has_eta": has_eta,
        "fd_scheme": fd_scheme,
        "second_derivs_from_source": second_from_source,
    }
    return NodalData(attrs=attrs, datasets=datasets, diagnostics=diagnostics)
