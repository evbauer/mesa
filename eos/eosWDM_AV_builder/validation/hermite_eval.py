"""Independent, vectorised numpy Hermite tensor-product evaluator for the
WDM_AV processed tables.

This module is written *from the mathematics* of the local Hermite
tensor-product interpolant (``eos_mixing_notes.tex`` Sec. 8, and the standard
quintic/cubic Hermite bases), not transcribed from the Fortran reader.  It is
used to validate that reader and the builder: evaluating a processed table at
its own nodes must reproduce the stored nodal data to roundoff, and the
analytic derivatives of the interpolant must satisfy the cell-corner identity
``integral(f_xitau) == corner differences of f``.

Coordinates are the natural logs ``xi = ln rho`` and ``tau = ln T``.  A
processed table is read with h5py per ``FORMAT.md`` Section 2: datasets are
stored ``(ntau, nxi)`` (tau is the slow axis), the grid is uniform with origin
``(xi_min, tau_min)`` and spacing ``(dxi, dtau)``.

Interpolant (one cell ``[xi_i, xi_{i+1}] x [tau_j, tau_{j+1}]``; local
coordinates ``t = (xi-xi_i)/h_xi``, ``u = (tau-tau_j)/h_tau``):

    f(t,u) = sum_{a,b in {L,R}} sum_{p,q} h_xi^p h_tau^q
             (d^{p+q} f / dxi^p dtau^q)_{a,b}  B^a_p(t) B^b_q(u)

with, for the quintic (order 5) scheme, p,q in {0,1,2} and the left/right 1-D
bases

    B^L_p(t) = phi_p(t),    B^R_p(t) = (-1)^p phi_p(1-t),

where phi_0, phi_1, phi_2 are the quintic Hermite basis polynomials
(``phi_k`` carries the k-th derivative at t=0).  The sign ``(-1)^p`` on the
right-node basis is the standard mirror of an odd-order Hermite weight under
``t -> 1-t`` (it is exactly the sign flip the Fortran ``qwts`` applies to the
``psi1`` terms).  The cubic (order 3) scheme is the same with p,q in {0,1} and
the cubic bases xpsi0, xpsi1.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

# ---------------------------------------------------------------------------
# 1-D Hermite basis polynomials and their first two derivatives (in t).
# Quintic basis: phi_k has unit k-th derivative at t=0, all other value / 1st /
# 2nd derivatives zero at both ends.  (eos_mixing_notes.tex Sec. 8.4.)
# ---------------------------------------------------------------------------


def _phi0(t):
    return 1.0 - 10.0 * t**3 + 15.0 * t**4 - 6.0 * t**5


def _phi1(t):
    return t - 6.0 * t**3 + 8.0 * t**4 - 3.0 * t**5


def _phi2(t):
    return 0.5 * t**2 - 1.5 * t**3 + 1.5 * t**4 - 0.5 * t**5


def _dphi0(t):
    return -30.0 * t**2 + 60.0 * t**3 - 30.0 * t**4


def _dphi1(t):
    return 1.0 - 18.0 * t**2 + 32.0 * t**3 - 15.0 * t**4


def _dphi2(t):
    return t - 4.5 * t**2 + 6.0 * t**3 - 2.5 * t**4


def _ddphi0(t):
    return -60.0 * t + 180.0 * t**2 - 120.0 * t**3


def _ddphi1(t):
    return -36.0 * t + 96.0 * t**2 - 60.0 * t**3


def _ddphi2(t):
    return 1.0 - 9.0 * t + 18.0 * t**2 - 10.0 * t**3


# Cubic basis (eos_mixing_notes.tex Sec. 8.3): xpsi0 carries the value, xpsi1
# the first derivative, at t=0.
def _xpsi0(t):
    return 2.0 * t**3 - 3.0 * t**2 + 1.0


def _xpsi1(t):
    return t**3 - 2.0 * t**2 + t


def _dxpsi0(t):
    return 6.0 * t**2 - 6.0 * t


def _dxpsi1(t):
    return 3.0 * t**2 - 4.0 * t + 1.0


def _ddxpsi0(t):
    return 12.0 * t - 6.0


def _ddxpsi1(t):
    return 6.0 * t - 4.0


def _quintic_weights(t: np.ndarray, h: float):
    """Return value / first- / second-derivative weight arrays, each shape
    ``(N, 6)``, for the 6 quintic nodal slots ordered

        (L value, L d1, L d2, R value, R d1, R d2),

    where the derivatives are with respect to the *physical* coordinate.  A dot
    of the value weights with the nodal-data vector gives f; using the first-
    (second-) derivative weights gives df/dx (d2f/dx2)."""
    t = np.asarray(t, dtype=float)
    mt = 1.0 - t
    h2 = h * h
    # value weights
    w0 = np.stack(
        [
            _phi0(t),
            h * _phi1(t),
            h2 * _phi2(t),
            _phi0(mt),
            -h * _phi1(mt),
            h2 * _phi2(mt),
        ],
        axis=-1,
    )
    # first derivative wrt physical coordinate (chain rule d/dx = (1/h) d/dt,
    # and d/dx of f(1-t) = -(1/h) f'(1-t))
    w1 = np.stack(
        [
            _dphi0(t) / h,
            _dphi1(t),
            h * _dphi2(t),
            -_dphi0(mt) / h,
            _dphi1(mt),
            -h * _dphi2(mt),
        ],
        axis=-1,
    )
    # second derivative wrt physical coordinate
    w2 = np.stack(
        [
            _ddphi0(t) / h2,
            _ddphi1(t) / h,
            _ddphi2(t),
            _ddphi0(mt) / h2,
            -_ddphi1(mt) / h,
            _ddphi2(mt),
        ],
        axis=-1,
    )
    return w0, w1, w2


def _cubic_weights(t: np.ndarray, h: float):
    """As :func:`_quintic_weights` but for the 4 cubic nodal slots ordered
    ``(L value, L d1, R value, R d1)``."""
    t = np.asarray(t, dtype=float)
    mt = 1.0 - t
    h2 = h * h
    w0 = np.stack([_xpsi0(t), h * _xpsi1(t), _xpsi0(mt), -h * _xpsi1(mt)], axis=-1)
    w1 = np.stack([_dxpsi0(t) / h, _dxpsi1(t), -_dxpsi0(mt) / h, _dxpsi1(mt)], axis=-1)
    w2 = np.stack(
        [_ddxpsi0(t) / h2, _ddxpsi1(t) / h, _ddxpsi0(mt) / h2, -_ddxpsi1(mt) / h],
        axis=-1,
    )
    return w0, w1, w2


# nodal field for each (xi-order p, tau-order q); the corners are gathered from
# these 2-D arrays.
_FIELD_Q = {
    (0, 0): "f",
    (1, 0): "f_xi",
    (2, 0): "f_xixi",
    (0, 1): "f_tau",
    (0, 2): "f_tautau",
    (1, 1): "f_xitau",
    (2, 1): "f_xixitau",
    (1, 2): "f_xitautau",
    (2, 2): "f_xixitautau",
}
_FIELD_C = {
    (0, 0): "f",
    (1, 0): "f_xi",
    (0, 1): "f_tau",
    (1, 1): "f_xitau",
}


@dataclass
class HermiteTable:
    """A processed WDM_AV table ready for evaluation."""

    order: int
    xi_min: float
    dxi: float
    nxi: int
    tau_min: float
    dtau: float
    ntau: int
    datasets: dict[str, np.ndarray]
    attrs: dict

    @property
    def norder(self) -> int:
        """Number of 1-D nodal slots per direction (3 cubic, 6 quintic... i.e.
        4 or 6)."""
        return 6 if self.order == 5 else 4


def table_from_arrays(attrs: dict, datasets: dict[str, np.ndarray]) -> HermiteTable:
    """Wrap in-memory nodal data (e.g. from the builder's ``NodalData``) as a
    HermiteTable, without a round trip through HDF5.  ``attrs`` and ``datasets``
    follow FORMAT.md Section 2.
    """
    return HermiteTable(
        order=int(attrs["interp_order"]),
        xi_min=float(attrs["xi_min"]),
        dxi=float(attrs["dxi"]),
        nxi=int(attrs["nxi"]),
        tau_min=float(attrs["tau_min"]),
        dtau=float(attrs["dtau"]),
        ntau=int(attrs["ntau"]),
        datasets={k: np.asarray(v, dtype=float) for k, v in datasets.items()},
        attrs=dict(attrs),
    )


def load_table(path: str | Path) -> HermiteTable:
    """Read a processed HDF5 table (FORMAT.md Sec. 2) into a HermiteTable.

    Only the fields the evaluator needs are loaded; everything is read straight
    from the file with h5py so the evaluator is independent of the builder.
    """
    path = Path(path)
    datasets: dict[str, np.ndarray] = {}
    attrs: dict = {}
    with h5py.File(path, "r") as fh:
        for k, v in fh.attrs.items():
            attrs[k] = v.decode() if isinstance(v, bytes) else v
        order = int(attrs["interp_order"])
        needed = ["f", "f_xi", "f_tau", "f_xitau"]
        if order == 5:
            needed += [
                "f_xixi",
                "f_tautau",
                "f_xixitau",
                "f_xitautau",
                "f_xixitautau",
            ]
        for name in needed:
            datasets[name] = np.asarray(fh[name][...], dtype=float)
        if "valid" in fh:
            datasets["valid"] = np.asarray(fh["valid"][...])
    return HermiteTable(
        order=order,
        xi_min=float(attrs["xi_min"]),
        dxi=float(attrs["dxi"]),
        nxi=int(attrs["nxi"]),
        tau_min=float(attrs["tau_min"]),
        dtau=float(attrs["dtau"]),
        ntau=int(attrs["ntau"]),
        datasets=datasets,
        attrs=attrs,
    )


def find_cell(table: HermiteTable, xi, tau):
    """Locate the cell containing each ``(xi, tau)`` and the local coordinates.

    Returns ``(ixi, itau, tx, ut)`` with 0-based lower-corner indices clamped to
    ``[0, n-2]`` (same convention as the Fortran ``find_cell``) and local
    coordinates ``tx, ut`` in ``[0, 1]`` for interior points.
    """
    xi = np.asarray(xi, dtype=float)
    tau = np.asarray(tau, dtype=float)
    ixi = np.floor((xi - table.xi_min) / table.dxi).astype(np.intp)
    ixi = np.clip(ixi, 0, table.nxi - 2)
    itau = np.floor((tau - table.tau_min) / table.dtau).astype(np.intp)
    itau = np.clip(itau, 0, table.ntau - 2)
    tx = (xi - (table.xi_min + ixi * table.dxi)) / table.dxi
    ut = (tau - (table.tau_min + itau * table.dtau)) / table.dtau
    return ixi, itau, tx, ut


def _gather_nodal(table: HermiteTable, ixi, itau, field_map, nslot):
    """Build the ``(N, nslot, nslot)`` nodal matrix D[a,b] for every query cell.

    Slot index is ``nord*corner + order`` where ``nord`` is the number of
    derivative orders per direction (2 cubic, 3 quintic), ``corner in {0,1}``
    is left/right along that axis, and ``order`` is the derivative order.
    """
    n = ixi.shape[0]
    nord = nslot // 2
    D = np.empty((n, nslot, nslot), dtype=float)
    for (p, q), name in field_map.items():
        arr = table.datasets[name]
        for ci in (0, 1):
            a = nord * ci + p
            ii = ixi + ci
            for cj in (0, 1):
                b = nord * cj + q
                jj = itau + cj
                D[:, a, b] = arr[jj, ii]
    return D


def evaluate(table: HermiteTable, xi, tau) -> dict[str, np.ndarray]:
    """Evaluate the Hermite interpolant and its derivatives at ``(xi, tau)``.

    ``xi`` and ``tau`` are broadcast together and flattened internally; the
    returned arrays have the broadcast shape.  Returns a dict with ``f``,
    ``f_xi``, ``f_tau``, ``f_xixi``, ``f_xitau``, ``f_tautau`` (all analytic
    derivatives of the same interpolant surface) plus ``ixi``, ``itau`` (the
    0-based lower-corner cell indices, for validity gating).

    For an order-3 table the second derivatives are the analytic second
    derivatives of the *bicubic* surface (built from f, f_xi, f_tau, f_xitau
    only); the stored f_xixi / f_tautau datasets are not used.
    """
    xi_b, tau_b = np.broadcast_arrays(np.asarray(xi, float), np.asarray(tau, float))
    shape = xi_b.shape
    xi_f = xi_b.ravel()
    tau_f = tau_b.ravel()

    ixi, itau, tx, ut = find_cell(table, xi_f, tau_f)

    if table.order == 5:
        wx0, wx1, wx2 = _quintic_weights(tx, table.dxi)
        wt0, wt1, wt2 = _quintic_weights(ut, table.dtau)
        field_map, nslot = _FIELD_Q, 6
    else:
        wx0, wx1, wx2 = _cubic_weights(tx, table.dxi)
        wt0, wt1, wt2 = _cubic_weights(ut, table.dtau)
        field_map, nslot = _FIELD_C, 4

    D = _gather_nodal(table, ixi, itau, field_map, nslot)

    def contract(wx, wt):
        # sum_{a,b} wx[:,a] D[:,a,b] wt[:,b]
        return np.einsum("na,nab,nb->n", wx, D, wt, optimize=True)

    out = {
        "f": contract(wx0, wt0),
        "f_xi": contract(wx1, wt0),
        "f_tau": contract(wx0, wt1),
        "f_xixi": contract(wx2, wt0),
        "f_xitau": contract(wx1, wt1),
        "f_tautau": contract(wx0, wt2),
    }
    out = {k: v.reshape(shape) for k, v in out.items()}
    out["ixi"] = ixi.reshape(shape)
    out["itau"] = itau.reshape(shape)
    return out


def cell_corner_twist_integral(
    table: HermiteTable, ixi: int, itau: int, ngauss: int = 8
):
    """Integrate the interpolant's ``f_xitau`` over cell ``(ixi, itau)`` by
    tensor Gauss-Legendre quadrature, and return ``(integral, corner_diff)``
    where ``corner_diff`` is the exact value the integral must equal:

        f(xi_{i+1},tau_{j+1}) - f(xi_i,tau_{j+1})
        - f(xi_{i+1},tau_j) + f(xi_i,tau_j).

    This is the fundamental-theorem-of-calculus identity for the mixed partial
    of the interpolant surface and holds to quadrature/roundoff for any nodal
    data.  (``ngauss`` points per direction integrate the degree<=10 integrand
    of the quintic scheme exactly up to roundoff.)
    """
    gx, gw = np.polynomial.legendre.leggauss(ngauss)
    # map nodes from [-1,1] to this cell
    xi0 = table.xi_min + ixi * table.dxi
    tau0 = table.tau_min + itau * table.dtau
    xq = xi0 + 0.5 * table.dxi * (gx + 1.0)
    tq = tau0 + 0.5 * table.dtau * (gx + 1.0)
    XX, TT = np.meshgrid(xq, tq, indexing="ij")
    vals = evaluate(table, XX.ravel(), TT.ravel())["f_xitau"].reshape(XX.shape)
    wx = 0.5 * table.dxi * gw
    wt = 0.5 * table.dtau * gw
    integral = float(np.einsum("i,ij,j->", wx, vals, wt))

    f = table.datasets["f"]
    corner_diff = float(
        f[itau + 1, ixi + 1] - f[itau + 1, ixi] - f[itau, ixi + 1] + f[itau, ixi]
    )
    return integral, corner_diff


__all__ = [
    "HermiteTable",
    "cell_corner_twist_integral",
    "evaluate",
    "find_cell",
    "load_table",
    "table_from_arrays",
]
