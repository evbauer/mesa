"""Analytic test EOS for a pure isotope (Z, A): fully ionized classical ideal
ions plus ideal nondegenerate electrons (spin degeneracy 2), with optional
blackbody radiation.

All thermodynamic relations are written in the log coordinates used by the
builder, xi = ln(rho), tau = ln(T).  The point of this module is to provide
*exact* nodal data (f and every log-coordinate derivative the processed format
needs, through fourth order, plus the two pure third derivatives) so the finite
-difference machinery in :mod:`wdm_av_builder.derivs` can be checked against a
closed-form answer.

Free energy per unit volume (task statement)::

    F/V = -kT n_i [ln(n_Q(m_i)/n_i) + 1]
          -kT n_e [ln(2 n_Q(m_e)/n_e) + 1]

with n_i = rho/(A m_u), n_e = Z n_i, m_i = A m_u, and the quantum
concentration n_Q(m) = (2 pi m k T / h^2)^{3/2}.  Dividing by rho and writing
a_i = n_i/rho = 1/(A m_u), a_e = n_e/rho = Z/(A m_u), the specific free energy
collapses to the simple form

    f(xi, tau) = -k e^{tau} [ C0 + C1 tau - C2 xi ],

    C0 = a_i b_i + a_e b_e,   C1 = (3/2)(a_i + a_e),   C2 = a_i + a_e,
    b_i = ln(C_{m_i}) - ln(a_i) + 1,   b_e = ln(2 C_{m_e}) - ln(a_e) + 1,
    C_m = (2 pi m k / h^2)^{3/2}.

Because f is (linear in xi, and e^{tau} times a linear function of tau) its
derivatives are elementary.  In particular every xi-derivative of order >= 2
vanishes for the gas part, so a nontrivial convergence test of the FD builder
along xi requires the radiation term.

Radiation adds f_rad = -a T^4/(3 rho) = -(a/3) e^{4 tau - xi}, whose mixed
derivative of order (m in xi, n in tau) is simply (-1)^m 4^n f_rad.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import constants as C

# Z -> element symbol, enough for the project isotopes and a bit beyond.
ELEMENTS: dict[int, str] = {
    1: "H",
    2: "He",
    3: "Li",
    4: "Be",
    5: "B",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    10: "Ne",
    11: "Na",
    12: "Mg",
    13: "Al",
    14: "Si",
    15: "P",
    16: "S",
    17: "Cl",
    18: "Ar",
    19: "K",
    20: "Ca",
    26: "Fe",
    28: "Ni",
}


@dataclass(frozen=True)
class Isotope:
    """A tabulated isotope: MESA name, nuclear charge, atomic mass (amu)."""

    name: str
    Z: int
    A: float

    @property
    def element(self) -> str:
        return ELEMENTS.get(self.Z, f"Z{self.Z}")


# The project isotope set (task "sample data set").
PROJECT_ISOTOPES: tuple[Isotope, ...] = (
    Isotope("h1", 1, 1.00782503),
    Isotope("he4", 2, 4.00260325),
    Isotope("c12", 6, 12.0),
    Isotope("o16", 8, 15.99491462),
    Isotope("ne20", 10, 19.99244018),
    Isotope("ne22", 10, 21.99138511),
    Isotope("mg24", 12, 23.98504170),
)

# derivative keys the processed format needs, plus the two pure third
# derivatives used for cross-checking the FD stencils.
DERIV_KEYS: tuple[str, ...] = (
    "f",
    "f_xi",
    "f_tau",
    "f_xixi",
    "f_tautau",
    "f_xitau",
    "f_xixitau",
    "f_xitautau",
    "f_xixitautau",
    "f_xixixi",
    "f_tautautau",
)


def quantum_concentration_coeff(m: float) -> float:
    """C_m = (2 pi m k / h^2)^{3/2} so that n_Q(m) = C_m T^{3/2}."""
    return (2.0 * C.pi * m * C.kerg / (C.planck_h * C.planck_h)) ** 1.5


def _gas_coeffs(Z: int, A: float) -> tuple[float, float, float]:
    """Return (C0, C1, C2) for the specific gas free energy."""
    a_i = 1.0 / (A * C.amu)  # ions per gram
    a_e = Z * a_i  # electrons per gram
    m_i = A * C.amu
    c_mi = quantum_concentration_coeff(m_i)
    c_me = quantum_concentration_coeff(C.me)
    b_i = np.log(c_mi) - np.log(a_i) + 1.0
    b_e = np.log(2.0 * c_me) - np.log(a_e) + 1.0
    c0 = a_i * b_i + a_e * b_e
    c1 = 1.5 * (a_i + a_e)
    c2 = a_i + a_e
    return c0, c1, c2


def free_energy(
    xi: np.ndarray,
    tau: np.ndarray,
    Z: int,
    A: float,
    *,
    radiation: bool = False,
) -> dict[str, np.ndarray]:
    """Specific free energy and all needed log-coordinate derivatives.

    Parameters
    ----------
    xi, tau : array_like
        Natural logs of density and temperature (broadcast together).
    Z, A : isotope charge and atomic mass (amu).
    radiation : add the blackbody term if True.

    Returns
    -------
    dict with the keys in :data:`DERIV_KEYS` plus the thermodynamic quantities
    ``P, e, s, cV, chi_rho, chi_T`` (all per gram / cgs), each a numpy array
    broadcast to the common shape of ``xi`` and ``tau``.
    """
    xi = np.asarray(xi, dtype=float)
    tau = np.asarray(tau, dtype=float)
    xi, tau = np.broadcast_arrays(xi, tau)
    k = C.kerg
    T = np.exp(tau)
    et = T  # e^{tau}

    c0, c1, c2 = _gas_coeffs(Z, A)
    u = c0 + c1 * tau - c2 * xi

    d: dict[str, np.ndarray] = {}
    # gas part
    d["f"] = -k * et * u
    d["f_xi"] = np.broadcast_to(k * c2 * et, xi.shape).copy()
    d["f_tau"] = -k * et * (u + c1)
    d["f_xixi"] = np.zeros_like(xi)
    d["f_tautau"] = -k * et * (u + 2.0 * c1)
    d["f_xitau"] = k * c2 * et
    d["f_xixitau"] = np.zeros_like(xi)
    d["f_xitautau"] = k * c2 * et
    d["f_xixitautau"] = np.zeros_like(xi)
    d["f_xixixi"] = np.zeros_like(xi)
    d["f_tautautau"] = -k * et * (u + 3.0 * c1)

    if radiation:
        f_rad = -(C.crad / 3.0) * np.exp(4.0 * tau - xi)
        # d^m/dxi^m d^n/dtau^n f_rad = (-1)^m 4^n f_rad
        d["f"] += f_rad
        d["f_xi"] += -f_rad
        d["f_tau"] += 4.0 * f_rad
        d["f_xixi"] += f_rad
        d["f_tautau"] += 16.0 * f_rad
        d["f_xitau"] += -4.0 * f_rad
        d["f_xixitau"] += 4.0 * f_rad
        d["f_xitautau"] += -16.0 * f_rad
        d["f_xixitautau"] += 16.0 * f_rad
        d["f_xixixi"] += -f_rad
        d["f_tautautau"] += 64.0 * f_rad

    # thermodynamics, assembled consistently from f and its derivatives so that
    # radiation is automatically included.
    rho = np.exp(xi)
    d["P"] = rho * d["f_xi"]
    d["e"] = d["f"] - d["f_tau"]
    d["s"] = -d["f_tau"] / T
    d["cV"] = (d["f_tau"] - d["f_tautau"]) / T
    d["chi_rho"] = 1.0 + d["f_xixi"] / d["f_xi"]
    d["chi_T"] = d["f_xitau"] / d["f_xi"]
    return d


def free_energy_direct(
    xi: np.ndarray, tau: np.ndarray, Z: int, A: float, *, radiation: bool = False
) -> dict[str, np.ndarray]:
    """P, e, s from the task's closed forms (not via f-derivatives).

    Used only as an independent cross-check in the tests.
    """
    xi = np.asarray(xi, dtype=float)
    tau = np.asarray(tau, dtype=float)
    xi, tau = np.broadcast_arrays(xi, tau)
    rho = np.exp(xi)
    T = np.exp(tau)
    a_i = 1.0 / (A * C.amu)
    n_i = a_i * rho
    n_e = Z * n_i
    P = (n_i + n_e) * C.kerg * T
    e = 1.5 * (n_i + n_e) * C.kerg * T / rho
    # f per gram from F/V divided by rho
    c_mi = quantum_concentration_coeff(A * C.amu)
    c_me = quantum_concentration_coeff(C.me)
    nq_i = c_mi * T**1.5
    nq_e = c_me * T**1.5
    fV = -C.kerg * T * n_i * (np.log(nq_i / n_i) + 1.0) - C.kerg * T * n_e * (
        np.log(2.0 * nq_e / n_e) + 1.0
    )
    f = fV / rho
    if radiation:
        P = P + C.crad * T**4 / 3.0
        e = e + C.crad * T**4 / rho
        f = f - C.crad * T**4 / (3.0 * rho)
    s = (e - f) / T
    return {"P": P, "e": e, "f": f, "s": s}


def make_source_columns(
    logrho: np.ndarray,
    logt: np.ndarray,
    Z: int,
    A: float,
    *,
    radiation: bool = False,
    with_second_derivs: bool = False,
    with_zbar: bool = False,
) -> dict[str, np.ndarray]:
    """Build 2-D source columns on the (logT outer, logRho inner) grid.

    Returns a dict of arrays shaped ``(nlogT, nlogRho)`` for ``P, E, F, S`` and,
    optionally, ``Zbar`` and the second-derivative columns
    ``dPdRho_T, dPdT_Rho, Cv``.  Also returns the exact nodal derivative dict
    under key ``_exact`` for test use (same shape per entry).
    """
    logrho = np.asarray(logrho, dtype=float)
    logt = np.asarray(logt, dtype=float)
    xi = C.ln10 * logrho[None, :]
    tau = C.ln10 * logt[:, None]
    d = free_energy(xi, tau, Z, A, radiation=radiation)
    rho = np.exp(xi)
    T = np.exp(tau)
    cols: dict[str, np.ndarray] = {
        "P": np.broadcast_to(d["P"], (logt.size, logrho.size)).copy(),
        "E": np.broadcast_to(d["e"], (logt.size, logrho.size)).copy(),
        "F": np.broadcast_to(d["f"], (logt.size, logrho.size)).copy(),
        "S": np.broadcast_to(d["s"], (logt.size, logrho.size)).copy(),
    }
    if with_zbar:
        cols["Zbar"] = np.full((logt.size, logrho.size), float(Z))
    if with_second_derivs:
        # (dP/drho)_T = f_xi + f_xixi ;  (dP/dT)_rho = rho f_xitau / T
        cols["dPdRho_T"] = np.broadcast_to(
            d["f_xi"] + d["f_xixi"], (logt.size, logrho.size)
        ).copy()
        cols["dPdT_Rho"] = np.broadcast_to(
            rho * d["f_xitau"] / T, (logt.size, logrho.size)
        ).copy()
        cols["Cv"] = np.broadcast_to(d["cV"], (logt.size, logrho.size)).copy()
    cols["_exact"] = d
    return cols
