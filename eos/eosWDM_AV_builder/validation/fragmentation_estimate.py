"""Electron-fragmentation error of additive-volume mixing, isolated.

Model: nonrelativistic ideal Fermi-Dirac electrons plus ideal ions, no Coulomb
term, no radiation.  At common (P, T) each species' electron gas sits at its own
n_e,i because the ion pressure per electron (kT/Z_i) differs, while the true
mixture has one electron gas at the mean n_e.  Additive volume is exact in both
the nondegenerate limit (Dalton) and the fully degenerate limit (pressure
equilibrium equalizes the n_e,i); in between it is not, and this script gives
the size of that error alone, for comparison with the WDM_AV-vs-reference maps
(validation/mixtures.py).

Usage: python3 fragmentation_estimate.py
"""

import numpy as np
from scipy.optimize import brentq

kB = 1.380649e-16
h = 6.62607015e-27
me = 9.1093837e-28
mu = 1.66053907e-24

_x = np.linspace(0.0, 400.0, 400001)


def fd(j: float, eta: float) -> float:
    f = _x**j / (np.exp(np.clip(_x - eta, -700, 700)) + 1.0)
    return float(np.trapezoid(f, _x))


def ne_of_eta(eta: float, T: float) -> float:
    return 4 * np.pi / h**3 * (2 * me * kB * T) ** 1.5 * fd(0.5, eta)


def pe_of_eta(eta: float, T: float) -> float:
    n = ne_of_eta(eta, T)
    return (2.0 / 3.0) * n * kB * T * fd(1.5, eta) / fd(0.5, eta)


def ee_of_eta(eta: float, T: float) -> float:
    """Electron kinetic energy per volume."""
    return 1.5 * pe_of_eta(eta, T)


def eta_of_ne(n: float, T: float) -> float:
    return brentq(lambda e: ne_of_eta(e, T) - n, -60, 400)


def species_state(P: float, T: float, Z: float) -> tuple[float, float]:
    """(eta_i, n_e,i) such that P_e(n) + (n/Z) kT = P (one ion per Z electrons)."""

    def g(e: float) -> float:
        return pe_of_eta(e, T) + ne_of_eta(e, T) / Z * kB * T - P

    e = brentq(g, -60, 400)
    return e, ne_of_eta(e, T)


def av_mixture(rho: float, T: float, X: list, Z: list, A: list) -> tuple:
    """Return (dP/P, de/e, eta_true, [n_e,i / n_e]) of additive volume vs exact."""
    ne_true = rho / mu * sum(Xi * Zi / Ai for Xi, Zi, Ai in zip(X, Z, A))
    nion_true = rho / mu * sum(Xi / Ai for Xi, Ai in zip(X, A))
    eta_t = eta_of_ne(ne_true, T)
    P_true = pe_of_eta(eta_t, T) + nion_true * kB * T
    e_true = (ee_of_eta(eta_t, T) + 1.5 * nion_true * kB * T) / rho

    def vol_resid(lnP: float) -> float:
        P = np.exp(lnP)
        v = 0.0
        for Xi, Zi, Ai in zip(X, Z, A):
            _, n = species_state(P, T, Zi)
            v += Xi / (n * mu * Ai / Zi)
        return v - 1.0 / rho

    lnP = brentq(vol_resid, np.log(P_true) - 1.0, np.log(P_true) + 1.0)
    P_av = np.exp(lnP)
    e_av = 0.0
    split = []
    for Xi, Zi, Ai in zip(X, Z, A):
        eta_i, n = species_state(P_av, T, Zi)
        rho_i = n * mu * Ai / Zi
        E_i = ee_of_eta(eta_i, T) + n / Zi * 1.5 * kB * T
        e_av += Xi * E_i / rho_i
        split.append(n / ne_true)
    return P_av / P_true - 1.0, e_av / e_true - 1.0, eta_t, split


MIXES = {
    # name: (X mass fractions, Z, A); stems match make_mesa_tables / mixtures.py
    "C/O 50/50": ([0.5, 0.5], [6, 8], [12.0, 16.0]),
    "O/Ne/Mg 55/35/10": ([0.55, 0.35, 0.10], [8, 10, 12], [16.0, 20.0, 24.0]),
    "He/C 50/50": ([0.5, 0.5], [2, 6], [4.0, 12.0]),
    "H/He 70/30": ([0.70, 0.30], [1, 2], [1.0, 4.0]),
    "C/O/22Ne 49/49/2": ([0.49, 0.49, 0.02], [6, 8, 10], [12.0, 16.0, 22.0]),
    "C/O 90/10": ([0.9, 0.1], [6, 8], [12.0, 16.0]),
    "C/O 99/1": ([0.99, 0.01], [6, 8], [12.0, 16.0]),
    "He/C 99/1": ([0.99, 0.01], [2, 6], [4.0, 12.0]),
    "He/C 10/90": ([0.10, 0.90], [2, 6], [4.0, 12.0]),
    "C/O/Ne 45/45/10": ([0.45, 0.45, 0.1], [6, 8, 10], [12.0, 16.0, 20.0]),
    "He/C/O 20/40/40": ([0.2, 0.4, 0.4], [2, 6, 8], [4.0, 12.0, 16.0]),
}
POINTS = [(0.0, 6.0), (1.0, 6.5), (2.0, 6.5), (2.0, 7.0), (3.0, 6.0), (4.0, 7.0)]

# target-gap grid for the --grid median (logRho in [-2,4], logT in [3,7])
GRID_LOGRHO = np.arange(-2.0, 4.01, 1.0)
GRID_LOGT = np.arange(3.0, 7.01, 1.0)


def grid_medians(X, Z, A):
    """Median |dP/P| and |de/e| over the target-gap grid (failed points skipped)."""
    dps, des = [], []
    for lr in GRID_LOGRHO:
        for lt in GRID_LOGT:
            try:
                dP, dE, _eta, _split = av_mixture(10**lr, 10**lt, X, Z, A)
            except (ValueError, RuntimeError):
                continue
            if np.isfinite(dP) and np.isfinite(dE):
                dps.append(abs(dP))
                des.append(abs(dE))
    if not dps:
        return None, None, 0
    return float(np.median(dps)), float(np.median(des)), len(dps)


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--grid",
        action="store_true",
        help="also report the median |dP/P| and |de/e| over the target-gap grid",
    )
    args = ap.parse_args()

    print(
        f"{'mixture':18s} {'logRho':>6s} {'logT':>5s} {'eta':>7s} "
        f"{'dP/P':>10s} {'de/e':>10s}  n_e,i/n_e"
    )
    for name, (X, Z, A) in MIXES.items():
        for lr, lt in POINTS:
            dP, dE, eta, split = av_mixture(10**lr, 10**lt, X, Z, A)
            print(
                f"{name:18s} {lr:6.1f} {lt:5.1f} {eta:7.2f} {dP:10.2e} {dE:10.2e}  "
                + " ".join(f"{s:.3f}" for s in split)
            )

    if args.grid:
        print()
        print(
            f"{'mixture':18s} {'med|dP/P|':>11s} {'med|de/e|':>11s} "
            f"{'n_pts':>6s}  (target gap: logRho[-2,4] x logT[3,7], 1-dex grid)"
        )
        for name, (X, Z, A) in MIXES.items():
            mdp, mde, npts = grid_medians(X, Z, A)
            if mdp is None:
                print(f"{name:18s} {'-':>11s} {'-':>11s} {0:6d}")
            else:
                print(f"{name:18s} {mdp:11.2e} {mde:11.2e} {npts:6d}")


if __name__ == "__main__":
    main()
