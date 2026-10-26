"""Synthetic analytic EOS: exact derivatives vs numerical differentiation of
the synthetic free energy, radiation scaling identity, and the direct vs
f-based thermodynamics."""

from __future__ import annotations

import numpy as np
import pytest
from wdm_av_builder import constants as C
from wdm_av_builder import synthetic as s


@pytest.mark.parametrize("Z,A", [(1, 1.00782503), (8, 16.0), (12, 23.9850417)])
def test_gas_derivatives_vs_finite_difference(Z, A):
    """Gas-only analytic derivatives match central differences of f."""
    xi0, tau0 = np.log(1.0e-2), np.log(3.0e5)
    h = 1e-4

    def fval(x, t):
        return float(s.free_energy(x, t, Z, A)["f"])

    d = s.free_energy(xi0, tau0, Z, A)
    fx = (fval(xi0 + h, tau0) - fval(xi0 - h, tau0)) / (2 * h)
    ft = (fval(xi0, tau0 + h) - fval(xi0, tau0 - h)) / (2 * h)
    ftt = (fval(xi0, tau0 + h) - 2 * fval(xi0, tau0) + fval(xi0, tau0 - h)) / h**2
    fxt = (
        fval(xi0 + h, tau0 + h)
        - fval(xi0 + h, tau0 - h)
        - fval(xi0 - h, tau0 + h)
        + fval(xi0 - h, tau0 - h)
    ) / (4 * h**2)
    assert fx == pytest.approx(float(d["f_xi"]), rel=1e-6)
    assert ft == pytest.approx(float(d["f_tau"]), rel=1e-6)
    assert ftt == pytest.approx(float(d["f_tautau"]), rel=1e-4)
    assert fxt == pytest.approx(float(d["f_xitau"]), rel=1e-4)
    # gas f is linear in xi -> second and higher xi derivatives vanish exactly
    assert float(d["f_xixi"]) == 0.0
    assert float(d["f_xixitau"]) == 0.0
    assert float(d["f_xixitautau"]) == 0.0


def test_radiation_scaling_identity():
    """Every radiation derivative equals (-1)^m 4^n f_rad."""
    xi = np.log(np.array([1e-4, 1.0, 1e3]))
    tau = np.log(np.array([1e4, 1e5, 1e6]))
    Z, A = 6, 12.0
    with_rad = s.free_energy(xi, tau, Z, A, radiation=True)
    no_rad = s.free_energy(xi, tau, Z, A, radiation=False)
    f_rad = -(C.crad / 3.0) * np.exp(4.0 * tau - xi)
    # (m in xi, n in tau)
    expect = {
        "f": (0, 0),
        "f_xi": (1, 0),
        "f_tau": (0, 1),
        "f_xixi": (2, 0),
        "f_tautau": (0, 2),
        "f_xitau": (1, 1),
        "f_xixitau": (2, 1),
        "f_xitautau": (1, 2),
        "f_xixitautau": (2, 2),
        "f_xixixi": (3, 0),
        "f_tautautau": (0, 3),
    }
    for key, (m, n) in expect.items():
        contrib = with_rad[key] - no_rad[key]
        target = (-1.0) ** m * 4.0**n * f_rad
        # for keys with a large gas part the subtraction loses precision;
        # rtol 1e-6 still rules out any wrong integer sign/factor.
        assert np.allclose(contrib, target, rtol=1e-6, atol=0.0)


@pytest.mark.parametrize("radiation", [False, True])
def test_direct_vs_fbased_thermo(radiation):
    """P, e, s from the task's closed forms match the f-derivative route."""
    xi = np.log(np.array([1e-3, 1.0]))
    tau = np.log(np.array([1e4, 1e6]))
    d = s.free_energy(xi, tau, 8, 16.0, radiation=radiation)
    dd = s.free_energy_direct(xi, tau, 8, 16.0, radiation=radiation)
    assert np.allclose(d["P"], dd["P"], rtol=1e-12)
    assert np.allclose(d["e"], dd["e"], rtol=1e-12)
    assert np.allclose(d["s"], dd["s"], rtol=1e-12)
    assert np.allclose(d["f"], dd["f"], rtol=1e-12)


def test_first_law_relations():
    """f_xi = P/rho, f_tau = f - e, c_V = (f_tau - f_tautau)/T hold."""
    xi = np.log(np.array([1e-2, 10.0]))
    tau = np.log(np.array([1e4, 1e7]))
    d = s.free_energy(xi, tau, 2, 4.0026, radiation=True)
    rho = np.exp(xi)
    T = np.exp(tau)
    assert np.allclose(d["f_xi"], d["P"] / rho, rtol=1e-12)
    assert np.allclose(d["f_tau"], d["f"] - d["e"], rtol=1e-12)
    assert np.allclose(d["cV"], (d["f_tau"] - d["f_tautau"]) / T, rtol=1e-12)
