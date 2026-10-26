"""Finite-difference machinery: Fornberg weights and apply_fd exact on
polynomials, and matching the classic 4th-order stencils."""

from __future__ import annotations

import math

import numpy as np
import pytest
from wdm_av_builder.derivs import apply_fd, fd_weights, fornberg_weights


@pytest.mark.parametrize("deriv", [1, 2, 3])
@pytest.mark.parametrize("npts", [5, 7, 9])
def test_fd_weights_exact_on_monomials(deriv, npts):
    """Weights reproduce derivatives of monomials up to degree npts-1.

    The derivative of x^deg at x=0 is deriv! when deg == deriv, else 0.
    """
    offsets = np.arange(npts) - npts // 2
    w = fd_weights(offsets, deriv)  # unit spacing
    for deg in range(npts):
        vals = offsets.astype(float) ** deg
        got = w @ vals
        exact = math.factorial(deriv) if deg == deriv else 0.0
        assert got == pytest.approx(exact, abs=1e-9)


def test_classic_fourth_order_first_derivative():
    """apply_fd matches the centered 4th-order first-derivative stencil."""
    x = np.linspace(0.0, 2.0, 21)
    h = x[1] - x[0]
    y = np.sin(3.0 * x)
    got = apply_fd(y, axis=0, deriv=1, h=h)
    # interior reference
    i = 10
    ref = (-y[i + 2] + 8 * y[i + 1] - 8 * y[i - 1] + y[i - 2]) / (12 * h)
    assert got[i] == pytest.approx(ref, rel=1e-12)
    assert got[i] == pytest.approx(3.0 * np.cos(3.0 * x[i]), rel=2e-3)


def test_classic_fourth_order_second_derivative():
    x = np.linspace(0.0, 2.0, 21)
    h = x[1] - x[0]
    y = np.cos(2.0 * x)
    got = apply_fd(y, axis=0, deriv=2, h=h)
    i = 10
    assert got[i] == pytest.approx(-4.0 * np.cos(2.0 * x[i]), rel=2e-3)


def test_fornberg_shape():
    c = fornberg_weights(0.0, np.arange(5.0) - 2.0, 2)
    assert c.shape == (5, 3)
    # 0th derivative weights at a node sum to 1
    assert c[:, 0].sum() == pytest.approx(1.0)


def test_apply_fd_two_dimensional_axis():
    """Derivative along each axis of a 2-D field matches analytic."""
    xi = np.linspace(-1, 1, 25)
    tau = np.linspace(0, 2, 17)
    X, T = np.meshgrid(xi, tau)  # (ntau, nxi)
    f = np.exp(0.5 * X) * np.sin(T)
    hxi = xi[1] - xi[0]
    htau = tau[1] - tau[0]
    fxi = apply_fd(f, axis=1, deriv=1, h=hxi)
    ftau = apply_fd(f, axis=0, deriv=1, h=htau)
    assert np.allclose(fxi[2:-2, 2:-2], 0.5 * f[2:-2, 2:-2], rtol=1e-4)
    ref = np.exp(0.5 * X) * np.cos(T)
    assert np.allclose(ftau[2:-2, 2:-2], ref[2:-2, 2:-2], rtol=1e-4)
