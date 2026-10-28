"""Tests for the blend-zone analysis helpers in ``blend.py``.

Synthetic only (no data files): the finite-difference operator, the ramp-centre
detection, and the one-sided trend fit used to measure the value step.
"""

from __future__ import annotations

import blend as bl
import numpy as np


def test_central_fd_linear():
    # value = m*coord + b  =>  d(value)/d(ln x) = m / ln10  (x = 10**coord)
    coord = np.linspace(3.0, 6.0, 61)
    m, b = 2.5, -1.0
    val = m * coord + b
    fd = bl.central_fd(val, coord)
    expected = m / bl.LN10
    assert np.allclose(fd[1:-1], expected, rtol=1e-10)
    assert np.isnan(fd[0]) and np.isnan(fd[-1])


def test_central_fd_quadratic_is_truncation_accurate():
    # for a smooth function central differences are accurate to O(step^2)
    coord = np.linspace(0.0, 2.0, 401)
    val = np.exp(coord)  # d/d(ln x) = exp(coord)/ln10
    fd = bl.central_fd(val, coord)
    exact = np.exp(coord) / bl.LN10
    rel = np.abs(fd[1:-1] - exact[1:-1]) / np.abs(exact[1:-1])
    assert rel.max() < 1e-4


def test_ramp_center_descending():
    coord = np.linspace(4.0, 5.0, 101)
    # frac_WDM goes 1 -> 0 linearly, crossing 0.5 at coord = 4.5
    frac = 1.0 - (coord - 4.0)
    c = bl.ramp_center(coord, frac)
    assert c is not None
    assert abs(c - 4.5) < 1e-6


def test_ramp_center_none_when_no_crossing():
    coord = np.linspace(0.0, 1.0, 11)
    frac = np.ones_like(coord)  # fully WDM_AV, never crosses 0.5
    assert bl.ramp_center(coord, frac) is None


def test_side_fit_extrapolates_line():
    coord = np.linspace(0.0, 1.0, 101)
    val = 3.0 * coord + 2.0
    mask = coord <= 0.4  # only the left side is "pure"
    # extrapolate the left-side line to center = 0.5 -> 3*0.5 + 2 = 3.5
    got = bl._side_fit(coord, val, mask, center=0.5, nfit=40, deg=2)
    assert got is not None
    assert abs(got - 3.5) < 1e-6


def test_side_fit_too_few_points():
    coord = np.array([0.0, 1.0])
    val = np.array([0.0, 1.0])
    mask = np.array([True, False])  # a single usable point, deg+1 = 3 required
    assert bl._side_fit(coord, val, mask, center=0.5) is None


def test_partial_name_matches_axis():
    assert bl.partial_name("rho_right_T6.0", "lnPgas") == "dlnPgas_dlnRho"
    assert bl.partial_name("T_top_R4.0", "lnE") == "dlnE_dlnT"
