"""Convergence of the FD-estimated nodal derivatives toward the analytic
synthetic values as the grid is refined (two resolutions)."""

from __future__ import annotations

import numpy as np
from wdm_av_builder import io as wio
from wdm_av_builder.derivs import build_nodal_data
from wdm_av_builder.synthetic import make_source_columns


def _build(ppd, radiation, lr, lt, Z=8, A=16.0):
    dl = 1.0 / ppd
    logrho = np.arange(lr[0], lr[1] + 1e-9, dl)
    logt = np.arange(lt[0], lt[1] + 1e-9, dl)
    cols = make_source_columns(logrho, logt, Z, A, radiation=radiation)
    exact = cols.pop("_exact")
    src = wio.SourceTable(
        header={"isotope": "o16", "Z": Z, "A": A, "includes_radiation": False},
        logrho=logrho,
        logt=logt,
        columns=cols,
    )
    return build_nodal_data(src, order=5), exact


def _rel_err(nd, exact, key, pad=8):
    sl = (slice(pad, -pad), slice(pad, -pad))
    num = np.max(np.abs(nd.datasets[key][sl] - exact[key][sl]))
    den = np.max(np.abs(exact[key][sl])) + 1e-300
    return num / den


def test_tau_derivatives_converge_fourth_order():
    """f_tautau, f_xitau converge at >= 4th order; f_xitautau faster (7-pt)."""
    lr, lt = (-2.0, 2.0), (4.0, 6.0)
    a = _build(10, False, lr, lt)
    b = _build(20, False, lr, lt)
    for key, floor in [("f_tautau", 3.5), ("f_xitau", 3.5), ("f_xitautau", 5.0)]:
        e1 = _rel_err(*a, key)
        e2 = _rel_err(*b, key)
        order = np.log2(e1 / e2)
        assert order > floor, f"{key} order {order:.2f} < {floor}"


def test_xi_derivatives_converge_fourth_order():
    """With radiation (nontrivial along xi) the xi second derivatives converge."""
    lr, lt = (-3.0, 3.0), (3.0, 5.0)
    a = _build(10, True, lr, lt)
    b = _build(20, True, lr, lt)
    for key, floor in [
        ("f_xixi", 3.5),
        ("f_xixitau", 3.5),
        ("f_xixitautau", 3.5),
    ]:
        e1 = _rel_err(*a, key)
        e2 = _rel_err(*b, key)
        order = np.log2(e1 / e2)
        assert order > floor, f"{key} order {order:.2f} < {floor}"


def test_estimated_derivatives_accurate_at_fine_grid():
    """Sanity: at 20 pts/decade interior relative errors are small."""
    nd, exact = _build(20, True, (-3.0, 3.0), (3.0, 5.0))
    for key in ("f_xixi", "f_tautau", "f_xitau"):
        assert _rel_err(nd, exact, key) < 1e-2
