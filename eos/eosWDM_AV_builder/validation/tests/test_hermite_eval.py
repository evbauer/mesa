"""Correctness tests for the independent Hermite evaluator.

These verify the evaluator against the *mathematical* properties it must have,
which is exactly how it validates the Fortran reader and the builder:

* at every node it reproduces the stored nodal data to roundoff
  (f, f_xi, f_tau for order 3; plus f_xixi, f_tautau, f_xitau for order 5);
* the cell-corner identity ``integral(f_xitau) == corner differences of f``
  holds for every cell;
* its analytic derivatives agree with finite differences of its own f surface.
"""

from __future__ import annotations

import hermite_eval as he
import numpy as np
import pytest
from wdm_av_builder import io as wio
from wdm_av_builder.cli import build_one
from wdm_av_builder.synthetic import Isotope, make_source_columns


def _make_processed(tmp_path, order, *, radiation_kept=True):
    """Build a processed HDF5 table from a synthetic source.

    ``radiation_kept=True`` keeps the radiation term in the table (the source is
    not flagged ``includes_radiation``), so the second and mixed derivatives are
    genuinely non-zero and the quintic basis is exercised non-trivially.
    """
    iso = Isotope("o16", 8, 15.99491462)
    logrho = np.linspace(-5.0, 3.0, 41)
    logt = np.linspace(4.0, 7.5, 33)
    cols = make_source_columns(
        logrho, logt, iso.Z, iso.A, radiation=True, with_second_derivs=True
    )
    cols.pop("_exact")
    header = {
        "isotope": iso.name,
        "Z": iso.Z,
        "A": iso.A,
        "element": iso.element,
        "includes_radiation": not radiation_kept,
        "model": "test",
    }
    src = wio.SourceTable(header=header, logrho=logrho, logt=logt, columns=cols)
    src_path = tmp_path / f"src_{order}.txt"
    wio.write_source_ascii(src, src_path)
    out = tmp_path / f"eosWDM_AV_o16_{order}.h5"
    build_one(src_path, out, order=order, write_audit=False)
    return he.load_table(out)


def _node_grid(t):
    xi = t.xi_min + np.arange(t.nxi) * t.dxi
    tau = t.tau_min + np.arange(t.ntau) * t.dtau
    return np.meshgrid(xi, tau)  # (ntau, nxi)


@pytest.mark.parametrize("order", [3, 5])
def test_nodal_reproduction(tmp_path, order):
    t = _make_processed(tmp_path, order)
    XX, TT = _node_grid(t)
    out = he.evaluate(t, XX, TT)
    keys = ["f", "f_xi", "f_tau", "f_xitau"]
    if order == 5:
        keys += ["f_xixi", "f_tautau"]
    for k in keys:
        ref = t.datasets[k]
        got = out[k]
        scale = np.abs(ref).max()
        assert np.max(np.abs(got - ref)) < 1e-11 * scale, k


def test_cell_indices_at_nodes(tmp_path):
    """Query points exactly on nodes land in a cell for which they are a corner:
    the local coordinate is at an endpoint (tx ~ 0 in cell i, or tx ~ 1 in cell
    i-1 when the exact division floors just below the integer) and the cell is
    in range.  Half-way points land strictly in the interior (tx ~ 0.5)."""
    t = _make_processed(tmp_path, 5)
    xi = t.xi_min + np.arange(t.nxi) * t.dxi
    ixi, _, tx, _ = he.find_cell(t, xi, np.full_like(xi, t.tau_min))
    assert np.all((ixi >= 0) & (ixi <= t.nxi - 2))
    at_endpoint = np.isclose(tx, 0.0, atol=1e-9) | np.isclose(tx, 1.0, atol=1e-9)
    assert np.all(at_endpoint)
    # the node is a corner of the chosen cell: reconstructed xi matches
    recon = t.xi_min + (ixi + tx) * t.dxi
    assert np.allclose(recon, xi)
    # midpoints of cells land in the interior with tx ~ 0.5
    xi_mid = t.xi_min + (np.arange(t.nxi - 1) + 0.5) * t.dxi
    ixim, _, txm, _ = he.find_cell(t, xi_mid, np.full_like(xi_mid, t.tau_min))
    assert np.all(ixim == np.arange(t.nxi - 1))
    assert np.allclose(txm, 0.5)


@pytest.mark.parametrize("order", [3, 5])
def test_cell_corner_twist_identity(tmp_path, order):
    t = _make_processed(tmp_path, order)
    f = t.datasets["f"]
    fscale = np.abs(f).max()
    rng = np.random.default_rng(0)
    cells = [
        (int(i), int(j))
        for i, j in zip(rng.integers(0, t.nxi - 1, 40), rng.integers(0, t.ntau - 1, 40))
    ]
    for ixi, itau in cells:
        integral, corner = he.cell_corner_twist_integral(t, ixi, itau)
        # absolute, scaled by the global f magnitude (corner diff can be ~0)
        assert abs(integral - corner) < 1e-9 * fscale


@pytest.mark.parametrize("order", [3, 5])
def test_fd_consistency(tmp_path, order):
    """Analytic derivatives of the interpolant match central finite differences
    of its own f surface, evaluated at cell centres (so the stencil stays inside
    one polynomial patch)."""
    t = _make_processed(tmp_path, order)
    # cell-centre query points over a block of interior cells
    ii = np.arange(5, t.nxi - 5, 3)
    jj = np.arange(5, t.ntau - 5, 3)
    II, JJ = np.meshgrid(ii, jj)
    xi = t.xi_min + (II + 0.5) * t.dxi
    tau = t.tau_min + (JJ + 0.5) * t.dtau
    ref = he.evaluate(t, xi, tau)

    dx = 1e-3 * t.dxi
    dt = 1e-3 * t.dtau

    def ev(a, b):
        return he.evaluate(t, xi + a * dx, tau + b * dt)["f"]

    f_xi = (ev(1, 0) - ev(-1, 0)) / (2 * dx)
    f_tau = (ev(0, 1) - ev(0, -1)) / (2 * dt)
    f_xixi = (ev(1, 0) - 2 * ref["f"] + ev(-1, 0)) / dx**2
    f_tautau = (ev(0, 1) - 2 * ref["f"] + ev(0, -1)) / dt**2
    f_xitau = (ev(1, 1) - ev(1, -1) - ev(-1, 1) + ev(-1, -1)) / (4 * dx * dt)

    def close(fd, key, rtol):
        a = ref[key]
        scale = np.abs(a).max()
        assert np.max(np.abs(fd - a)) < rtol * scale, key

    close(f_xi, "f_xi", 1e-6)
    close(f_tau, "f_tau", 1e-6)
    close(f_xitau, "f_xitau", 1e-5)
    close(f_xixi, "f_xixi", 1e-4)
    close(f_tautau, "f_tautau", 1e-4)


def test_biquintic_accuracy_vs_analytic(tmp_path):
    """With exact nodal data, the biquintic interpolant reproduces the smooth
    analytic free energy between nodes to high accuracy."""
    from wdm_av_builder.synthetic import free_energy

    t = _make_processed(tmp_path, 5)
    # mid-cell points, away from the edges
    ii = np.arange(5, t.nxi - 5)
    jj = np.arange(5, t.ntau - 5)
    II, JJ = np.meshgrid(ii, jj)
    xi = t.xi_min + (II + 0.5) * t.dxi
    tau = t.tau_min + (JJ + 0.5) * t.dtau
    out = he.evaluate(t, xi, tau)
    # evaluate the analytic free energy directly at the scattered points
    d = free_energy(xi, tau, 8, 15.99491462, radiation=True)
    rel = np.abs(out["f"] - d["f"]) / (np.abs(d["f"]) + 1e-300)
    assert np.median(rel) < 1e-7
    # P = rho f_xi accurate too
    relP = np.abs(out["f_xi"] - d["f_xi"]) / (np.abs(d["f_xi"]) + 1e-300)
    assert np.median(relP) < 1e-6
