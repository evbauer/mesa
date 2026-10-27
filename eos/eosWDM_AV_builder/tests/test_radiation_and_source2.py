"""Radiation stripping recovers the no-radiation nodal data, and the
source-supplied second-derivative path matches the analytic values."""

from __future__ import annotations

import numpy as np
from wdm_av_builder import io as wio
from wdm_av_builder.derivs import build_nodal_data
from wdm_av_builder.synthetic import make_source_columns


def _source(radiation, second=False, Z=8, A=16.0):
    logrho = np.linspace(-4.0, 2.0, 61)
    logt = np.linspace(4.0, 6.5, 51)
    cols = make_source_columns(
        logrho, logt, Z, A, radiation=radiation, with_second_derivs=second
    )
    cols.pop("_exact")
    return wio.SourceTable(
        header={"isotope": "o16", "Z": Z, "A": A, "includes_radiation": radiation},
        logrho=logrho,
        logt=logt,
        columns=cols,
    )


def test_radiation_stripping_recovers_norad():
    """Building from a table that *includes* radiation, after analytic
    stripping, reproduces the table built without radiation."""
    with_rad = build_nodal_data(_source(True), order=5)
    no_rad = build_nodal_data(_source(False), order=5)
    # global scale for absolute tolerance (some fields are identically zero for
    # the ideal gas, e.g. f_xixi, so a per-field relative test is not enough)
    scale = np.abs(no_rad.datasets["f"]).max()
    for key in ("f", "f_xi", "f_tau", "f_xixi", "f_tautau", "f_xitau"):
        a = with_rad.datasets[key]
        b = no_rad.datasets[key]
        assert np.allclose(a, b, rtol=1e-8, atol=1e-6 * scale), key


def test_second_derivs_from_source_exact():
    """When dPdRho_T, dPdT_Rho, Cv are supplied, f_xixi/f_xitau/f_tautau equal
    the analytic values (no FD error), and the flag is set."""
    src = _source(False, second=True)
    exact = make_source_columns(src.logrho, src.logt, 8, 16.0, with_second_derivs=True)[
        "_exact"
    ]
    nd = build_nodal_data(src, order=5)
    assert nd.attrs["second_derivs_from_source"] == 1
    # f_xixi = (dP/drho)_T - f_xi : for the ideal gas this is identically 0,
    # computed as a difference of ~1e14 numbers so only ~roundoff survives.
    assert np.max(np.abs(nd.datasets["f_xixi"])) < 1e-10 * np.max(
        np.abs(nd.datasets["f_xi"])
    )
    assert np.allclose(nd.datasets["f_xitau"], exact["f_xitau"], rtol=1e-12)
    assert np.allclose(nd.datasets["f_tautau"], exact["f_tautau"], rtol=1e-12)


def test_ignore_source_derivs_forces_fd():
    """``ignore_source_derivs`` takes the finite-difference path even when the
    dPdRho_T/dPdT_Rho/Cv columns are present: the flag is cleared, the scheme is
    annotated, and the second derivatives match a source with those columns
    dropped entirely (pure FD)."""
    src = _source(False, second=True)
    nd_src = build_nodal_data(src, order=5)
    nd_fd = build_nodal_data(src, order=5, ignore_source_derivs=True)
    assert nd_src.attrs["second_derivs_from_source"] == 1
    assert nd_fd.attrs["second_derivs_from_source"] == 0
    assert "forced_fd_ignore_source_derivs" in nd_fd.attrs["fd_scheme"]
    # same result as building from a source lacking the second-derivative columns
    src_noder = _source(False, second=False)
    nd_noder = build_nodal_data(src_noder, order=5)
    for key in ("f_xixi", "f_tautau", "f_xitau", "f_xixitau", "f_xitautau"):
        assert np.allclose(nd_fd.datasets[key], nd_noder.datasets[key]), key


def test_second_derivs_fxixi_relation_no_rho_factor():
    """Guard the corrected relation f_xixi = (dP/drho)_T - f_xi (no rho factor):
    for the ideal gas (dP/drho)_T == f_xi so f_xixi must be ~0, which the
    rho-factor version would badly violate."""
    src = _source(False, second=True)
    nd = build_nodal_data(src, order=5)
    f_xi = nd.datasets["f_xi"]
    # the (wrong) rho-factor form would give rho*dPdRho_T - f_xi ~ rho*f_xi,
    # i.e. enormous; the correct form is ~0.
    assert np.max(np.abs(nd.datasets["f_xixi"])) < 1e-3 * np.max(np.abs(f_xi))
