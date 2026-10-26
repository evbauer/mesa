"""I/O round trips: source ASCII and HDF5, and processed-file conformance to
FORMAT.md including the guaranteed nodal relations."""

from __future__ import annotations

import numpy as np
import pytest
from wdm_av_builder import constants as C
from wdm_av_builder import io as wio
from wdm_av_builder.cli import build_one
from wdm_av_builder.derivs import build_nodal_data
from wdm_av_builder.synthetic import Isotope, make_source_columns


def _make_source(iso: Isotope, *, radiation=False, second=False, zbar=False):
    logrho = np.linspace(-6.0, 2.0, 41)
    logt = np.linspace(4.0, 8.0, 31)
    cols = make_source_columns(
        logrho,
        logt,
        iso.Z,
        iso.A,
        radiation=radiation,
        with_second_derivs=second,
        with_zbar=zbar,
    )
    cols.pop("_exact")
    header = {
        "isotope": iso.name,
        "Z": iso.Z,
        "A": iso.A,
        "element": iso.element,
        "includes_radiation": radiation,
        "model": "test",
    }
    return wio.SourceTable(header=header, logrho=logrho, logt=logt, columns=cols)


def test_source_ascii_round_trip(tmp_path):
    iso = Isotope("o16", 8, 15.99491462)
    tbl = _make_source(iso, zbar=True)
    path = tmp_path / "src.txt"
    wio.write_source_ascii(tbl, path)
    back = wio.read_source(path)
    assert back.name == "o16"
    assert back.Z == 8
    assert back.A == pytest.approx(15.99491462)
    assert np.allclose(back.logrho, tbl.logrho)
    assert np.allclose(back.logt, tbl.logt)
    for key in ("P", "E", "F", "S", "Zbar"):
        assert np.allclose(back.columns[key], tbl.columns[key], rtol=1e-12)


def test_source_hdf5_round_trip(tmp_path):
    import h5py

    iso = Isotope("c12", 6, 12.0)
    tbl = _make_source(iso)
    path = tmp_path / "src.h5"
    with h5py.File(path, "w") as f:
        f.attrs["isotope"] = iso.name
        f.attrs["Z"] = np.int32(iso.Z)
        f.attrs["A"] = np.float64(iso.A)
        f.attrs["includes_radiation"] = np.int32(0)
        f.create_dataset("logRho", data=tbl.logrho)
        f.create_dataset("logT", data=tbl.logt)
        for k, v in tbl.columns.items():
            f.create_dataset(k, data=v)
    back = wio.read_source(path)
    assert back.Z == 6
    assert np.allclose(back.columns["P"], tbl.columns["P"])


def test_processed_conforms_and_relations(tmp_path):
    iso = Isotope("ne20", 10, 19.99244018)
    src = _make_source(iso, zbar=True)
    src_path = tmp_path / "ne20.txt"
    wio.write_source_ascii(src, src_path)
    out = tmp_path / "eosWDM_AV_ne20.h5"
    build_one(src_path, out, order=5)

    assert wio.validate_processed(out) == []

    data = wio.read_processed_hdf5(out)
    attrs, dsets = data["attrs"], data["datasets"]
    # attribute types
    assert isinstance(attrs["Z"], int) and attrs["Z"] == 10
    assert isinstance(attrs["A"], float)
    assert isinstance(attrs["isotope"], str)
    assert attrs["format_version"] == 1
    assert attrs["includes_radiation"] == 0
    assert attrs["nxi"] == src.nlogrho and attrs["ntau"] == src.nlogt
    # shapes (ntau, nxi)
    for name in wio.REQUIRED_DATASETS:
        assert dsets[name].shape == (src.nlogt, src.nlogrho)
    # guaranteed relations at nodes
    rho = np.exp(C.ln10 * src.logrho)[None, :]
    T = np.exp(C.ln10 * src.logt)[:, None]
    assert np.allclose(dsets["f_xi"], src.columns["P"] / rho, rtol=1e-12)
    assert np.allclose(dsets["f_tau"], -T * src.columns["S"], rtol=1e-12)
    assert np.allclose(dsets["f"], src.columns["F"], rtol=1e-12)
    # grid metadata
    assert attrs["xi_min"] == pytest.approx(C.ln10 * src.logrho[0])
    assert attrs["dxi"] == pytest.approx(C.ln10 * (src.logrho[1] - src.logrho[0]))
    # Zbar carried through
    assert attrs["has_Zbar"] == 1
    assert np.allclose(dsets["Zbar"], 10.0)


def test_bicubic_order3_zeros_higher_terms(tmp_path):
    iso = Isotope("he4", 2, 4.0026)
    src = _make_source(iso)
    nd = build_nodal_data(src, order=3)
    for name in ("f_xixi", "f_tautau", "f_xixitau", "f_xitautau", "f_xixitautau"):
        assert np.all(nd.datasets[name] == 0.0)
    # f_xitau is still meaningful for bicubic
    assert np.any(nd.datasets["f_xitau"] != 0.0)
    assert nd.attrs["interp_order"] == 3
