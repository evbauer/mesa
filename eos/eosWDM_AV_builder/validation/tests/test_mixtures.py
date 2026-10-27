"""Tests for the reference-correction helpers in ``mixtures.py``.

Covers the three reference corrections applied before the WDM_AV mixture
comparison (brief WP4b task 4):

* ``delta_s_ideal``  -- Skye double-counts it (subtracted), HELM omits it (added);
* ``ion_offset_dE``  -- Skye's ion offset is non-additive (added to the Skye E).

Plus a synthetic "the correction recovers additivity" check and a reader check.
"""

from __future__ import annotations

import mixtures as mx
import numpy as np


def test_delta_s_ideal_pure_is_zero():
    # a single species has no mixing entropy
    assert mx.delta_s_ideal({"c12": 1.0}) == 0.0


def test_ion_offset_dE_pure_is_zero():
    # a pure species: number and mass weighting coincide -> no non-additivity
    assert mx.ion_offset_dE({"c12": 1.0}) == 0.0
    assert mx.ion_offset_dE({"he4": 1.0}) == 0.0


def test_ion_offset_dE_co_matches_manager_value():
    # C/O 50/50: manager-established 7.0e13 erg/g (report gives 6.98e13)
    dE = mx.ion_offset_dE({"c12": 0.5, "o16": 0.5})
    assert np.isclose(dE, 6.98e13, rtol=5e-3)


def test_ion_offset_dE_hec_matches_manager_value():
    # He/C 50/50: manager-established 2.29e14 erg/g
    dE = mx.ion_offset_dE({"he4": 0.5, "c12": 0.5})
    assert np.isclose(dE, 2.29e14, rtol=5e-3)


def test_delta_s_ideal_co_matches_report_value():
    ds = mx.delta_s_ideal({"c12": 0.5, "o16": 0.5})
    assert np.isclose(ds, 4.141e6, rtol=5e-3)


def test_ion_offset_correction_recovers_additivity():
    """The Skye mixture ion offset uses number fractions; the mass-weighted sum
    of the pure-species offsets (what WDM_AV carries) is recovered by adding
    ``ion_offset_dE``.
    """
    comp = {"c12": 0.5, "o16": 0.5}
    n = {k: comp[k] / mx.atomic_weight(k) for k in comp}
    ntot = sum(n.values())
    y = {k: n[k] / ntot for k in comp}
    scale = mx.EV2ERG / mx.C.amu
    mixture_offset = sum(y[k] * mx.ION_EV[mx.charge(k)] for k in comp) * scale
    massweighted_offset = sum(comp[k] * mx.ION_EV[mx.charge(k)] for k in comp) * scale
    assert np.isclose(
        mixture_offset + mx.ion_offset_dE(comp), massweighted_offset, rtol=1e-12
    )


def test_helm_entropy_correction_recovers_additivity():
    """HELM carries no ideal mixing entropy, so adding ``delta_s_ideal`` to a
    purely-additive synthetic HELM reference S reproduces the WDM_AV convention
    (additive + one mixing term).
    """
    comp = {"c12": 0.5, "o16": 0.5}
    s_additive = 1.234e8  # arbitrary mass-weighted pure-species entropy
    ds = mx.delta_s_ideal(comp)
    s_helm_ref = s_additive  # HELM: no mixing entropy
    assert np.isclose(s_helm_ref + ds, s_additive + ds)


def test_reader_roundtrip(tmp_path):
    """read_header / read_table parse a synthetic reference-format file."""
    p = tmp_path / "mix_test.txt"
    p.write_text(
        'composition          "c12=0.50 o16=0.50"\n'
        "nlogRho              2     logRho_min    -8.0000    logRho_max    10.0000\n"
        "nlogT                1     logT_min       3.0000    logT_max       9.5000\n"
        "# a comment line\n"
        + " ".join(["0.0"] * len(mx._REF_NAMES))
        + "\n"
        + " ".join(["1.0"] * len(mx._REF_NAMES))
        + "\n"
    )
    hdr = mx.read_header(p)
    assert mx.parse_composition(hdr["composition"]) == {"c12": 0.5, "o16": 0.5}
    arr = mx.read_table(p, len(mx._REF_NAMES))
    assert arr.shape == (2, len(mx._REF_NAMES))
