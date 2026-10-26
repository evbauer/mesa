"""Physical constants (cgs), copied verbatim from MESA so that the Python
builder and the Fortran reader agree to the last digit.

Source: ``mesa/const/public/const_def.f90``.  Line numbers refer to that file
as of the WDM_AV work (2026-10-26).  Only the constants the builder needs are
reproduced here; do not edit the values without re-copying from the source.

All quantities cgs.  Natural logs are used internally (xi = ln rho,
tau = ln T); ``ln10`` converts to/from the base-10 logs at the file boundary.
"""

from __future__ import annotations

import math

# --- straight copies from const_def.f90 ------------------------------------
# line 44:  pi = 3.1415926535897932384626433832795028841971693993751d0
pi = 3.1415926535897932384626433832795028841971693993751
# line 52:  ln10 = 2.3025850929940455_dp  ! = log(10d0)
ln10 = 2.3025850929940455
# line 71:  avo = 6.02214076d23  ! Avogadro constant (mole^-1)
avo = 6.02214076e23
# line 72:  amu = 1d0/avo  ! atomic mass unit (g)
amu = 1.0 / avo
# line 73:  clight = 2.99792458d10  ! speed of light in vacuum (cm s^-1)
clight = 2.99792458e10
# line 75:  kerg = 1.380649d-16  ! Boltzmann constant (erg K^-1)
kerg = 1.380649e-16
# line 77:  planck_h = 6.62607015d-27  ! Planck constant (erg s)
planck_h = 6.62607015e-27
# line 78:  hbar = planck_h/(2*pi)
hbar = planck_h / (2.0 * pi)
# line 101: me = 9.1093837015d-28  ! electron mass (g)
me = 9.1093837015e-28

# --- derived, matching const_def.f90 exactly -------------------------------
# line 86: boltz_sigma = (pi^4 kerg^4)/(60 hbar^3 c^2)  (Stefan-Boltzmann)
boltz_sigma = (pi * pi * kerg * kerg * kerg * kerg) / (
    60.0 * hbar * hbar * hbar * clight * clight
)
# line 87: crad = boltz_sigma*4/clight  ! radiation density constant "a"
#          (erg cm^-3 K^-4);  Prad = crad * T^4 / 3
crad = boltz_sigma * 4.0 / clight

# convenience
iln10 = 1.0 / ln10  # const_def.f90 line 53 gives 0.43429448190325187

__all__ = [
    "amu",
    "avo",
    "clight",
    "crad",
    "hbar",
    "iln10",
    "kerg",
    "ln10",
    "me",
    "pi",
    "planck_h",
]


def _self_check() -> None:
    """Cheap sanity check that the derived constants match the source value."""
    # const_def.f90 effectively yields crad ~ 7.5657e-15 erg cm^-3 K^-4
    assert math.isclose(crad, 7.565733250945076e-15, rel_tol=1e-10), crad
    assert math.isclose(ln10, math.log(10.0), rel_tol=1e-15)


_self_check()
