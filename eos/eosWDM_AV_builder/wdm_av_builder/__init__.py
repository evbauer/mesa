"""wdm_av_builder: Python table builder for the MESA WDM_AV (warm dense
matter, additive-volume) equation of state.

The package reads single-isotope source tables (P, E, F, [S], ... on a uniform
log10 rho-T grid), estimates the Hermite nodal data in log coordinates
(xi = ln rho, tau = ln T), audits consistency, and writes processed HDF5
tables for the Fortran reader.  See ``FORMAT.md`` for the file contracts and
``claude_md_notes/`` for the design.

Modules
-------
constants   physical constants copied from MESA's const_def.f90
synthetic   analytic ideal-gas test EOS with exact log-coordinate derivatives
io          read/write source (ASCII/HDF5) and processed (HDF5) tables
derivs      nodal-data builder: exact first derivatives + FD higher derivatives
audit       table-consistency diagnostics (Euler-Maclaurin, Maxwell, jumps)
plots       matplotlib maps of audit quantities
cli         ``python -m wdm_av_builder synth|build|audit|check``
"""

from __future__ import annotations

__version__ = "0.1.0"
