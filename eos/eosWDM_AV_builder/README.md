# eosWDM_AV_builder

Tools for the WDM_AV (warm dense matter, additive volume) EOS component of MESA:
single-isotope free-energy tables combined on the fly into mixtures.

- `FORMAT.md` - the file-format contract (source tables in, processed HDF5 out).
- `wdm_av_builder/` - Python package: read source tables, audit them, estimate the
  nodal derivatives, write processed HDF5 tables, make synthetic test tables, plots.
- `mesa_tables/` - Fortran program that writes source-format tables for each isotope
  from MESA's own HELM (electrons + ideal ions) and Skye, used for validation.
- `test_fortran/` - standalone Fortran test programs for the eos_wdm_av_* modules
  (built like `eos/plotter`, linking the already-built MESA libraries).

Design notes and status live outside the repository in the project's
`claude_md_notes/` and the scientific writeup in `claude_tex_notes/`.

Build environment: `export MESA_DIR=<this checkout>`; MESA must be built (`./install`
or `make eos` from the top level).
