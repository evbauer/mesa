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

## Using WDM_AV in the eosDT cascade

The WDM_AV component is integrated as a level of MESA's `eosDT` cascade, between the
CMS and Skye levels (Decision D11). It is **disabled by default**; set the `&eos`
control `use_WDM_AV = .true.` to enable it. When enabled it takes precedence over
Skye/PC/HELM inside its blend region; outside that region, or at any point its tables
do not cover, it declines and the cascade falls through to the lower levels.

### Where tables live

Processed per-isotope HDF5 tables are read from `WDM_AV_table_dir` (default
`mesa_data_dir/eosWDM_AV_data`), optionally under a table-set subdirectory
`WDM_AV_table_set`, as `eosWDM_AV_<isotope>.h5`. For example

```
&eos
   use_WDM_AV = .true.
   WDM_AV_table_set = 'skye'   ! -> data/eosWDM_AV_data/skye/eosWDM_AV_<iso>.h5
/
```

reads e.g. `data/eosWDM_AV_data/skye/eosWDM_AV_c12.h5`. Build processed tables from
source tables with the Python builder, e.g.

```
python3 -m wdm_av_builder build \
   data/eosWDM_AV_data/source_skye/eosWDM_AV_o16.txt \
   data/eosWDM_AV_data/skye/eosWDM_AV_o16.h5 --order 5
```

### Species policy

For each isotope in the net with mass fraction above
`mass_fraction_limit_for_WDM_AV` (default `1d-4`):

- an exact-isotope table (`eosWDM_AV_<name>.h5`) is used if present;
- otherwise, if `WDM_AV_rescale_missing_isotopes = .true.` (default), the nearest-A
  same-Z table is used with the exact density rescaling (Decision D6);
- otherwise the species is uncovered. Species below the mass-fraction limit are
  dropped and the mass fractions renormalised over the covered species. If the total
  uncovered mass fraction exceeds the limit, the component declines the point.

`WDM_AV_mixing_entropy` selects the ideal mixing-entropy model
(`'none'`, `'ions'`, or `'ions_electrons'`, default `'ions_electrons'`), evaluated
with the actual additive-volume partial volumes (Decision D8). Composition
derivatives are not provided in this phase (`d_dxa = 0`, Decision D17).

### Blend controls

The blend region is a simple trapezoidal box in `logRho` and `logT` set by the
`..._for_any_WDM_AV` / `..._for_all_WDM_AV` controls: outside the "any" limits the
component is not used, inside the "all" limits it is fully used, with linear ramps
between (`Do_Blend` applies the quintic smoothstep). The defaults target the low-T
"hole" between FreeEOS and Skye (Decision D20) and are **provisional**:

```
logRho_min_for_any_WDM_AV = -3.0   logRho_min_for_all_WDM_AV = -2.5
logRho_max_for_all_WDM_AV =  4.5   logRho_max_for_any_WDM_AV =  5.0
logT_min_for_any_WDM_AV   =  2.9   logT_min_for_all_WDM_AV   =  3.0
logT_max_for_all_WDM_AV   =  7.0   logT_max_for_any_WDM_AV   =  7.5
```

Solver controls: `WDM_AV_newton_tol` (default `1d-12`), `WDM_AV_max_iter`
(default `50`), and `WDM_AV_debug` (default `.false.`; when true, solver failures
return `ierr` instead of skipping and a Gibbs-space cross-check is printed). See
`eos/defaults/eos.defaults` for the full list with documentation.

