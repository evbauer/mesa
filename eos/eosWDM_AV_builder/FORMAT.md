# WDM_AV table formats (contract between the Python builder and the Fortran reader)

Format version: **1**. Any change to names, shapes, units, or semantics bumps the
version; the Fortran reader refuses files whose `format_version` it does not know.

All quantities are cgs. Per-gram quantities are per gram of the pure isotope.
Natural logarithms are written `ln`; base-10 logarithms `log10`.

## 1. Source tables (input to the builder)

What collaborators deliver. Either HDF5 or ASCII; the builder reads both. One file per
isotope. The ASCII form (see also `claude_tex_notes/wdm_av_table_requests.tex`):

```
# WDM_AV source table, version 1
isotope              ne22          # MESA isotope name (h1, he4, c12, o16, ne20, ne22, mg24, ...)
Z                    10
A                    21.991385     # atomic mass used in the model (amu); integer mass number acceptable
grid                 log10 uniform
nlogRho              321     logRho_min  -8.0    logRho_max   8.0
nlogT                201     logT_min     3.0    logT_max    13.0
includes_radiation   F
energy_zero          "<text>"
entropy_zero         "<text>"
constants            "<text>"
model                "<text>"
tolerance            1e-12
columns              logRho logT P E F [S] [Zbar] [dPdRho_T] [dPdT_Rho] [Cv] [eta] [nheavy]
# rows: logT outer loop, logRho inner loop; NaN for undefined entries
```

Column meanings: `P` dyn/cm^2; `E`, `F` erg/g; `S` erg/g/K; `Zbar` free electrons per
nucleus; `dPdRho_T` (dP/drho)_T; `dPdT_Rho` (dP/dT)_rho; `Cv` erg/g/K; `eta`
electron degeneracy parameter; `nheavy` heavy particles per nucleus. At least one of
`F`, `S` must be present (both preferred). Keys are case-insensitive; `#` starts a
comment. Any additional header key is kept as provenance.

HDF5 source form: the same keys as root attributes, datasets named as the columns with
shape `(nlogT, nlogRho)` in h5py (logT is the slow axis), plus 1-D `logRho`, `logT`.

Synthetic source tables produced by the builder (`synth` command) and the MESA-generated
tables (`eosWDM_AV_builder/mesa_tables/`) use exactly this ASCII form.

## 2. Processed tables (output of the builder, input to MESA)

HDF5, one file per isotope: `data/eosWDM_AV_data/eosWDM_AV_<isotope>.h5`,
e.g. `eosWDM_AV_ne22.h5`. Read in Fortran through `forum` (`hdf5io_t`).

### Root attributes

| name | type | meaning |
|---|---|---|
| `format_version` | int | 1 |
| `isotope` | string | MESA isotope name, e.g. `ne22` |
| `element` | string | e.g. `Ne` |
| `Z` | int | nuclear charge |
| `A` | double | atomic mass (amu) used for number densities: n = X rho/(A m_u) |
| `interp_order` | int | 3 (bicubic: only `f, f_xi, f_tau, f_xitau` are meaningful) or 5 (biquintic) |
| `nxi`, `ntau` | int | grid sizes |
| `xi_min`, `dxi` | double | ln(rho) of the first node and the uniform spacing |
| `tau_min`, `dtau` | double | ln(T) of the first node and the uniform spacing |
| `includes_radiation` | int | must be 0 (the builder strips radiation) |
| `e_offset`, `f_offset` | double | constants added to the source E and F (erg/g), zero if none |
| `has_Zbar`, `has_eta` | int | 0/1 presence flags for optional datasets |
| `fd_scheme` | string | how the estimated derivatives were formed (free text, versioned by builder) |
| `second_derivs_from_source` | int | 1 if `f_xixi, f_tautau, f_xitau` came from supplied dPdRho_T, dPdT_Rho, Cv |
| `source_file`, `source_model`, `builder_version`, `created` | string | provenance |

### Datasets (all `double`)

Written from h5py with shape **`(ntau, nxi)`** (tau is the slow axis). Fortran/forum
then receives arrays dimensioned **`(nxi, ntau)`**, i.e. `f(ixi, itau)`. (This is the
same reversal the AESOPUS opacity reader relies on.)

Required, all erg/g, derivatives with respect to xi = ln rho and tau = ln T:

```
f  f_xi  f_tau  f_xixi  f_tautau  f_xitau  f_xixitau  f_xitautau  f_xixitautau
```

For `interp_order = 3` the last five are still present but may be zero.

Optional (shape `(ntau, nxi)` unless noted):

```
Zbar  Zbar_xi  Zbar_tau  Zbar_xitau      # bicubic nodal data for mean ionization
eta   eta_xi   eta_tau   eta_xitau       # bicubic nodal data for degeneracy parameter
valid                                    # int8, 1 = usable node, 0 = unusable (P<=0, NaN, ...)
```

Optional group `audit/` with diagnostic maps (not read by MESA): `maxwell_residual`,
`euler_maclaurin_xi`, `euler_maclaurin_tau`, `chi_rho_jump`, `cv_jump`, `P_rho_fd`,
`cv_fd`.

### Relations the builder guarantees at every valid node

```
f_xi  = P/rho                 (exact from the source P)
f_tau = -T s = f - e          (from S if present, else from F - E)
```

and, when second derivatives are supplied by the source,

```
f_xixi   = (dP/drho)_T - P/rho = P_rho - f_xi
f_xitau  = T (dP/dT)_rho / rho
f_tautau = f_tau - T c_V
```

(Derivation: f_xixi = d(P/rho)/dxi = rho d(P/rho)/drho = P_rho - P/rho.)

Otherwise these three and the three higher cross terms come from finite differences
described in `fd_scheme`.

## 3. Isotope naming and lookup

Files are found by MESA isotope name: `chem_isos% name(chem_id)` -> `eosWDM_AV_<name>.h5`
in `mesa_data_dir/eosWDM_AV_data/` (overridable by the `WDM_AV_table_dir` control).
Species mapping policy is Decision D6 in `claude_md_notes/04_DECISIONS.md`.
