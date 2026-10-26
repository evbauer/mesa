# make_mesa_tables (WDM_AV WP4)

A small standalone Fortran program that links the built MESA `eos` library and
writes **WDM_AV source-format ASCII tables** (see
`../FORMAT.md` Section 1) for single, pure isotopes, evaluated directly from
MESA's analytic EOS components (HELM and Skye). It also writes **reference
mixture tables** (the full MESA result vector plus first partials) for later
validation of the additive-volume mixture against Skye/HELM run directly on the
mixture.

These tables have exactly known analytic answers and serve as validation data
for the WDM_AV builder (WP1) and runtime (WP2/WP3).

## What it produces

Under `$(MESA_DIR)/data/eosWDM_AV_data/` (git-ignored), with
`which_eos = eosDT_get_component(...)` for the **pure isotope**
(`species = 1, xa = 1`) at every grid node:

| directory | component | contents |
|---|---|---|
| `source_helm/eosWDM_AV_<iso>.txt` | `i_eos_HELM` | electrons/positrons + ideal ions (+ Coulomb per HELM defaults) |
| `source_skye/eosWDM_AV_<iso>.txt` | `i_eos_Skye` | Skye: ions + electrons + Coulomb (liquid/solid) |
| `reference_skye_mixtures/mix_*.txt` | `i_eos_Skye` | full result vector for the reference mixtures |
| `reference_helm_mixtures/mix_*.txt` | `i_eos_HELM` | full result vector for the reference mixtures |

Default isotopes: `h1 he4 c12 o16 ne20 ne22 mg24`.
Reference mixtures (mass fractions):
1. `c12=0.50 o16=0.50`
2. `o16=0.55 ne20=0.35 mg24=0.10`
3. `he4=0.50 c12=0.50`
4. `h1=0.70 he4=0.30`
5. `c12=0.49 o16=0.49 ne22=0.02`

### Source-table columns (per isotope)
```
logRho logT P E F S Zbar dPdRho_T dPdT_Rho Cv eta
```
All cgs, gas-only (radiation stripped), per gram of the pure isotope.
`P = Pgas`, `E`, `F = E - T S`, `S` gas-only; `Zbar = Z` (full ionization);
`dPdRho_T = (dPgas/dRho)_T`, `dPdT_Rho = (dPgas/dT)_rho`; `Cv` gas-only; `eta`
is the electron degeneracy parameter. Rows: `logT` outer loop, `logRho` inner;
failed nodes are written as `NaN`.

### Reference-mixture columns
```
logRho logT lnPgas lnE lnS chiRho chiT Cv Cp grad_ad gamma1 gamma3 eta mu
lnfree_e dlnPgas_dlnRho dlnPgas_dlnT dlnE_dlnRho dlnE_dlnT dlnS_dlnRho
dlnS_dlnT dCv_dlnRho dCv_dlnT dgrad_ad_dlnRho dgrad_ad_dlnT
```
These are the **raw** MESA component outputs (`res`, `d_dlnd`, `d_dlnT`); no
radiation stripping is applied (so the AV mixture can be compared against Skye
on the same footing).

## Build and run

```bash
export MESA_DIR=/path/to/this/mesa     # MESA must already be built under $MESA_DIR/build
cd eos/eosWDM_AV_builder/mesa_tables
make                                   # builds build/make-mesa-tables/bin/make-mesa-tables
make tables                            # build + run (reads inlist_mesa_tables from CWD)
```
or run the binary directly from this directory (it reads `inlist_mesa_tables`
from the current directory and uses `$MESA_DIR` for data):
```bash
../../../build/make-mesa-tables/bin/make-mesa-tables
```
Runtime is a few minutes for the default grid. The program prints per-table
diagnostics (failed-node counts and their bounding box, S<=0 count,
P-monotonicity, the ideal-gas check and the finite-difference `dP/drho` check).
MESA prints one `failed in Test_one_eosDT_component` line per failed node to
stdout; redirect to a log and filter it out when reading diagnostics.

## Configuration (`inlist_mesa_tables`)

Two namelists:
- `&eos` — read by `alloc_eos_handle_using_inlist`; sets
  `include_radiation = .false.`, `use_Skye`, `Skye_use_ion_offsets`, etc.
- `&table_controls` — grid (`logRho_min/max`, `logT_min/max`), resolution
  (`pts_per_decade` for isotopes = 20, `pts_per_decade_mix` for mixtures = 10),
  the isotope list (`n_iso`, `iso_names(:)`), and `output_root` (blank =
  `$(MESA_DIR)/data/eosWDM_AV_data`).

## Radiation handling (important)

The output source tables are gas-only (`includes_radiation = F`). The two
components honour `include_radiation` differently in the current source:

- **HELM** respects `include_radiation` (`helm_sum_totals.dek`): with it
  `.false.`, P, E, S, chiRho, chiT and Cv are all gas-only. No subtraction
  needed.
- **Skye ignores `include_radiation` entirely.** `skye_eos` / `pack_for_export`
  (`skye_thermodynamics.f90`) always write `res(i_lnPgas)` as the gas pressure
  (good) but always add radiation to E, S, Cv, and compute chiRho/chiT from the
  **total** pressure. For Skye the program therefore subtracts radiation
  analytically for the source tables:
  ```
  e_rad  = crad*T^4/rho
  s_rad  = 4*crad*T^3/(3*rho)
  cv_rad = 4*crad*T^3/rho
  ```

Because Skye's `chiRho`/`chiT` refer to the *total* pressure, the gas-only
pressure derivatives are computed from the gas-pressure log-derivatives for
**both** components:
```
dPdRho_T = Pgas * d_dlnd(i_lnPgas) / rho
dPdT_Rho = Pgas * d_dlnT(i_lnPgas) / T
```
For HELM (radiation off) this equals `chiRho*P/rho`; for Skye it is the correct
gas-only value (whereas `chiRho*P/rho` would wrongly use the total-pressure
susceptibility). The program's finite-difference check confirms this agrees
with a central difference of `Pgas` across neighbouring density nodes
(relative error ~1e-9).

## Energy / entropy zero conventions (recorded in each table header)

- **HELM** — electron/positron kinetic energy (no electron rest mass) + ideal
  ion thermal energy + Coulomb (HELM defaults); no ionization offset; radiation
  excluded.
- **Skye** — ideal ion + electron/positron + Coulomb free energy;
  `Skye_use_ion_offsets = .true.` (default) adds the total electron binding
  (full-ionization) energy per gram so the fully-ionized ground state is the
  zero of energy (consistency with FreeEOS); radiation stripped analytically.

## Expected failures

Forcing a single component over the whole plane means each component fails
outside its domain of validity, which is expected and written as `NaN`:
- HELM fails in the cold regime (`helmeos2` errors when gas energy/entropy/
  pressure fall below 1e-20 — the HELM electron energy without rest mass goes
  negative at low T), i.e. exactly where MESA normally blends to Skye/PC/FreeEOS.
- Skye fails where it returns `s<0`, `e<0` or `pgas<0` (extreme corners).
