# WDM_AV table validation

Independent validation of the WDM_AV Hermite EOS tables: a from-scratch numpy
evaluator (checked against the nodal data and the cell-corner identity, i.e.
against the Fortran reader's conventions), and a held-out-node interpolation
accuracy study on the real Skye- and HELM-derived source tables.

Everything here is pure Python and reads/writes HDF5 per `../FORMAT.md`. It does
not touch the Fortran or the build system.

## Files

- `hermite_eval.py` — independent, vectorised numpy implementation of the
  bicubic (order 3) and biquintic (order 5) Hermite tensor-product evaluation of
  a processed table. Written from the mathematics (`eos_mixing_notes.tex` Sec. 8
  and the standard quintic/cubic Hermite bases), **not** transcribed from the
  Fortran. Returns `f, f_xi, f_tau, f_xixi, f_xitau, f_tautau` at arbitrary
  arrays of `(xi, tau) = (ln rho, ln T)`. Key entry points:
  - `load_table(path)` — read a processed HDF5 table into a `HermiteTable`.
  - `table_from_arrays(attrs, datasets)` — wrap in-memory builder output.
  - `evaluate(table, xi, tau)` — evaluate the interpolant and its derivatives.
  - `find_cell(table, xi, tau)` — cell lookup (matches the Fortran convention).
  - `cell_corner_twist_integral(table, ixi, itau)` — the quadrature check used by
    the cell-corner identity test.
  For an order-3 table the second derivatives are the analytic second
  derivatives of the *bicubic* surface (built from `f, f_xi, f_tau, f_xitau`
  only), exactly as the Fortran `hc(..., cd2x)` does.

- `heldout.py` + CLI — the held-out-node accuracy study (see below).

- `tests/` — pytest tests for the evaluator:
  - nodal reproduction (every node reproduces the stored nodal data to roundoff:
    `f, f_xi, f_tau` for order 3; plus `f_xixi, f_tautau, f_xitau` for order 5);
  - the cell-corner identity `∫∫ f_xitau = corner differences of f`;
  - finite-difference consistency of the analytic derivatives with the f surface;
  - accuracy against the analytic synthetic free energy.

- `results/` — machine-readable held-out summaries:
  - `heldout_summary.json` — full nested summary (both regions, all six
    quantities, worst-error locations);
  - `heldout_summary.csv` — one row per (table, variant, quantity, region);
  - `heldout_summary.md` — median/p99/max table for the five primary
    quantities, one block each for the skye_valid, skye_valid_no_pair, and
    target regions.

Error maps are written to `claude_tex_notes/figures/heldout/` as
`<table>_<variant>.png` (6 panels: P, s, c_V, chi_rho, chi_T, e).

## Held-out method

For a source table on the fine grid (20 points/decade), `heldout.py`:

1. coarsens to every other node in both directions (10 points/decade);
2. builds a processed table from the coarse grid, in three variants:
   - `a_o5` biquintic with second derivatives from the source
     `dPdRho_T, dPdT_Rho, Cv` columns;
   - `b_o5` biquintic with second derivatives by finite differences
     (`build --ignore-source-derivs`);
   - `a_o3` bicubic;
3. evaluates each coarse table with `hermite_eval` at the omitted fine nodes;
4. compares the derived quantities against the source values there:

   | quantity | from the interpolant | source truth |
   |---|---|---|
   | P       | `rho f_xi`              | `P` |
   | s       | `-f_tau/T`              | `S` |
   | c_V     | `(f_tau - f_tautau)/T`  | `Cv` |
   | chi_rho | `1 + f_xixi/f_xi`       | `rho dPdRho_T / P` |
   | chi_T   | `f_xitau/f_xi`          | `T dPdT_Rho / P` |
   | e       | `f - f_tau`             | `E` |

A fine node is counted only if it is omitted, its coarse cell is fully valid
(all four coarse corners valid), and the source truth is finite and non-zero.
Statistics (median, p90, p99, max) are reported over four regions:

- `skye_valid` — `logT >= 6.5 OR logRho >= 3.7` (Decision D20);
- `skye_valid_no_pair` — the above, excluding the high-T / low-ρ pair-radiation
  corner (`NOT (logT > 8.0 AND logRho < 0)`), where P is nearly flat in ρ and
  chi_rho cancels;
- `target` — the WDM_AV coverage hole, `logRho in [-2, 4]` and `logT in [3, 7]`;
- `whole_valid` — the whole valid plane.

## Rerunning

From this directory (`eos/eosWDM_AV_builder/validation`):

```bash
# tests
python3 -m pytest tests/
python3 -m pytest ../tests/            # builder tests (incl. --ignore-source-derivs)

# held-out study for all seven Skye isotopes + c12/o16 from HELM
python3 heldout.py \
  ../../../data/eosWDM_AV_data/source_skye/eosWDM_AV_{h1,he4,c12,o16,ne20,ne22,mg24}.txt \
  ../../../data/eosWDM_AV_data/source_helm/eosWDM_AV_{c12,o16}.txt
# -> results/heldout_summary.{json,csv,md} and
#    claude_tex_notes/figures/heldout/<table>_<variant>.png
```

Options: `--variants a_o5 b_o5 a_o3`, `--figdir`, `--resultsdir`, `--label`
(single source). Paths above are relative; absolute paths also work.

## Full-resolution tables (integration-test inputs)

The non-coarsened processed tables (order 5, source second derivatives) are
built straight from the builder CLI:

```bash
cd ..    # eos/eosWDM_AV_builder
for iso in h1 he4 c12 o16 ne20 ne22 mg24; do
  python3 -m wdm_av_builder build \
    ../../data/eosWDM_AV_data/source_skye/eosWDM_AV_${iso}.txt \
    ../../data/eosWDM_AV_data/skye/eosWDM_AV_${iso}.h5 --order 5
  python3 -m wdm_av_builder build \
    ../../data/eosWDM_AV_data/source_helm/eosWDM_AV_${iso}.txt \
    ../../data/eosWDM_AV_data/helm/eosWDM_AV_${iso}.h5 --order 5
done
```
