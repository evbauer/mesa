# New shifted native edge-cell bounds — experiment 035

Exclusive sidecar: code only here, all generated artifacts only in
`freeeos_biquintic_progress/experiments/035_edge_bounds/`. It imports read-only
`positivity_bounds/{bounds.py,run_bounds.py,bound_split.py}` and
`edge_control/{edge_control.py,clean_floor.py}`. No existing math/kernel/default,
source table, runtime cache, old experiment, or other-worker output is changed.
No fit, compile, MESA call, debugger, approval, commit, or push is involved.

## Reproduce / resource policy

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s mesa/eos/eosFreeEOS_builder/biquintic/edge_bounds -p 'test_*.py' -v
timeout 190s env PYTHONDONTWRITEBYTECODE=1 python3 mesa/eos/eosFreeEOS_builder/biquintic/edge_bounds/run_control.py
PYTHONDONTWRITEBYTECODE=1 python3 mesa/eos/eosFreeEOS_builder/biquintic/edge_bounds/verify_artifacts.py
```

Two candidates only: **shifted quartic native jets, radiation split, source floor
3.2 default / 3.1 experimental**. Each selects only newly complete cells versus
the **same fixed centered radius2-mean floor3.2 split baseline**, not versus a
centered 3.1 baseline. Require all four candidate corners native valid and
complete. The 20,000-new-cell cap per candidate **refuses** excess populations;
no unannounced subsample/truncation. Depth defaults/maxes at 3, visited nodes per
field/cell at 85, batch size 512 (cap1024), whole audit budget180s with a default
**daily local 20:00** finish fence. `--finish-by` accepts an optional timezone-aware
ISO timestamp; on 2026-10-28 the effective fence cannot exceed 20:00 PDT even with
a later override. Later days are not disabled by the old date. Both requested
and effective invocation fences are recorded. The enclosing190s timeout also bounds imports,
source preparation and final file writes. A runtime cap leaves explicit
unresolved/unprocessed counts; incomplete calibration suppresses all positive
certificates for that candidate. `--overwrite` replaces only 035-owned output.

Read the pinned original source SHA
`e6ca7ce27bee8197bb257a8c4886a06366e6951e58958fbed3735c877ee68502`.
The existing `build_control`/`build_jets` builds both candidate masks and native
jets. Its recursive support guards are preserved; additionally archive/check
the selected four-corner original donor bounds against the declared source
floor and absolute **3.1 minimum**. Native flags, printed-coordinate validity,
and all other exclusions remain intact. No below3.1 repair or extrapolation.

## Sign polynomial and restored-core safeguards

Reuse longdouble absolute-Hermite assembly, physical chain-rule derivatives,
closed-cell stretching to actual stored endpoints, and Bernstein convex-hull
bounds. `T*Cv = f_tau - 2*f_u - f_tautau + 4*f_utau - 4*f_uu` and
`P_rho = f_u + f_uu`. The analytic radiation Cv numerator is
`4*CRAD*10**(2*logT-logQ+12)`, monotonically bounded on each leaf. The analytic
radiation Prho contribution is zero, **but the restored-core subtraction
cancellation guard is not zero**. Both polynomial/core margins and existing
analytic/restoration guards are used unchanged. No weak-margin fallback.

Root-positive cells need no subdivision. Ambiguous roots use the existing
bounded `adaptive_split`, recomputing analytic radiation bounds per leaf.
Only an actual finite negative **restored-core** Cv or Prho from a forced-cell
query is a negative witness; a negative coefficient or lower bound is not.
Four restored-core queries/cell use the exact024 phases and saved cell ownership;
compare polynomial-plus-analytic values to the full effective margin, and save
all values/coordinates/margins. Negative calibration samples also supply actual
witnesses, not coefficient-only failures. A certificate that contradicts a
restored-core query fails the run. Any nonfinite core query or out-of-margin
comparison suppresses positive certificates. Fixed-input numerical scope only:
not directed-rounding formal interval arithmetic, source accuracy, whole-star
stability, or runtime authorization.

Classify `full_one_branch_own_use` **exactly with the existing024 `domain_masks`
helper**. This is a closed-cell geometric routing predicate, not measured
activity or an inferred actual MESA fraction; the source-floor3.1 cells below3.2
cannot satisfy its unchanged temperature condition. Material corner masks are
archived as context, not silently substituted for the selected population.

## Artifacts / tests

Top-level `results.json` and concise `README.md` table; per-floor `results.json`
and `cell_bounds_and_core_queries.npz`. They retain source/API/sidecar hashes,
read-only before/after inventories, exact original cell indices, candidate and
baseline masks, corner-valid flags, donor bounds/floor metadata, statuses,
adaptive resources, explicit unresolved cells, negative-witness coordinates,
and four restored-core comparisons per cell. Tests reject excessive resources,
baseline overlap, invalid source corners, below-floor donors, nonfinite jets,
and cancellation-guard omissions, and verify the closed-cell polynomial.
No jets are installed or exported as a runtime cache/source table.

The final recorded-only verifier checks exact new-cell masks, four native
corners, donor floor, unchanged024 geometric own-use classification, all four
core queries/margins, resource caps, and finite-negative witness counts against
the saved NPZ. It writes `validation_results.json`, fingerprints the frozen
sidecar, and performs no source/jet rebuild, MESA call, or cache access.

## Observed bounded gate

Both candidates finish in **3.63 seconds**. Floor3.2 has **9,464** newly complete
native cells; floor3.1 has **11,894**, relative to the same centered3.2 baseline.
Both contain **5,069** full-one-branch geometric own-use new cells, all positive
for Cv numerator and Prho. Outside that geometric subset, each candidate has
**72 Cv and 16 Prho negative-witness cells**; unresolved counts are **zero**.
Each cell receives four restored-core queries; all comparisons pass the full
effective guard margin. The extra experimental low-floor coverage is not an
actual-MESA activity measurement or a runtime-safe fallback authorization.

The recorded-only finalizer additionally computes exact original-coordinate
set overlap and floor-only populations, with closed-cell logT extents and each
field's positive/witness/unresolved counts. `floor_comparison.npz` retains all
available per-cell field arrays paired on overlap (including query values,
guards, statuses and resources), and floor-only records. An unsupported other
floor has status -1, never unresolved or invented numeric predictions. Lowering
the source floor can change overlap donor support and jets; no equality is
assumed. `floor_comparison.json` records input/archive hashes, observed paired
changes and geometric024 own-use counts, not actual MESA fractions. A test
checks reordered coordinate identity, disjoint cells and duplicate rejection.
