# Bounded standalone performance control (experiment 017)

This **sidecar only** benchmarks the existing metadata-aware `freeeos_table`
adapter in one persistent Fortran process per replicate. It neither installs
that adapter in MESA nor changes a source table, cache, kernel, module, runtime
route, default, or another worker's analysis/report. All generated objects,
executables, logs, query copies, runtime cache links, and results live in
`freeeos_biquintic_progress/experiments/017_performance_control/`. Code changes
are confined to this directory. `mesa_probe/` is a read-only design reference,
not a build/run dependency. No MESA or table module recompilation is attempted.

## Reproduction and hard bounds

From the workspace root with the existing Python/numpy installation:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s mesa/eos/eosFreeEOS_builder/biquintic/performance_control -p 'test_*.py' -v
PYTHONDONTWRITEBYTECODE=1 python3 mesa/eos/eosFreeEOS_builder/biquintic/performance_control/run_benchmark.py --smoke
PYTHONDONTWRITEBYTECODE=1 python3 mesa/eos/eosFreeEOS_builder/biquintic/performance_control/run_benchmark.py
PYTHONDONTWRITEBYTECODE=1 python3 mesa/eos/eosFreeEOS_builder/biquintic/performance_control/capture_preflight.py
PYTHONDONTWRITEBYTECODE=1 python3 mesa/eos/eosFreeEOS_builder/biquintic/performance_control/finalize_results.py
```

The wrapper builds **only the harness** with the local gfortran-15 SDK and
installed libraries. It links the already-built `table_fortran/build/`
`freeeos_table.o`, `freeeos_biquintic.o`, `helm_polynomials.o` and module files.
Missing objects/libraries are fatal; do not rebuild dependencies to satisfy
this control. Harness flags are `-O2`, no fast math, `-ffp-contract=off`, and
OpenMP enabled but one thread. Existing objects' original flags are not
retrospectively inferred from the current adjacent Makefile. Their hashes,
the current source/Makefile hashes, compiler version, and actual build command
are recorded; installed EOS library flags may differ. A changed harness is
rebuilt by make without writing to the adjacent build directory.

The default wrapper wall bound is **240 seconds including build and every
process**, capped at 280 seconds. It also stops launching/running tasks at the
user's **2026-10-28 20:18 PDT** fence. Each subprocess receives only the remaining
whole-run time; a timeout kills it and records partial results. No debugger,
dependency download, permissions request, commit, or push is used. The explicit
date fence intentionally prevents later reruns without revisiting the deadline.

Defaults: 3 fresh processes/backend, 7 warmed timing batches/mode/process,
0.35-second batch target, all 13,577 archived 014 faces. Pass counts are
calibrated per mode/process and capped at 10,000; at least one complete corpus
pass is necessary even if it exceeds the target. A fixed seed 20261028 shuffles
the corpus once; every mode uses the identical order repeatedly. CPU affinity
is the first allowed CPU unless `--cpu N` selects another allowed CPU. Both
executables are pinned identically; neither priority nor other workers changes.

For reduced controls use `--limit 1024 --process-repeats 1 --repeats 3
--target-seconds .02 --time-limit 120`. `--no-mesa` retains the required four
standalone paths. `--smoke` uses the **existing** `ideal_rad.cache` fixture,
five analytic queries, and short timing controls, writing into `017/smoke/`.
It generates no source table or cache, and checks real/AD gas+radiation state
values and Jacobians against analytic formulas/central differences. Unit tests
cover coordinate conversion, deterministic corpus order, parsing, timing units,
and the analytic gas limit. These are numerical smoke checks, not a new source
accuracy or stability assessment.

## Paths and measurements

- `real_eval`: physical natural-log `(ln rho, ln T)` input, metadata-aware cell
  lookup/mask, real kernel, coordinate chain rule, declared radiation restoration;
  returns ten total-potential partials through third order.
- `ad_eval`: same table/lookup/radiation with `auto_diff_real_2var_order3`,
  identity seeds rho=`d1val2`, T=`d1val1`, returning the same ten physical partials.
- `real_state`: `real_eval` plus existing `derive_state`, yielding eleven physical
  states and their 22 first derivatives (43 scalars including potential partials).
- `ad_state`: extracts the ten AD partials into `potential_partials`, then calls
  the **same** `derive_state` with the loaded metadata's `crad`. This is not a
  separate AD-over-all-thermodynamic-formulas implementation or MESA Jacobian.
- Optional `default_eosDT`: new Fortran harness initializes installed MESA,
  calls `eosDT_get` directly in process, and consumes all basic fields, rho/T
  derivatives, and abundance derivatives. Composition, blend fractions, and
  output count are recorded. `inlist_benchmark` disables only WDM_AV, as in 014;
  this is otherwise the typical **installed** default dispatch, not a claim about
  every MESA release's defaults.

Query preparation is untimed: 014's saved `u_q/u_t`, `tau_q/tau_t` are base-10
`logQ/logT`, converted via `ln rho=(logQ+2logT-12)*ln10`, `ln T=logT*ln10`.
For MESA, matching linear and base-10 coordinates are also prepared once.
The wrapper snapshots the archive bytes before loading and records its hash
and selected original indices; it never writes back to 014. The optional
source-archive-change flag permits a concurrent lead update without attributing
new bytes to the saved snapshot. All four standalone paths must succeed on the
**same full selected corpus**; no rejection-based selection or physical clipping
is silently used. The minimum implementation needs only the 014 archive;
003 saved queries are not regenerated or necessary here.

Startup/cache load, query read, full-domain correctness preflight, explicit
complete warm passes, pass-count calibration, and repeated warmed loops are
separate report records. The preflight exposes any MESA lazy component loading;
after-load and warmed RSS distinguish those phases. Timing uses int64
`system_clock` wall time and `cpu_time`, not Python pipe/startup throughput.
No formatted query output occurs in timed loops. Each successful call packs
outputs and accumulates `sum(value/field_scale)`, with one scale/field obtained
untimed from the maximum absolute output or 1. Real and AD share the exact same
scales. Status branches, output packing, checksum arithmetic, and callback
dispatch are included and disclosed, not subtracted as fictitious EOS work.

`results.json` retains every batch, wall/CPU times, checksums/pass, per-process
startup/memory/validation statistics, raw field sums, and all 43 per-field
maximum real/AD differences scaled by `max(1,abs(real),abs(AD))`. Rotating mode
order reduces consistent first/last-mode bias. Median, p10/p90, min/max, relative
sample deviation, and an approximate process-level bootstrap median interval
describe scatter. Three fresh processes do **not** justify a precise statistical
confidence claim: batches correlate and shared-node load, frequency, thermal,
compiler, and cache effects remain systematic uncertainties.

## Memory and cache safety

The standalone text cache is `008_standalone_table/solar_split.cache`, loaded
once per process with expected X=.70/Z=.02 and its declared source hash.
The full cache bytes and linked inputs are hashed before/after. The module
retains provenance/radiation metadata; header hashes alone do not authenticate
unread source-table bytes. Numeric payload is `nq*nt*(9*sizeof(real(dp)) +
sizeof(logical)) + (nq+nt)*sizeof(real(dp))`: invalid nodes remain allocated too.
This excludes descriptors, metadata, allocator overhead, and other process
memory. `/proc/self/status` supplies RSS and high-water RSS in KiB; neither
resident memory nor peak RSS should be called an exact allocated cache size.

MESA gets private 017 runtime/temp cache roots, with symlinks to **existing**
installed EOS bins and temporary-cache safety required. Missing installed bins
are not generated by this wrapper. Installed size/mtime inventories and private
link integrity are audited even after failure/timeout; any replaced link or
temporary file fails the control. Installed caches are not rehashed in full:
the inventory check is not proof against a same-size/same-mtime external edit.
Total installed disk cache bytes include all compositions/components, not just
what a query loads; resident MESA memory includes more than EOS tables. Neither
OS page cache nor the installed caches are flushed, so fresh-process startup
is **not cold disk-cache performance**.

## Fair interpretation and artifacts

Real versus AD table modes have matched coordinates, masks, radiation,
derivative order, outputs, and state derivation. Their throughput ratios are
fair **for these prebuilt standalone paths and this corpus**. MESA has extra
composition interpolation, routing/blending, outputs, and different installed
optimization flags; standalone paths are fixed X/Z without dispatch. The
face-query corpus is not a production stellar-track frequency distribution.
The full face corpus exposes raw **own-column** relative third-partial/Jacobian
tails. A small derivative near cancellation can have a large relative difference
without a large difference on its parent response's scale. These raw diagnostics
are **not** project response-parent-scaled errors, physical/source errors, or a
blanket field-scaled failure gate. `capture_preflight.py` builds the harness in
017's separate `preflight/build/` and evaluates each of the 13,577 recorded
queries once on each state path. It uses `--preflight-only`, writes all 43 paired
columns, runs no warm/calibration/timed loops, preserves the original timed
binary, checks the same recorded cache SHA before/after, and has a 120-second
whole-run bound plus the 20:18 fence. It generates no table/cache. Once capture
finishes, all subsequent finalization is recorded-only, and 008 may be reexported.

`finalize_results.py` uses captured raw arrays to compute each state derivative's
denominator `max(1,abs(real_state),abs(AD_state),abs(real_partial),abs(AD_partial))`.
The zero-based offsets are potential=0 (10 fields), state values=10 (11),
rho derivatives=21 (11), T derivatives=32 (11). State values self-scale; the
first ten potential columns retain their raw own-column diagnostic, clearly
separated from the project state-response comparison. It records all 43 scales,
differences, maxima, and worst actual values/parent responses. The former
own-column `Cp_x` maximum .0767 becomes about **2.475e-9** on its `Cp` parent
scale; the maximum parent-scaled state Jacobian over the whole corpus is about
**4.898e-7**. Neither normalization constitutes a source-accuracy test. Exact
numerical identity and integrated production speedup remain unclaimed.

**Lead-review resolution:** the former interpretation of raw `Cp_x` .0767 as
a project response-scaled failure is withdrawn. At that same query the correct
parent-scaled difference is **2.474902236e-9**; the full-corpus maximum is
**4.898390989e-7**. Raw own-column tails remain recorded as cancellation-sensitive
diagnostics, not physical-error estimates. This worker's benchmark/capture work
is complete; no further benchmark or solar-cache read is pending, and 008 may
be reexported without changing the frozen measured-cache manifest.

Report MESA throughput as context, **never an exact production speedup or an
authorization to replace eosDT_get**. To claim such a speedup would require
matched production outputs, dispatch, composition work, build options,
representative workload, and end-to-end measurement in an integrated backend.

`cache_manifest.json` separates the exact benchmarked cache SHA and its declared
preprocessing SHA from the **current** `freeeos_potential.py` core source SHA at
finalization. The strict-header core change reported by the user postdates the
old export's provenance; unchanged numerical jets are user-reported, not a new
export comparison here. A subsequent 008 reexport must not replace the measured
cache identity in this manifest. Finalization reads only recorded 017 artifacts
and the current core source; it never rereads the solar cache or blocks reexport.

Artifacts: `RESULTS.md`, `results.json`, `cache_manifest.json`,
`response_parent_scaling.json`, `preflight/capture_manifest.json`,
`preflight/state_trace.stdout`, `preflight/response_parent_traces.npz`, `queries.txt`,
`selected_saved_indices.npy`, `build.stdout/stderr`, `table_NN.stdout/stderr`,
optional `mesa_NN.stdout/stderr`, `build/`, private `runtime_cache/` and
`runtime_temp_cache/`, and `smoke/`. Logs retain raw per-query validation outside
timed loops. Main runs replace only 017's own results; smoke uses its subdirectory.
No lead report/PDF or active-face analysis is edited.
