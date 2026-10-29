# Native shifted-stencil coverage control

Optional, standalone degree-4 local polynomial derivative estimator: five contiguous valid donors containing the target, centered first, then target positions `[1,3,0,4]`. No gap bridging, extrapolated query cells, source-mask changes, low-temperature bug repair, or runtime fallback integration.

```bash
PYTHONDONTWRITEBYTECODE=1 python3 mesa/eos/eosFreeEOS_builder/biquintic/edge_control/test_edge_control.py -v
PYTHONDONTWRITEBYTECODE=1 python3 mesa/eos/eosFreeEOS_builder/biquintic/edge_control/edge_control.py
```

Outputs exclusively `freeeos_biquintic_progress/experiments/020_edge_control/`. Default is the total-potential control; optional `--split-radiation` writes a separate `split_radiation/` subdirectory and uses the existing declared analytic-reference APIs, not a verified material EOS or historical radiation convention. The recorded primary experiment is total-potential only. Eight OpenBLAS threads; no runtime probe, debugger, dependency install, approval, commit, or push.

The six low jets and native/finite/positive/density/logT>=3.2 source mask remain core-generated. Missing 21/12 jets use both redundant derivative routes; 22 uses both second-stage routes. Window starts, target positions, stage masks, and recursive conservative source-support bounds are archived. Source donors below 3.2 and invalid rows cannot be used; short valid runs remain rejected.

Coverage classifications reuse experiment019's exact rejected-query/default-MESA fraction archive. Those fractions are not rerun or assigned to this candidate. The sidecar checks core overlap equality, recovered-query positivity/value/Jacobian distributions, strictly withheld stride2 low-edge material nodes against their printed values, and exact shared faces near logT3.2–3.36. Holdouts are table fidelity, not independent source accuracy.

The globally even stride2 control has first valid donor logT3.22; withheld logT3.2 rows remain excluded instead of extrapolating across an invalid coarse cell. Face populations separately distinguish newly covered from centered-overlap faces but do not establish active runtime fractions. Recovered-state magnitude distributions include signed extrema and coordinates. Holdout value departures normalize by absolute printed values; face Jacobian jumps normalize by the maximum of both one-sided derivative magnitudes and response magnitudes. Jacobians use physical natural-log density/temperature coordinates. Gamma3 references the printed Gamma3, not Gamma3-minus1.

Current evidence recovers all 9,216 archived active rejects with clean source cells, including 8,496 in logT3.2–3.28; the 15,952 active rejects below3.2 deliberately remain excluded. This is not a positivity guarantee: broader recovered queries include negative Cv/Cp/compressibility, while sampled active recoveries do not. Coarse holdouts and near-edge material faces have large response tails even when basic values agree closely. See the dated, hashed numerical record rather than promoting this estimator as a safe runtime fallback.

## Actual-MESA classification appendix

`active_quality.py` reconstructs the identical main-run accepted material holdouts and faces, verifies all unclassified numerical metrics exactly, and classifies the combined 3,399 queries in one actual-MESA probe batch. It redirects the existing cache-safe probe launcher in memory to experiment020-owned caches/metadata; no core edits, source builds, or experiment014 output writes. The `SUMMARY.md` appendix and `active_quality.json` distinguish active/pure/blended/inactive populations and explicitly compute the separate `run_audit.own_use` predicate. The original 584 new material holdouts are not an own-use-prefiltered population; the q4.52/T3.56 worst Cv departure has zero actual FreeEOS fraction.

```bash
PYTHONDONTWRITEBYTECODE=1 python3 mesa/eos/eosFreeEOS_builder/biquintic/edge_control/test_active_quality.py -v
PYTHONDONTWRITEBYTECODE=1 python3 mesa/eos/eosFreeEOS_builder/biquintic/edge_control/active_quality.py
PYTHONDONTWRITEBYTECODE=1 python3 mesa/eos/eosFreeEOS_builder/biquintic/edge_control/active_quality.py --reuse-probe
```

The final command reanalyzes the saved identical probe batch without another MESA call, checking original read-only inputs and recorded coordinates/output. Original probe-driver and current analysis-driver hashes are distinguished. Raw probe output and queries are archived under `active_quality_*`. No fractions are invented for failed probes. No independent-accuracy or causality claim follows from these unlike new/overlap populations.

## Explicit 3.1 clean-floor control

`clean_floor.py --source-floor 3.1` uses the lead's optional source-constraint floor API, without changing its default3.2. Outputs exclusively `freeeos_biquintic_progress/experiments/027_clean_floor/`; an optional split-reference run has its own `split_radiation/` directory. The recorded run is total-potential only. Native/finite/positive/flag/density exclusions are unchanged, and all recursive donor support remains at or above the explicitly tested floor. Values below3.1 stay excluded, with no bug repair.

The runner checks native material printed-row closure over3.1–3.2, globally even coarse-only fits versus odd/odd withheld material rows over3.1–3.4, and exact material faces selected near3.1–3.22. One batch classifies all newly accepted archived queries plus these diagnostics by actual baseline-MESA fractions. `controls.npz` records masks, exact queries, both face cells, stencils, source bounds, tested-floor metadata and raw probe rows. `027_clean_floor/INTERPRETATION.md` distinguishes fresh classifications from the inherited019 remaining-rejection census and gives concrete quality limits. Lowering the experimental floor changes previously shifted edge twists; only the original centered complete jets are verified bitwise unchanged.

```bash
PYTHONDONTWRITEBYTECODE=1 python3 mesa/eos/eosFreeEOS_builder/biquintic/edge_control/test_clean_floor.py -v
PYTHONDONTWRITEBYTECODE=1 python3 mesa/eos/eosFreeEOS_builder/biquintic/edge_control/clean_floor.py --source-floor 3.1
```

One-sided quartic derivative weights have larger noise/truncation sensitivity. Coverage is not accuracy or safe fallback behavior; C2 shared jets still permit response-Jacobian jumps. Read `020_edge_control/SUMMARY.md` and `results.json` for failures, residual exclusions, and numerical scope.
