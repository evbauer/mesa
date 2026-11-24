# TASK057: matched-age temporal refinement

Explicitly cleared at19:48PDT on29October2026, exactly as the frozen052 `TASK057_PROPOSAL.md`. This is a temporal-sensitivity experiment at frozen composition, not solver retuning, burning or a continuum-order claim. Four sequential runs only: coarse baseline/candidate, then half baseline/candidate after both coarse reproduction gates pass. No further resolution or automatic rerun is authorized.

The executable, its Fortran guard/callback driver, potential kernel and coherent private libraries are read-only052artifacts. No build occurs. `timestep.py` imports frozen052Python helpers read-only with bytecode writing disabled; it verifies the551-artifact freeze at preparation. All new outputs are in `057_stellar_timestep/`; no052source or evidence is changed.

## Inputs and endpoint

Every run starts from the exact049common seed SHA256 `cd5b74e7ccd55d8253d045abc5bb95cba1506ba42fba70482bccf5f29cf9eb3e`, using052binary SHA256 `70c1d3f50d36c22519093675438aead92a4eeb18d4602f93a4d280cc240871df`. The only four whitelisted time/stop keys are `max_model_number=300`, `max_age=1089056.4362468508d0`, initial dt100/50yr and maximum dt1000/500yr for coarse/half. Other input files match052byte-for-byte. At each resolution, backend inputs match exactly; only the hook-selector environment value differs.

Requested successful termination is `max_age`, with the final accepted event and saved profile at or above the target and within2binary64ULPs. Model-limit, minimum-dt, timeout, raw-cap or other termination is a failure. No endpoint interpolation or control retuning is permitted. All frozen20composition, full27-row hybrid policy, current/previous/retry/final guards, strict support/no fallback, physics and nonlinear tolerances are unchanged. The same deliberate model2retry is retained.

`timestep_audit.py` checks streamed component/final sequence/state/phase, normalization and relative composition, fractions/WDM, completion counters, matrices, endpoint, iterations/retries and accepted-profile native crossings. Before either half run, coarse pre-endpoint event/profile/matrix streams must reproduce052exactly. Final endpoint clamping may differ only at the declared floating-point scale: final dt within2target-age ULPs in seconds and final selected-profile scaled differences at most5e-11, with matching iterations. These are external fail-closed audit limits, not modifications to solver tolerances.

## Resources and commands

Run only on inherited rzwhippet14/job5528075. One numerical thread and8GiB per sequential run;300s timeout,2GiB raw-log cap excluding copied immutable caches.28GiB disk reserve is required at preparation. No new run after20:15; stellar work stops20:20. Failed attempts and STOP.json remain evidence, and block further launches.

From workspace root, with `PYTHONDONTWRITEBYTECODE=1`, all numerical-library thread caps1, and `ulimit -v 8388608`:

```
python3 mesa/eos/test/biquintic/stellar_timestep/timestep.py prepare
python3 mesa/eos/test/biquintic/stellar_timestep/timestep.py run coarse_baseline
python3 mesa/eos/test/biquintic/stellar_timestep/timestep.py run coarse_candidate
python3 mesa/eos/test/biquintic/stellar_timestep/timestep.py run half_baseline
python3 mesa/eos/test/biquintic/stellar_timestep/timestep.py run half_candidate
python3 mesa/eos/test/biquintic/stellar_timestep/analyze_refinement.py
```

Each launch requires successful prior audits; half runs additionally require both coarse gates. Never reuse an existing run directory or bypass STOP.json. Per-run metadata/logs retain command, resource caps, exact inputs, before/after hashes, status, all matrices/queries/profiles/photos and final model. The primary endpoint analysis compares within-backend half-minus-coarse shifts, EOS contrasts at coarse and half resolution, and their interaction on the identical mass grid. Any interior time-alignment plot must distinguish exact ages from interpolation; the primary endpoint analysis does not interpolate.

## Explicit continuation dispositions

The original coarse-candidate final-profile screen fails and remains failed: the3.906ms endpoint clamp changes eps_grav by6.593e-9erg/g/s, above the worker-selected5e-11scaled reproduction ceiling. All pre-endpoint records, endpoint age and runtime guards pass. The ceiling was not a solver tolerance or physical-accuracy bound. `STOP.json` and `analysis/coarse_candidate.json` remain pinned and unchanged.

The lead explicitly authorizes the two half runs as a separate exploratory continuation in `CONTINUATION.md/json`. The launcher recognizes only that specific authorization ID, allowed run names and fixed original-stop/audit digests; it does not convert the original audit to a pass or provide a generic bypass.

Half baseline then completes188steps normally, but offline analysis discovers the unchanged MESA default `max_num_profile_models=100` reuses profile filenames. The original `CONTINUATION_STOP_half_baseline.json` remains intact. At20:14 the lead authorizes only index-decoder correction and the remaining candidate, in `DECODER_CONTINUATION.md/json`, with launch cutoff20:25, numerical cutoff20:30 and report20:45. No output-retention input, physical control or executable changes.

`retained_profiles` now reads `profiles.index` and verifies each file's model metadata; final model188 is in `profile89.data`. Half-run models2–89 are not retained, while model1and90–188remain. Full query/matrix/event streams and photos remain. Crossings count only consecutive retained model pairs; the1to90gap is excluded. A saved endpoint is available without rerunning or reconstructing missing profiles. The original failed coarse screen and decoder stop stay in the final provenance.

`analyze_refinement.py` requires the specific reviewed continuation and all unaffected/new runtime gates, then reports endpoint contrasts against the empirical coarse-clamp sensitivity. Two resolutions do not prove convergence; contrasts near that measured scale are conservatively unresolved. New genuine failures stop the remaining sequence in separate artifacts, never overwriting an original failure record.
