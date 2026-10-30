# TASK052 exploratory stellar A/B

This is the explicitly approved eight-step frozen-composition stellar relaxation experiment. It is not a burning calculation, production acceptance or a global rough-face stress test. The two retained tight-1e-13 task051 stress failures are not erased by this experiment. Do not extend these runs without separate lead clearance.

At19:31PDT the lead additionally cleared the named `TASK052_same_seed_100` extension after the first-pair audit passed. It uses separate `baseline_100` and `candidate_100` directories, the exact same common seed and frozen executable, and changes only `max_model_number=100` in both inlists. The original pair is preserved. No further extension is implied.

## Inputs and builds

Both runs copy the exact six files from experiment049 `control_legacy/`, including the preconditioned `seed.mod` with SHA256 `cd5b74e7ccd55d8253d045abc5bb95cba1506ba42fba70482bccf5f29cf9eb3e`. The twenty-isotope reaction-free network, uniform mixture, standard time-centered `eps_grav`, disabled composition-changing processes, timestep policy, eight-model limit, output controls and deliberate model2 retry are unchanged. Both also receive the same authenticated047 `solar.identity` fixture.

`Makefile` uses only the coherent private052dependency tree. Rates, kap, net, atm, star_data and star are privately rebuilt packages supplied by the lead; EOS resolves to the reviewed private050build. Other dependencies are verified links to installed packages. The launcher checks GNU link-map archive paths and hashes before freezing the single executable. Never prepend a new EOS library to the old installed star library.

The frozen047adapter and unchanged potential core sources are compiled into `052_stellar_api/runs/driver_build/`, with no shared/library rebuild or install. The build record and `freeze.json` preserve commands, sources, actual linked archives and all private/installed dependency identities. The build and both runs verify these hashes remain unchanged.

## Callback and composition policy

The driver derives from049's tested `run_star_extras` logic. `extras_controls` initializes the immutable047adapter from the authenticated008cache before stellar model-load EOS queries. Both backend modes perform the same initialization. The sole physics selector is `TASK052_BACKEND=BASELINE` or `CANDIDATE`, which disables/enables `other_FreeEOS_results` on the same EOS handle, before the existing blend. No input file is changed between modes.

The component wrapper calls `runtime_query` without reimplementing potential mathematics. It verifies that all27legacy abundance derivative rows remain finite and unchanged, and that auxiliary mu/lnfree_e/eta values and rho/T partials are exactly preserved. This explicitly hybrid policy relies on049's measured irrelevance to the frozen structural residual; it is not a composition-consistency claim. A nonzero adapter status or policy violation is logged and terminates the run immediately. There is no fallback, clipping, control adjustment or silent legacy return.

The independent final-EOS guard runs even when FreeEOS is inactive. It and every current/previous/accepted-state check enforce the047relative5e-13per-isotope criterion, full species/mapping identity, finite positive masses and absolute normalization tolerance. Its reference is normalized with the same Fortran procedure as047. It does not change any stellar abundance. The component adapter additionally validates the complete X/Z/abar/zbar tuple.

Previous-state replay occurs at the first matrix callback per attempt, after MESA's `save_start_values`. Every matrix callback checks both current and saved composition, active four-variable structural block dimensions and fixed physics. The energy equation's composition term remains disabled; replays are diagnostics only. One deliberate retry at model2 is identical in both modes.

## Artifacts and reproduction

Run locally on inherited job5528075/rzwhippet14. One numerical thread,8GiB address-space cap, up to8driver compilation processes and15-minute run timeouts are enforced. The hard stop remains21:20PDT. Existing EOS binary caches are copied into each run's own cache directory. Only the worker-owned source directory,052`runs/`, `candidate_summary/` and task report paths are written; lead-owned build/checkpoint/support preflight remain read-only.

From workspace root:

```sh
python3 mesa/eos/eosFreeEOS_builder/biquintic/stellar_candidate/harness.py build
python3 mesa/eos/eosFreeEOS_builder/biquintic/stellar_candidate/harness.py prepare baseline_00
python3 mesa/eos/eosFreeEOS_builder/biquintic/stellar_candidate/harness.py prepare candidate_00
python3 mesa/eos/eosFreeEOS_builder/biquintic/stellar_candidate/harness.py run baseline_00 BASELINE
python3 mesa/eos/eosFreeEOS_builder/biquintic/stellar_candidate/harness.py run candidate_00 CANDIDATE
python3 mesa/eos/eosFreeEOS_builder/biquintic/stellar_candidate/analyze.py
```

Existing run directories are never overwritten. Use the archived results for the completed experiment; any new numerical extension needs lead clearance. Set `PYTHONDONTWRITEBYTECODE=1` for analysis, and keep numerical-library thread caps at1.

Each run preserves the command, process identity, before/after dependency hashes, outcome, all successful final-EOS queries/fractions/relative-composition errors, full candidate-component inputs/status/hybrid checks, backend identity record, Newton iteration norms, active residual/Jacobian streams, attempts/retries/timesteps, every-step profiles/photos and final model. A failed candidate attempt remains evidence and must not be silently rerun with altered controls.

`analyze.py` verifies complete run/callback accounting, compares the private opt-out run with049, and measures profile differences on the identical Lagrangian grid. It also counts native-cell changes between accepted profiles1through8 and distances to the frozen24face manifest. Endpoint cell changes give minimum crossings under a continuous path assumption; they do not recover all within-step or rejected-trial crossings, or prove a selected rough face was crossed. If accepted timesteps differ, profile differences at the same model number also include different physical ages.

## Cleared100-step extension

`extend.py` reuses the frozen binary without changing its Fortran driver, adapter, kernel or dependencies. It independently verifies the one permitted inlist change, original input hashes, frozen source hashes and private/installed dependency inventory. Both runs retain one numerical thread,8GiB, the deliberate model2 retry, full guards, profiles/photos and raw audit streams. Each run has a900-second wall limit and a1GiB raw-output monitor cap (copied existing EOS caches excluded); a cap/error ends and preserves the attempt, not a silent retry. Private cache copies require approximately4.49GB per run.

The recorded execution sequence, authorized only for this named pair, is:

```sh
python3 mesa/eos/eosFreeEOS_builder/biquintic/stellar_candidate/extend.py prepare baseline_100
python3 mesa/eos/eosFreeEOS_builder/biquintic/stellar_candidate/extend.py prepare candidate_100
python3 mesa/eos/eosFreeEOS_builder/biquintic/stellar_candidate/extend.py run baseline_100
python3 mesa/eos/eosFreeEOS_builder/biquintic/stellar_candidate/extend.py run candidate_100
python3 mesa/eos/eosFreeEOS_builder/biquintic/stellar_candidate/analyze_extension.py
python3 mesa/eos/eosFreeEOS_builder/biquintic/stellar_candidate/plot_profiles.py --extension
```

Set the same one-thread library caps and8GiB limit for analysis. `analyze_extension.py` checks each component/final record's state, phase and sequence; all stored composition/fraction/lifecycle guards; full input/binary identity; initial eight-profile/matrix-prefix reproduction; iteration/retry/dt accounting; and native-cell transitions. It quantifies common physical-age overlap. If timesteps differ, only exact common ages or explicitly labeled bracketed per-zone temporal interpolation are compared, never extrapolation or an assumption that model numbers imply equal ages. The final profile plot refuses unequal ages.

First-pair analysis remains in `candidate_summary/comparison.json`; the separate extended result is `extension_comparison.json`. `HANDOFF.json` retains the first audit and original proposal, while final handoff text records the subsequent clearance/outcome. The source of `harness.py`, `run_star_extras.f90` and `Makefile` remains pinned by the original `freeze.json` throughout.
