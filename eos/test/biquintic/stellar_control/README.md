# Frozen-composition stellar control (TASK049)

This experiment prepares a short **LEGACY-only** stellar control against the installed, mutually compatible MESA modules and libraries. Candidate evolution is blocked in both the Python launcher and the Fortran driver. A successful legacy run is not an API acceptance result or a candidate stellar comparison.

## Physical route and exact composition

The starting structural seed is the existing 366-zone, one-solar-mass, age-1-Myr pre-main-sequence model in `star/test_suite/simplex_solar_calibration/standard_late_pre_zams_1.0M.mod`. Its nuclear luminosity is only about 3e-5 Lsun, so a short nonburning contraction test has a physical interpretation. It is not a main-sequence model with an artificially disabled dominant energy source. Its original eight-isotope mixture is not an acceptable candidate composition.

`harness.py prepare` makes a traceable derived saved model: it preserves the structural columns and replaces the network and every cell's abundance vector with the exact audited pattern. Hydrogen and helium are 0.70 and 0.28; the 18 metal weights are multiplied by `0.02/sum(weights)` once. The original metal file is authenticated by SHA256. The local `frozen20.net` uses `add_isos` only, with no reactions:

```
h1 he4 c12 n14 o16 ne20 na23 mg24 al27 si28
p31 s32 cl35 ar40 ca40 ti48 cr52 mn55 fe56 ni58
```

This is an EOS proxy mixture, including Ar40, not a general scaled-solar stellar network. No extra species, even at zero abundance, are accepted. A derived seed must first pass legacy structural adjustment. It is not an equilibrium model merely because its structural columns came from a stellar model.

## Frozen residual policy

- `dxdt_nuc_factor=0`, `eps_nuc_factor=0`, and `mix_factor=0`; reaction-free network; no adaptive or automatic network extension.
- No diffusion, convective premixing, predictive mixing, remeshing, accretion/winds, or rotation. Convection still transports heat through the same MLT physics.
- `energy_eqn_option='eps_grav'`, `use_time_centered_eps_grav=true`, `include_composition_in_eps_grav=false`, `fix_d_eos_dxa_partials=false`; WDM_AV off and strict EOS error handling.
- Source inspection establishes that `struct_burn_mix` chooses `nvar_hydro`, excluding all abundance unknowns, when both burning and mixing flags are false. Runtime matrix hooks check this directly. `star_utils` drops abundance-column entries beyond the active variable count.
- The final EOS callback can replace **all** returned abundance derivative rows by `+/-1e6*(row+species_index)` while retaining all values and density/temperature derivatives. Matched replay and full raw residual/Jacobian-block output test whether those placeholders influence the intended solve. This is a dependency test, not a claim that the placeholders are physically correct.

## Runtime audit

The installed final-EOS callback checks the complete isotope vector and forward/reverse mappings at every successful density/temperature query, including EOS inversion subqueries. It rejects nonfinite/negative masses, normalization errors, wrong species identity, extra species, or drift greater than 5e-13. It checks WDM_AV is inactive and records query coordinates, phase, model/attempt and FreeEOS fraction. Candidate integration must additionally retain lead/047's **pre-return component guard of the full callback tuple**; this final callback does not replace it.

The 049 threshold is absolute; 047 uses 5e-13 relative to each pinned isotope. The completed stellar runs have zero drift from the stored049 vector. A matching-SDK Fortran audit of 047's normalization finds maximum relative difference 2.1746673870457e-16 between that stored vector and its reference, safely below 047's criterion. Preserve the stricter component guard in a candidate integration.

`other_before_struct_burn_mix` runs before `save_start_values`; it checks current composition and tags each attempt. At the first matrix callback of each attempt, the saved previous-state vector is available. All current and saved previous compositions are checked there; previous rho/T/composition queries are explicitly replayed through the final API and tagged separately. The composition term in the actual residual is disabled; these replay queries are diagnostic, not an additional energy term. Accepted-step and startup checks are also logged. Rejected trial states remain in the query log.

`iterations.tsv` records correction and residual norms, `events.tsv` records attempts/retries/timestep/age, `matrices.bin` stores full active residual and three Jacobian blocks before factorization, and MESA history/profiles and photos are saved every accepted step. Matrix records contain four native int32 values `(model, attempt, nvar_hydro, nz)`, followed by Fortran-order float64 arrays `equ`, `dblk`, `lblk`, `ublk`. This local binary artifact requires the recorded compiler/platform; the summaries are portable JSON.

## Reproduction and limits

Run on inherited allocation 5528075/rzwhippet14, before the original 21:20 PDT hard stop. The script enforces placement and time. The lead cleared an isolated build with at most eight compile processes; numerical work uses one thread and an 8 GiB address-space limit. Installed dependencies and the protected plotter files are hashed before and after builds/runs. No shared module, library, source table or existing test case is modified. Runtime caches are private copies of existing binary caches.

From the workspace root:

```sh
python3 mesa/eos/test/biquintic/stellar_control/harness.py build
python3 mesa/eos/test/biquintic/stellar_control/harness.py prepare legacy_control
python3 mesa/eos/test/biquintic/stellar_control/harness.py run legacy_control
```

Use a new run name every time. `--abundance-mode 1` or `-1` on `run` selects the placeholder stress test. The baseline is bounded to eight accepted models, 100-yr initial/1000-yr maximum timesteps, 12 retries and a 15-minute external timeout. Logs, photos, saved models and command/provenance JSON remain under experiment 049. An external timeout can interrupt an incomplete step; the latest completed photo/model remains the restart artifact.

The completed handoff is `control_legacy/seed.mod`, SHA256 `cd5b74e7ccd55d8253d045abc5bb95cba1506ba42fba70482bccf5f29cf9eb3e`, alongside its six matched inputs. `prepare NAME --seed PATH` creates the derived copy from a warmed model. `run NAME --force-retry` requests the same single deliberate model2 retry used by all three matched controls. `analyze.py --compare control_legacy control_plus control_minus` checks the complete residual/Jacobian streams, physical outputs and input/executable identities. `test_guards.py --output-name guard_checks_relative` exercises seven actual installed-API cases and the normalization comparison; use another fresh output name to repeat it.

Final evidence: the eight-step warmup has zero retries; three eight-step derivative controls each have one deliberate retry and two Newton iterations per accepted step. Their matrices, histories, all profiles and final models are byte-identical. All logged FreeEOS states have candidate mask support. Candidate evolution remains unperformed.

## Candidate integration point and gates

The lead owns the component callback and adapter. Replace the explicitly closed backend installation branch in `extras_controls` only in a coordinated integration, using a consumer dependency closure rebuilt against the new EOS interface. Appended derived-type fields change object size/array stride: an old star library plus a newly linked EOS library is not assumed ABI safe.

Before enabling candidate evolution, require explicit lead API clearance, compatible builds, authenticated cache/mixture identities, strict status propagation, final routing/blend/fallback accounting, and a common composition-validated legacy-adjusted seed. Preserve the independent final query guard/monitor. Repeat abundance-placeholder invariance using the candidate path as well. Every potential-routed current/previous/trial state must have complete centered source-floor-3.2 cache support; a stellar photosphere/envelope may violate support even when interior zones are safe. Do not introduce silent legacy fallback, temperature clipping or selector changes to obtain a run. A bounded nonzero fallback policy would require a separately declared lead decision and separate populations.

Only matched runs from the same adjusted seed, with identical physics, timestep policies and tolerances, qualify as the eventual stellar A/B. The preparation/legacy adjustment is distinct. A run that does not reach the previously identified 24 solver faces cannot settle their convergence impact.
