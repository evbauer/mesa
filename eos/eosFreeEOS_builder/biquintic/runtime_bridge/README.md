# Opt-in FreeEOS component/API experiment

This lead-owned harness tests a specific FreeEOS component callback against
the complete fixed-composition adapter in `../runtime_adapter/`. It is not
enabled in production defaults. The legacy component still supplies auxiliary
fields and explicitly hybrid abundance derivatives; the callback replaces
thermodynamic fields before the existing patchwork blends them.

`use_other_FreeEOS_results` defaults false. When enabled, the EOS handle must
provide `other_FreeEOS_results` and set `okay_to_convert_ierr_to_skip=false`.
Missing callback (`-7901`) or unsafe error-to-skip configuration (`-7902`) is
rejected at both component and final density/temperature API entry, before
the selector can turn it into a different EOS route. Errors from the supplied
callback then propagate through the strict selector. The original path and
selection/blending rules are unchanged when the option is off.

The callback receives the full composition tuple and all 27 legacy value,
state-partial and internal abundance-partial rows. It is called only after
successful legacy FreeEOS evaluation and phase/fraction setup. Failed calls
must not be interpreted as valid returned arrays. FreeEOS identity is retained;
no missing-composition-partial repair or composition-consistency claim follows.

The new procedure pointer and control are appended to `EoS_General_Info`,
not inserted among its existing members. Consumers must nevertheless use
matching module/library builds. The experiment builds EOS in a private
directory against existing dependencies; it does not replace installed
libraries, regenerate installed caches, or extract MESA data archives.

Acceptance remains pending until the component/final APIs, malformed inputs,
strict statuses, mapping/derivatives, auxiliary/hybrid policy, opt-out parity,
blends and controlled consumers have been checked. A compiled callback alone
is not runtime or stellar acceptance.
