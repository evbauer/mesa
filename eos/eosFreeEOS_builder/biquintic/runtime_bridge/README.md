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

Experiment050 passes the component/final mapping, strict-error, hybrid,
opt-out, stronger-blend and derivative checks. Eight additional transport-stub
checks harden the Python client against wrong callback/fraction relations and
overflow in reconstructed physical states. Three additional pipe checks cover
successful lines, EOF and the60-second per-query deadline (20 cases total).
The client is serial; it is not the stellar callback implementation.

051 repeats the bounded consumers through both APIs with original forcing:
2542/2544 converge, with the same candidate scalar1e-13stress case failing
twice. The standard-tolerance results survive; an independent precision
diagnostic locates the exception within a cell. It is not silently replaced
by a successful easier-tolerance result. Candidate stellar testing requires
separate coherent-build, global-composition and hybrid-irrelevance gates.
Direct adapter malformed-input checks do not sanitize public-MESA pointers
before its own composition processing. No unrestricted acceptance is claimed.
