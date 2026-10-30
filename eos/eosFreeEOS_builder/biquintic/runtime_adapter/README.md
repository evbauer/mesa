# Fixed-composition FreeEOS runtime adapter

This isolated module fills the full MESA27-slot component return using the
unchanged standalone biquintic cache/evaluator. It owns no routing, EOS
handles, global options, fallback, table generation, or installed libraries.
The integration interface is in experiment047 `API.md`; the accompanying
`solar.identity` pins the existing experiment008 solar split cache.

## Initialize and query

After the caller initializes MESA const/math/chem, populate `runtime_identity`
from a trusted manifest and call
`runtime_initialize(adapter, cache_path, identity, pattern_path, ierr)`.
Publish a successfully initialized adapter to worker threads only after
initialization returns. Never reinitialize that object concurrently with a
query. A failed initialization leaves it uninitialized.

The adapter computes the actual cache SHA256 before and after loading, then
checks the complete declared source/preprocessing/grid/composition/reference
metadata against the trusted identity. The byte digest covers every node,
validity mask, and trailing byte. The existing V1 header alone is insufficient.
The source/preprocessing digests identify the audited source and algorithm;
the caller must supply the expected identity from that audit, not derive it
from the cache under examination. Inputs must be stable during initialization;
the two digest reads detect ordinary concurrent changes but are not an
adversarial filesystem snapshot protocol. Successful queries do no I/O or
hashing and use immutable, privately held loaded arrays and composition data.

Only radius2/mean/radiation-split/logT-floor3.2 is accepted. The cache radiation
constant must agree with installed MESA within32 machine epsilons; the full
loaded metadata must equal the trusted header exactly. The floor is explicitly
guarded on every query; its preprocessing/mask provenance is covered by the
trusted algorithm and payload identities. There is no cold-table repair.

`runtime_query` receives the complete legacy component result and full
composition tuple. It replaces the13 potential-derived values and both state
partials, including all six partials of `dE_dRho,dS_dT,dS_dRho`. Pressure,
energy, and entropy API values are natural logarithms; pressure subtracts
MESA radiation, while energy, entropy, and response fields use the restored
total potential. Gamma3 is Gamma3 itself. Phase/latent fields and their state
partials become zero; FreeEOS fraction becomes one, all other fractions zero.

The actual legacy `mu,lnfree_e,eta` values and both state partials are retained.
All27-by20 legacy abundance rows must be finite and are retained exactly.
`runtime_hybrid_policy` explicitly tags this policy; the fixed-composition
potential has no source composition derivatives. This is not a
composition-consistent EOS and must not be used for a normal burning run.
The eventual fixed-composition consumer must prove that it ignores abundance
derivatives; disabling a control alone is not such proof.

Every query checks the exact audited20-isotope proxy at X=.7,Z=.02. Canonical
IDs are looked up once in initialized MESA chemistry, using the unchanged
audited solar metal file and its SHA256. Metals are normalized once using
the existing probe procedure; queries never renormalize. Every mass must be
positive and agree within5e-13 relative, with total normalization within5e-13
absolute. All four supplied summaries must agree within5e-13 relative with
the vector recomputed from captured isotope mass numbers and charges. Array
lengths, all forward/reverse mappings, duplicate IDs, and nonzero pattern are
checked. Network ordering may change consistently; extra zero species are
rejected. Inconsistent rho/logRho or T/logT is rejected at64 machine epsilons
in decimal log coordinates. The adapter adds a strict outer domain guard so
the underlying evaluator's endpoint-roundoff allowance never clips a query.

All output arrays remain unchanged after any failed query. On success every
slot is finite. Callers must check status before using arrays. No fallback or
counter mutation occurs. Underlying table errors10--18 are propagated;
adapter errors are100 uninitialized,101 identity,102 composition,103 input,
104 uninitialized/nonfinite legacy data,105 unsupported policy. Error104 can
verify finite contents, not prove the provenance of arbitrary finite inputs;
the caller remains responsible for obtaining the actual legacy result first.

## Bounded private validation

Run only on the inherited allocated compute node, with two library threads
and8GiB virtual memory. From the workspace:

```sh
bash freeeos_biquintic_progress/experiments/047_runtime_adapter/run_validation.sh
```

The Makefile compiles unchanged `helm_polynomials`, `freeeos_biquintic`, and
`freeeos_table` into a private build directory, then the adapter/test driver.
It links installed const/chem/math/utils/auto_diff dependencies and does not
build or install EOS. Synthetic caches are analytic test fixtures, not new
FreeEOS source tables. Numerical outputs, command transcripts, identities,
and results stay in experiment047. Component/final-API acceptance is owned
by the lead's separate runtime bridge and is not supplied by these unit tests.
