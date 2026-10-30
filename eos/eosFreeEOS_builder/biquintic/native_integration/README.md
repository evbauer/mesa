# Native biquintic FreeEOS integration

This is a default-off experimental backend using the existing FreeEOS source
tables, not a new EOS calculation. Enable through the ordinary `&eos` inlist:

```fortran
FreeEOS_use_biquintic_interpolation = .true.
FreeEOS_biquintic_table_dir = '/absolute/path/to/the/generated/family'
```

The optional directory defaults to `mesa_data_dir/eosFreeEOS_biquintic`.
Once caches exist there, only the logical toggle is necessary. Caches are
generated separately by `../family_cache/`; the backend never modifies source
tables or generates new FreeEOS calculations. Development runs use an explicit
private directory rather than installing or replacing data.

The native backend follows legacy X/Z selection/weights but combines potential
jets first, then derives thermodynamic fields and their composition partials.
Legacy auxiliary mu, free-electron and eta quantities remain in use. This
switch cannot be combined with the earlier `use_other_FreeEOS_results` hook.
Missing/unsupported/nonphysical candidate queries return error -7950 even when
`okay_to_convert_ierr_to_skip` is true; there is no silent legacy fallback.
Detailed underlying status is logged by the family runtime.

Run `bash build_private.sh eos` or `bash build_private.sh all` on the inherited
single compute node. Outputs stay in experiment065. Installed archives are
hash checked before and after builds. The all stage rebuilds EOS and every
direct/indirect stellar consumer of its changed control type: rates, kap, net,
atm, star_data and star. Do not link new eos_def modules with old stellar objects.

The legacy switch-off path should remain unchanged. Whole-domain support,
composition derivative controls, and ordinary stellar evolution are distinct
validation requirements; build success alone does not establish them.
