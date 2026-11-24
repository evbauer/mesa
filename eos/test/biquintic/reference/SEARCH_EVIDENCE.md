# Saved reference-dependency evidence

## Scope and disposition

Recorded 2026-10-28 from the initial local inspection and brief resumed network
attempt; no new broad search was performed to complete this document.
Current scope is **existing tabulations only**. The compile-only direct-call
scaffold is deferred and unlinked.

No FreeEOS archive, source tree, library, executable, or network response body
was obtained. No tables were generated and no software was installed.
The reference inventory contains the scaffold, its object, three validation
logs, and documentation only. After the scope update,
`ps -C curl -o pid,etime,args` returned only its header: no active download.
The interrupted escalated release-list request had no execution session ID.

## Instructions and repository

Initial checks covered `AGENTS.md` at `/g`, `/g/g14`, `/g/g14/bauer39`,
`/g/g14/bauer39/work`, `/g/g14/bauer39/work/WD_Gbar`, and the workspace root,
and `find mesa -name AGENTS.md -print`; no applicable files were found.
Later scoped checks under `eosFreeEOS_builder/biquintic` also returned none.

The workspace root is not a Git repository. `git -C mesa rev-parse HEAD`
returned `7b19fe283fc5315a7aaeb02f353e5ee4e611f3ed`.
Pre-existing plotter edits and the builder plan were not changed.

## Local searches completed before the pause

These are historical commands/results, not instructions to repeat expensive
searches in the table-only task. Case-insensitive patterns were `*freeeos*`,
`*free_eos*`, `*free-eos*`, and `*mesasdk*` as applicable.

| Command / locations | Result |
| --- | --- |
| `find . -iname '*freeeos*' -o -iname '*free_eos*'` | Builder, archives, installed tables/caches, interpolation artifacts, reports, branch metadata; no solver dependency located. |
| `find /g/g14/bauer39 -path /g/g14/bauer39/work/WD_Gbar/evbauer -prune -o \( -iname '*freeeos*' -o -iname '*free_eos*' -o -iname '*free-eos*' -o -iname '*mesasdk*' \) -print` | No matching names outside workspace. |
| `find /opt /usr/local /usr/lib64 /usr/lib /usr/include /usr/tce/packages -iname '*freeeos*' -o -iname '*free_eos*' -o -iname '*free-eos*' -o -iname '*mesasdk*'` | Completed, exit 1 from permission errors; no candidate printed. Default find does not recursively follow directory symlinks. |
| `find /usr/tce/packages -maxdepth 4` with those four name patterns | No candidates. |
| `find /usr/src /usr/local/src -maxdepth 4` with the three FreeEOS patterns | No candidates. |
| `find /usr/lib64/pkgconfig /usr/share/pkgconfig /usr/local -iname '*eos*' -o -iname '*mesa*'` | Only `/usr/lib64/pkgconfig/osmesa.pc`, not FreeEOS. |
| `find mesa/build -name '*.a' -print0 \| xargs -0 nm -A --defined-only 2>/dev/null \| rg -i ' [Tt] free_eos_?$'` | No direct solver text definition in inspected archives. |
| `nm -A --defined-only mesa/build/eos/lib/libeos.a`, filtered for FreeEOS | MESA metadata/control symbols and `free_eos_handle`, not the solver. |
| Fortran search for `subroutine free_eos`, `function free_eos`, `call free_eos` | Only builder external call at `src/free_eos_table.f90:124`; no implementation/declaration found. |

Inaccessible paths reported: GCC-toolset private `root/root` directories
(10–15), GCC-toolset-10 `var/lib/stap-server`, `/opt/intel/oneapi/modulefiles`,
`/opt/jupyterhub-batchspawner`, `/opt/FalconForensics`, `/usr/lib64/redis`,
and Nagios event-handler directories `redundancy-scenario1` and
`distributed-monitoring`. This is not a filesystem-complete absence proof.

Absent inspected paths: `/g/g14/bauer39/mesa`, `/g/g14/bauer39/mesasdk`,
`/g/g14/bauer39/.local`, `/opt/mesa`, `/opt/mesasdk`, `/usr/local/mesa`,
`/usr/local/mesasdk`, `/usr/tce/packages/mesasdk`,
`mesa/build/free-eos-builder`, and legacy `mesa/lib`.
`/usr/tce/packages/mesa` has `mesa-19.0.1`/`mesa-23.3.6` directories but the
bounded inspection found no FreeEOS candidate. Their names do not establish
the scientific FreeEOS solver's presence.

## Environment, build metadata, and checks

- `MESASDK_ROOT`, `MESA_DIR`, `FREE_EOS_DIR`, `PKG_CONFIG_PATH` were unset;
  `LD_LIBRARY_PATH` only named `/usr/tce/packages/texlive/texlive-20220321/lib`.
- `pkg-config --variable pc_path pkg-config` returned
  `/usr/lib64/pkgconfig:/usr/share/pkgconfig`.
- `pkg-config --modversion free-eos` and `--libs free-eos` both reported the
  package missing. `ldconfig -p` filtered for `free.?eos` returned no entry.
- Builder Makefile declares external dependency `free-eos`; MESA
  `make/setup-depends.mk` uses pkg-config. The legacy README's `make/makefile`
  instruction does not match this layout. Installed `mesa-eos.pc` lists MESA
  internal dependencies, not a direct solver.
- GNU Fortran 10.3.1 20210422 (Red Hat 10.3.1-1); nm/ar/make/python3 available.
- `make check`: compile-only success with `-O0 -g -Wall -Wextra -fcheck=all`,
  no diagnostics; `logs/compile_check.log`.
- `gfortran build/direct_freeeos.o -o build/dependency_probe`: exit 1,
  `undefined reference to free_eos_`; `logs/dependency_probe.log`.
- Default `make`: refuses missing `FREEEOS_LIBS`, make exit 2;
  `logs/dependency_guard.log`. No direct query or source-accuracy test ran.

Build/run/dependency directories are ignored locally; validation logs remain
tracked evidence. No generated tables/caches exist here.

## Resumed network attempt, stopped

| Request | Outcome |
| --- | --- |
| Normal `curl -IL --connect-timeout 15 --max-time 40 https://freeeos.sourceforge.net/` | Exit 6: could not resolve host. |
| Same with `require_escalated` approval request | Exit 6: could not resolve host. |
| Normal HEAD request to `https://sourceforge.net/projects/freeeos/files/` | Exit 6: could not resolve host. |
| Same release-list request with `require_escalated` | Interrupted by user; no response/session ID. |

No alternate host/IP/proxy bypass was attempted. No official release filename,
hash, ABI success semantics, or upstream radiation/gauge conventions were
confirmed. FreeEOS 2.2.1 EOS1 is the existing documentation's table-family
attribution, not retrieved release provenance. Download/build pursuit stops.
