"""Command-line interface: ``python -m wdm_av_builder synth|build|audit|check``.

* ``synth``  write synthetic source-format ASCII tables for a list of isotopes.
* ``build``  read a source table, estimate nodal data, write processed HDF5.
* ``audit``  run the consistency checks on a source table and print a summary.
* ``check``  re-validate a processed HDF5 file against FORMAT.md.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from . import __version__
from . import io as wio
from .audit import PROCESSED_AUDIT_KEYS, format_summary, run_audit
from .derivs import build_nodal_data
from .synthetic import PROJECT_ISOTOPES, Isotope, make_source_columns


def _grid(args) -> tuple[np.ndarray, np.ndarray]:
    nrho = round((args.logrho_max - args.logrho_min) * args.ppd) + 1
    nt = round((args.logt_max - args.logt_min) * args.ppd) + 1
    logrho = np.linspace(args.logrho_min, args.logrho_max, nrho)
    logt = np.linspace(args.logt_min, args.logt_max, nt)
    return logrho, logt


def _resolve_isotopes(names: list[str] | None) -> list[Isotope]:
    by_name = {iso.name: iso for iso in PROJECT_ISOTOPES}
    if not names:
        return list(PROJECT_ISOTOPES)
    out = []
    for n in names:
        if ":" in n:  # name:Z:A
            nm, z, a = n.split(":")
            out.append(Isotope(nm, int(z), float(a)))
        elif n in by_name:
            out.append(by_name[n])
        else:
            raise SystemExit(
                f"unknown isotope {n!r}; use name:Z:A or one of {sorted(by_name)}"
            )
    return out


def cmd_synth(args) -> int:
    logrho, logt = _grid(args)
    out_dir = Path(args.out_dir)
    for iso in _resolve_isotopes(args.isotopes):
        cols = make_source_columns(
            logrho,
            logt,
            iso.Z,
            iso.A,
            radiation=args.radiation,
            with_second_derivs=args.second_derivs,
            with_zbar=args.zbar,
        )
        cols.pop("_exact", None)
        header = {
            "isotope": iso.name,
            "Z": iso.Z,
            "A": iso.A,
            "element": iso.element,
            "includes_radiation": args.radiation,
            "energy_zero": "ideal: ionized ideal ions + nondeg electrons, no rest mass",
            "entropy_zero": "ideal classical Sackur-Tetrode",
            "model": "synthetic analytic ideal EOS (wdm_av_builder.synthetic)",
            "constants": "MESA const_def.f90",
            "tolerance": 1e-12,
        }
        table = wio.SourceTable(header=header, logrho=logrho, logt=logt, columns=cols)
        path = out_dir / f"eosWDM_AV_{iso.name}.txt"
        wio.write_source_ascii(table, path)
        print(f"wrote {path}  ({table.nlogt} x {table.nlogrho})")
    return 0


def build_one(
    source_path: Path,
    out_path: Path,
    *,
    order: int = 5,
    fxt_route: int = 1,
    zero_twist: bool = False,
    e_offset: float = 0.0,
    f_offset: float = 0.0,
    ignore_source_derivs: bool = False,
    write_audit: bool = True,
) -> dict:
    """End-to-end: read source, build nodal data, (audit,) write processed."""
    source = wio.read_source(source_path)
    nodal = build_nodal_data(
        source,
        order=order,
        fxt_route=fxt_route,
        zero_twist=zero_twist,
        e_offset=e_offset,
        f_offset=f_offset,
        ignore_source_derivs=ignore_source_derivs,
    )
    attrs = dict(nodal.attrs)
    attrs.update(
        wio.default_provenance(
            source_file=str(source_path),
            source_model=str(source.header.get("model", "")),
        )
    )
    audit_maps = None
    summary = None
    if write_audit:
        maps, summary = run_audit(source, nodal)
        audit_maps = {k: maps[k] for k in PROCESSED_AUDIT_KEYS if k in maps}
    wio.write_processed_hdf5(
        out_path, attrs=attrs, datasets=nodal.datasets, audit=audit_maps
    )
    return {"attrs": attrs, "summary": summary}


def cmd_build(args) -> int:
    res = build_one(
        Path(args.source),
        Path(args.out),
        order=args.order,
        fxt_route=args.fxt_route,
        zero_twist=args.zero_twist,
        e_offset=args.e_offset,
        f_offset=args.f_offset,
        ignore_source_derivs=args.ignore_source_derivs,
        write_audit=not args.no_audit,
    )
    print(
        f"wrote {args.out}  (order {args.order}, "
        f"fd_scheme: {res['attrs']['fd_scheme']})"
    )
    if res["summary"] is not None:
        print(format_summary(res["attrs"]["isotope"], res["summary"]))
    return 0


def cmd_audit(args) -> int:
    source = wio.read_source(args.source)
    nodal = build_nodal_data(source, order=args.order)
    _, summary = run_audit(source, nodal)
    print(format_summary(source.name, summary))
    if args.plot:
        from .plots import plot_audit_maps

        p = plot_audit_maps(source, args.plot, nodal)
        print(f"wrote {p}")
    return 0


def cmd_check(args) -> int:
    problems = wio.validate_processed(args.processed)
    if problems:
        print(f"{args.processed}: NONCONFORMING")
        for p in problems:
            print(f"  - {p}")
        return 1
    print(f"{args.processed}: conforms to FORMAT.md version {wio.FORMAT_VERSION}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="python -m wdm_av_builder")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    ps = sub.add_parser("synth", help="write synthetic source ASCII tables")
    ps.add_argument("--out-dir", required=True)
    ps.add_argument("--isotopes", nargs="*", help="names or name:Z:A (default: all)")
    ps.add_argument("--logrho-min", type=float, default=-10.0)
    ps.add_argument("--logrho-max", type=float, default=6.0)
    ps.add_argument("--logt-min", type=float, default=3.0)
    ps.add_argument("--logt-max", type=float, default=9.0)
    ps.add_argument("--ppd", type=float, default=20.0, help="points per decade")
    ps.add_argument("--radiation", action="store_true")
    ps.add_argument(
        "--second-derivs",
        action="store_true",
        help="include dPdRho_T, dPdT_Rho, Cv columns",
    )
    ps.add_argument("--zbar", action="store_true", help="include Zbar=Z column")
    ps.set_defaults(func=cmd_synth)

    pb = sub.add_parser("build", help="build a processed HDF5 table")
    pb.add_argument("source")
    pb.add_argument("out")
    pb.add_argument("--order", type=int, choices=(3, 5), default=5)
    pb.add_argument("--fxt-route", type=int, choices=(1, 2, 3), default=1)
    pb.add_argument("--zero-twist", action="store_true")
    pb.add_argument(
        "--ignore-source-derivs",
        action="store_true",
        help="ignore supplied dPdRho_T/dPdT_Rho/Cv; use finite differences for "
        "the second derivatives (second_derivs_from_source=0)",
    )
    pb.add_argument("--e-offset", type=float, default=0.0)
    pb.add_argument("--f-offset", type=float, default=0.0)
    pb.add_argument("--no-audit", action="store_true")
    pb.set_defaults(func=cmd_build)

    pa = sub.add_parser("audit", help="audit a source table")
    pa.add_argument("source")
    pa.add_argument("--order", type=int, choices=(3, 5), default=5)
    pa.add_argument("--plot", help="write audit-map PNG to this path")
    pa.set_defaults(func=cmd_audit)

    pc = sub.add_parser("check", help="validate a processed HDF5 file")
    pc.add_argument("processed")
    pc.set_defaults(func=cmd_check)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
