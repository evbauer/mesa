"""Reading source tables (ASCII per FORMAT.md Section 1, and the HDF5 source
form) and reading/writing processed HDF5 tables (FORMAT.md Section 2).

All I/O keeps the array layout ``(nlogT, nlogRho)`` == ``(ntau, nxi)`` in h5py
(tau/logT is the slow axis), which is exactly what the Fortran/forum reader
expects to receive transposed as ``(nxi, ntau)``.

Base-10 logs live only at the file boundary (source tables and the derived
grid metadata); everything in :mod:`wdm_av_builder.derivs` works in natural
logs xi = ln rho, tau = ln T.
"""

from __future__ import annotations

import datetime
import shlex
from dataclasses import dataclass, field
from pathlib import Path

import h5py
import numpy as np

from . import __version__
from . import constants as C
from .synthetic import ELEMENTS

FORMAT_VERSION = 1

# processed-file dataset/attribute contract (FORMAT.md Section 2)
REQUIRED_DATASETS: tuple[str, ...] = (
    "f",
    "f_xi",
    "f_tau",
    "f_xixi",
    "f_tautau",
    "f_xitau",
    "f_xixitau",
    "f_xitautau",
    "f_xixitautau",
)

INT_ATTRS: tuple[str, ...] = (
    "format_version",
    "Z",
    "interp_order",
    "nxi",
    "ntau",
    "includes_radiation",
    "has_Zbar",
    "has_eta",
    "second_derivs_from_source",
)
DOUBLE_ATTRS: tuple[str, ...] = (
    "A",
    "xi_min",
    "dxi",
    "tau_min",
    "dtau",
    "e_offset",
    "f_offset",
)
STR_ATTRS: tuple[str, ...] = (
    "isotope",
    "element",
    "fd_scheme",
    "source_file",
    "source_model",
    "builder_version",
    "created",
)


# ---------------------------------------------------------------------------
# source tables
# ---------------------------------------------------------------------------
@dataclass
class SourceTable:
    """A single-isotope source table on a uniform log10 (rho, T) grid.

    Attributes
    ----------
    header : dict
        Raw header keys (lower-cased), plus parsed ``Z`` (int), ``A`` (float),
        ``isotope`` (str), ``includes_radiation`` (bool).
    logrho, logt : 1-D arrays of the base-10 log grid.
    columns : dict of 2-D arrays shaped ``(nlogT, nlogRho)`` for each data
        column (P, E, F, S, Zbar, ...).  NaN marks undefined entries.
    """

    header: dict = field(default_factory=dict)
    logrho: np.ndarray = field(default_factory=lambda: np.empty(0))
    logt: np.ndarray = field(default_factory=lambda: np.empty(0))
    columns: dict[str, np.ndarray] = field(default_factory=dict)

    @property
    def nlogrho(self) -> int:
        return self.logrho.size

    @property
    def nlogt(self) -> int:
        return self.logt.size

    @property
    def name(self) -> str:
        return str(self.header["isotope"])

    @property
    def Z(self) -> int:
        return int(self.header["Z"])

    @property
    def A(self) -> float:
        return float(self.header["A"])

    @property
    def includes_radiation(self) -> bool:
        return bool(self.header.get("includes_radiation", False))

    @property
    def element(self) -> str:
        return str(self.header.get("element", ELEMENTS.get(self.Z, f"Z{self.Z}")))

    def xi(self) -> np.ndarray:
        """Natural-log density grid (1-D)."""
        return C.ln10 * self.logrho

    def tau(self) -> np.ndarray:
        """Natural-log temperature grid (1-D)."""
        return C.ln10 * self.logt


def _parse_bool(tok: str) -> bool:
    return tok.strip().lower() in ("t", "true", ".true.", "1", "yes")


def read_source_ascii(path: str | Path) -> SourceTable:
    """Parse the ASCII source format of FORMAT.md Section 1."""
    path = Path(path)
    header: dict = {}
    columns_names: list[str] = []
    data_rows: list[list[float]] = []
    with path.open() as fh:
        for raw in fh:
            line = raw.split("#", 1)[0].rstrip("\n")
            if not line.strip():
                continue
            toks = line.split()
            key = toks[0].lower()
            if key == "columns":
                columns_names = toks[1:]
                break
            if key in ("nlogrho", "nlogt"):
                # e.g. "nlogRho 321 logRho_min -8.0 logRho_max 8.0"
                header[key] = int(toks[1])
                rest = toks[2:]
                for j in range(0, len(rest) - 1, 2):
                    header[rest[j].lower()] = float(rest[j + 1])
            elif key in ("z",):
                header["Z"] = int(toks[1])
            elif key == "a":
                header["A"] = float(toks[1])
            elif key == "tolerance":
                try:
                    header["tolerance"] = float(toks[1])
                except ValueError:  # e.g. "analytic"
                    header["tolerance"] = line.split(None, 1)[1].strip()
            else:
                # value is the remainder after the key (robust to leading
                # whitespace and to long / quoted free-text header strings)
                value = line.split(None, 1)[1].strip() if len(toks) > 1 else ""
                try:
                    parts = shlex.split(value) if value else []
                except ValueError:
                    parts = [value.strip("\"'")]  # unbalanced quotes etc.
                if key == "includes_radiation":
                    header[key] = _parse_bool(parts[0]) if parts else False
                elif len(parts) == 1:
                    header[key] = parts[0]
                else:
                    header[key] = value.strip("\"'")
        # data rows
        for raw in fh:
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            data_rows.append([float(x) for x in line.split()])
    if not columns_names:
        raise ValueError(f"{path}: no 'columns' line found")
    header.setdefault("isotope", header.get("isotope", path.stem))

    data = np.array(data_rows, dtype=float)
    ncol = len(columns_names)
    if data.shape[1] != ncol:
        raise ValueError(f"{path}: {data.shape[1]} data columns, {ncol} column names")
    nlogrho = int(header["nlogrho"])
    nlogt = int(header["nlogt"])
    if data.shape[0] != nlogrho * nlogt:
        raise ValueError(
            f"{path}: {data.shape[0]} rows, expected nlogRho*nlogT={nlogrho * nlogt}"
        )
    by_name = {name: data[:, j] for j, name in enumerate(columns_names)}
    # logT outer loop, logRho inner loop -> reshape (nlogT, nlogRho)
    logrho = by_name["logRho"].reshape(nlogt, nlogrho)[0, :].copy()
    logt = by_name["logT"].reshape(nlogt, nlogrho)[:, 0].copy()
    columns = {
        name: by_name[name].reshape(nlogt, nlogrho).copy()
        for name in columns_names
        if name not in ("logRho", "logT")
    }
    return SourceTable(header=header, logrho=logrho, logt=logt, columns=columns)


def read_source_hdf5(path: str | Path) -> SourceTable:
    """Read the HDF5 source form (same keys as attributes, column datasets
    shaped ``(nlogT, nlogRho)`` plus 1-D ``logRho``, ``logT``)."""
    path = Path(path)
    header: dict = {}
    columns: dict[str, np.ndarray] = {}
    with h5py.File(path, "r") as f:
        for k, v in f.attrs.items():
            header[k.lower() if k not in ("Z", "A") else k] = _attr_to_py(v)
        logrho = np.asarray(f["logRho"][...], dtype=float)
        logt = np.asarray(f["logT"][...], dtype=float)
        for name, dset in f.items():
            if name in ("logRho", "logT"):
                continue
            columns[name] = np.asarray(dset[...], dtype=float)
    header.setdefault("isotope", path.stem)
    if "includes_radiation" in header:
        header["includes_radiation"] = bool(header["includes_radiation"])
    return SourceTable(header=header, logrho=logrho, logt=logt, columns=columns)


def read_source(path: str | Path) -> SourceTable:
    """Dispatch on extension: ``.h5``/``.hdf5`` -> HDF5, else ASCII."""
    path = Path(path)
    if path.suffix.lower() in (".h5", ".hdf5"):
        return read_source_hdf5(path)
    return read_source_ascii(path)


def write_source_ascii(
    table: SourceTable,
    path: str | Path,
    *,
    column_order: list[str] | None = None,
) -> None:
    """Write a SourceTable in the ASCII source format."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if column_order is None:
        preferred = ["P", "E", "F", "S", "Zbar", "dPdRho_T", "dPdT_Rho", "Cv"]
        column_order = [c for c in preferred if c in table.columns]
        column_order += [c for c in table.columns if c not in column_order]
    nlogt, nlogrho = table.nlogt, table.nlogrho
    hdr = table.header
    lines = ["# WDM_AV source table, version 1"]
    lines.append(f"isotope              {hdr.get('isotope', table.name)}")
    lines.append(f"Z                    {table.Z}")
    lines.append(f"A                    {float(table.A)!r}")
    lines.append("grid                 log10 uniform")
    lines.append(
        f"nlogRho              {nlogrho}     logRho_min  "
        f"{float(table.logrho[0])!r}    logRho_max   {float(table.logrho[-1])!r}"
    )
    lines.append(
        f"nlogT                {nlogt}     logT_min     "
        f"{float(table.logt[0])!r}    logT_max    {float(table.logt[-1])!r}"
    )
    rad = "T" if table.includes_radiation else "F"
    lines.append(f"includes_radiation   {rad}")
    for key in ("energy_zero", "entropy_zero", "constants", "model"):
        if key in hdr:
            lines.append(f'{key:<20} "{hdr[key]}"')
    if "tolerance" in hdr:
        lines.append(f"tolerance            {hdr['tolerance']!r}")
    lines.append("columns              logRho logT " + " ".join(column_order))

    # rows: logT outer, logRho inner
    xx, yy = np.meshgrid(table.logrho, table.logt)  # (nlogT, nlogRho)
    flat = {"logRho": xx.ravel(), "logT": yy.ravel()}
    for c in column_order:
        flat[c] = table.columns[c].ravel()
    names = ["logRho", "logT", *column_order]
    body = np.column_stack([flat[n] for n in names])
    with path.open("w") as fh:
        fh.write("\n".join(lines) + "\n")
        np.savetxt(fh, body, fmt="%.16e")


# ---------------------------------------------------------------------------
# processed tables
# ---------------------------------------------------------------------------
def _attr_to_py(v):
    if isinstance(v, bytes):
        return v.decode()
    if isinstance(v, np.generic):
        return v.item()
    return v


def write_processed_hdf5(
    path: str | Path,
    *,
    attrs: dict,
    datasets: dict[str, np.ndarray],
    audit: dict[str, np.ndarray] | None = None,
) -> None:
    """Write a processed HDF5 table per FORMAT.md Section 2, validating first.

    ``attrs`` must contain every entry of INT/DOUBLE/STR_ATTRS; ``datasets``
    must contain every REQUIRED_DATASETS entry with shape ``(ntau, nxi)``.
    Optional datasets (Zbar*, eta*, valid) are written if present.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ntau = int(attrs["ntau"])
    nxi = int(attrs["nxi"])
    shape = (ntau, nxi)

    _validate_before_write(attrs, datasets, shape)

    with h5py.File(path, "w") as f:
        for k in INT_ATTRS:
            f.attrs.create(k, np.int32(attrs[k]))
        for k in DOUBLE_ATTRS:
            f.attrs.create(k, np.float64(attrs[k]))
        for k in STR_ATTRS:
            f.attrs[k] = str(attrs[k])
        for name in REQUIRED_DATASETS:
            f.create_dataset(name, data=np.ascontiguousarray(datasets[name], float))
        for name in ("Zbar", "Zbar_xi", "Zbar_tau", "Zbar_xitau"):
            if name in datasets:
                f.create_dataset(name, data=np.ascontiguousarray(datasets[name], float))
        for name in ("eta", "eta_xi", "eta_tau", "eta_xitau"):
            if name in datasets:
                f.create_dataset(name, data=np.ascontiguousarray(datasets[name], float))
        if "valid" in datasets:
            f.create_dataset(
                "valid", data=np.ascontiguousarray(datasets["valid"], np.int8)
            )
        if audit:
            grp = f.create_group("audit")
            for name, arr in audit.items():
                grp.create_dataset(name, data=np.ascontiguousarray(arr, float))


def _validate_before_write(attrs: dict, datasets: dict, shape: tuple) -> None:
    missing = [k for k in (*INT_ATTRS, *DOUBLE_ATTRS, *STR_ATTRS) if k not in attrs]
    if missing:
        raise ValueError(f"processed table missing attributes: {missing}")
    if int(attrs["format_version"]) != FORMAT_VERSION:
        raise ValueError(f"format_version must be {FORMAT_VERSION}")
    if int(attrs["includes_radiation"]) != 0:
        raise ValueError("processed table must have includes_radiation == 0")
    if int(attrs["interp_order"]) not in (3, 5):
        raise ValueError("interp_order must be 3 or 5")
    for name in REQUIRED_DATASETS:
        if name not in datasets:
            raise ValueError(f"missing required dataset {name!r}")
        if datasets[name].shape != shape:
            raise ValueError(
                f"dataset {name!r} shape {datasets[name].shape} != {shape}"
            )
    valid = datasets.get("valid")
    if valid is not None:
        mask = np.asarray(valid, bool)
    else:
        mask = np.ones(shape, bool)
    for name in ("f", "f_xi", "f_tau"):
        if not np.all(np.isfinite(datasets[name][mask])):
            raise ValueError(f"dataset {name!r} has non-finite valid entries")


def read_processed_hdf5(path: str | Path) -> dict:
    """Read a processed table into ``{'attrs': {...}, 'datasets': {...},
    'audit': {...}}``."""
    path = Path(path)
    out: dict = {"attrs": {}, "datasets": {}, "audit": {}}
    with h5py.File(path, "r") as f:
        for k, v in f.attrs.items():
            out["attrs"][k] = _attr_to_py(v)
        for name, item in f.items():
            if isinstance(item, h5py.Group) and name == "audit":
                for an, ad in item.items():
                    out["audit"][an] = np.asarray(ad[...])
            elif isinstance(item, h5py.Dataset):
                out["datasets"][name] = np.asarray(item[...])
    return out


def validate_processed(path: str | Path) -> list[str]:
    """Re-validate a processed file against FORMAT.md.  Returns a list of
    problem strings (empty == conforms)."""
    problems: list[str] = []
    data = read_processed_hdf5(path)
    attrs, dsets = data["attrs"], data["datasets"]

    for k in INT_ATTRS:
        if k not in attrs:
            problems.append(f"missing int attr {k}")
        elif not isinstance(attrs[k], int):
            problems.append(f"attr {k} not int (got {type(attrs[k]).__name__})")
    for k in DOUBLE_ATTRS:
        if k not in attrs:
            problems.append(f"missing double attr {k}")
        elif not isinstance(attrs[k], float):
            problems.append(f"attr {k} not float (got {type(attrs[k]).__name__})")
    for k in STR_ATTRS:
        if k not in attrs:
            problems.append(f"missing string attr {k}")
        elif not isinstance(attrs[k], str):
            problems.append(f"attr {k} not str")
    if attrs.get("format_version") != FORMAT_VERSION:
        problems.append("format_version != 1")
    if attrs.get("includes_radiation") != 0:
        problems.append("includes_radiation != 0")
    if attrs.get("interp_order") not in (3, 5):
        problems.append("interp_order not in (3, 5)")

    ntau = attrs.get("ntau")
    nxi = attrs.get("nxi")
    shape = (ntau, nxi)
    for name in REQUIRED_DATASETS:
        if name not in dsets:
            problems.append(f"missing dataset {name}")
        elif dsets[name].shape != shape:
            problems.append(f"dataset {name} shape {dsets[name].shape} != {shape}")
    if "valid" in dsets and dsets["valid"].dtype != np.int8:
        problems.append("valid dataset not int8")
    return problems


def default_provenance(source_file: str = "", source_model: str = "") -> dict:
    """Common provenance attributes."""
    return {
        "source_file": source_file,
        "source_model": source_model,
        "builder_version": __version__,
        "created": datetime.datetime.now(datetime.UTC).isoformat(timespec="seconds"),
    }
