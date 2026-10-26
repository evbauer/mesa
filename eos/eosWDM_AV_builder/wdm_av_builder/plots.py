"""Matplotlib maps of the audit quantities over (logRho, logT), plus a
held-out-node interpolation-error helper.  Kept deliberately simple.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .audit import run_audit
from .io import SourceTable


def _pcolor(ax, logrho, logt, z, title, *, symmetric=False, log=False):
    z = np.asarray(z, dtype=float)
    if log:
        data = np.log10(np.abs(z) + 1e-300)
        im = ax.pcolormesh(logrho, logt, data, shading="auto", cmap="viridis")
    elif symmetric:
        vmax = np.nanmax(np.abs(z)) or 1.0
        im = ax.pcolormesh(
            logrho, logt, z, shading="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax
        )
    else:
        im = ax.pcolormesh(logrho, logt, z, shading="auto", cmap="viridis")
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("logRho")
    ax.set_ylabel("logT")
    plt.colorbar(im, ax=ax)


def plot_audit_maps(source: SourceTable, path: str | Path, nodal=None) -> Path:
    """Write a grid of audit-map PNGs for one source table."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    maps, _ = run_audit(source, nodal)
    logrho, logt = source.logrho, source.logt
    panels = [
        ("precision", True, False),
        ("euler_maclaurin_xi", True, False),
        ("euler_maclaurin_tau", True, False),
        ("maxwell_residual", True, False),
        ("chi_rho_jump", True, False),
        ("cv_jump", True, False),
        ("P_rho_fd", False, True),
        ("cv_fd", False, True),
    ]
    ncol = 2
    nrow = (len(panels) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(11, 3.2 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for ax, (key, symmetric, logscale) in zip(axes, panels):
        label = key if not (symmetric or logscale) else f"{key} (log10|.|)"
        _pcolor(ax, logrho, logt, maps[key], label, symmetric=False, log=True)
    for ax in axes[len(panels) :]:
        ax.axis("off")
    fig.suptitle(f"WDM_AV audit maps: {source.name}")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path


def plot_held_out_error(
    logrho: np.ndarray,
    logt: np.ndarray,
    error: np.ndarray,
    path: str | Path,
    *,
    title: str = "held-out-node interpolation error",
) -> Path:
    """Plot a precomputed held-out-node error map (log10 of |relative error|)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 4.5))
    _pcolor(ax, logrho, logt, error, title, log=True)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path


def held_out_node_error(
    source: SourceTable, key: str = "f_xitau"
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build nodal data on the full grid and on the half grid (every other
    node), then compare the coarse-grid estimate of ``key`` at the shared nodes
    to the full-grid estimate.  Returns (logRho_half, logT_half, rel_error)."""
    from .derivs import build_nodal_data

    full = build_nodal_data(source)
    half = SourceTable(
        header=dict(source.header),
        logrho=source.logrho[::2].copy(),
        logt=source.logt[::2].copy(),
        columns={k: v[::2, ::2].copy() for k, v in source.columns.items()},
    )
    coarse = build_nodal_data(half)
    fine = full.datasets[key][::2, ::2]
    rel = (coarse.datasets[key] - fine) / (np.abs(fine) + 1e-300)
    return half.logrho, half.logt, rel


__all__ = [
    "held_out_node_error",
    "plot_audit_maps",
    "plot_held_out_error",
]
