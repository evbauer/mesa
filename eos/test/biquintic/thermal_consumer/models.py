"""Read-only adapters for the verified potential and independent-Cv cache."""

import sys
from pathlib import Path

import numpy as np


CODE = Path(__file__).resolve().parents[1]
WORKSPACE = CODE.parents[3]
sys.path.insert(0, str(CODE))
sys.path.insert(0, str(WORKSPACE / "freeeos_report"))

import audit_freeeos as legacy
from freeeos_potential import (
    LN10, checksum, construct_twists, evaluate, query_cells, read_table,
    source_constraints, thermodynamic_partials, valid_cells,
)


SOURCE = WORKSPACE / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
CACHE = WORKSPACE / "mesa/data/eosDT_data/cache/mesa-FreeEOS_02z70x.bin"
SOURCE_SHA = "e6ca7ce27bee8197bb257a8c4886a06366e6951e58958fbed3735c877ee68502"
CORE_SHA = "40b513e48b512e44ea2c020cd29c000b60128661260f38c0c0b7365c2ae10a7d"


class TableModels:
    def __init__(self):
        if checksum(SOURCE) != SOURCE_SHA or checksum(CODE / "freeeos_potential.py") != CORE_SHA:
            raise ValueError("Frozen source/kernel identity changed")
        header, self.logq, self.logt, rows = read_table(SOURCE)
        if abs(header[1] - .7) > 1e-12 or abs(header[2] - .02) > 1e-12:
            raise ValueError("Unexpected source composition")
        low, valid, _ = source_constraints(rows, self.logq, self.logt,
                                           split_radiation=True, minimum_logT=3.2)
        self.nodes, complete, _ = construct_twists(
            low, valid, (self.logq[1] - self.logq[0]) * LN10,
            (self.logt[1] - self.logt[0]) * LN10, radius=2, route="mean")
        self.cells = valid_cells(complete)
        self.complete = complete
        _, self.cache = legacy.read_cache(CACHE, len(self.logq), len(self.logt))

    def valid(self, density_log, temperature_log):
        sample_t = temperature_log / LN10
        sample_q = density_log / LN10 - 2 * sample_t + 12
        density_index, temperature_index, inside = query_cells(
            self.logq, self.logt, np.array([sample_q]), np.array([sample_t]))
        return bool(inside[0] and self.cells[density_index[0], temperature_index[0]])

    def at(self, name, density_log, temperature_log, cell=None):
        sample_t = np.array([temperature_log / LN10])
        sample_q = np.array([density_log / LN10 - 2 * sample_t[0] + 12])
        if name == "candidate":
            indices = None if cell is None else tuple(np.array([index]) for index in cell)
            native = evaluate(self.nodes, self.complete, self.logq, self.logt, sample_q,
                              sample_t, cell_indices=indices, add_radiation=True)
            state = thermodynamic_partials(native, sample_q, sample_t)
            return float(state["Cv"][0]), float(state["Cv_t"][0])
        if name == "legacy":
            if cell is not None:
                raise ValueError("Legacy evaluator does not support forced cells")
            values = legacy.evaluate_spline(self.cache, self.logq, self.logt, sample_q, sample_t)
            density_partial = legacy.evaluate_spline(self.cache, self.logq, self.logt,
                                                     sample_q, sample_t, (1, 0))
            temperature_partial = legacy.evaluate_spline(self.cache, self.logq, self.logt,
                                                         sample_q, sample_t, (0, 1))
            field = legacy.FIELDS["Cv"]
            return float(values[0, field]), float((temperature_partial[0, field]
                         - 2 * density_partial[0, field]) / LN10)
        raise ValueError("Unknown table model")


class SmoothModel:
    def __init__(self, cv):
        self.cv = float(cv)

    def __call__(self, temperature_log):
        return self.cv, 0.
