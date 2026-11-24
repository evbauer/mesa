"""Read-only candidate, independent-column legacy, and ideal-gas adapters."""

import sys
from pathlib import Path

import numpy as np


CODE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE))

from freeeos_potential import CRAD, LN10, evaluate, thermodynamic_partials
from thermal_consumer.models import TableModels, legacy


RESPONSE_FIELDS = ("Cv", "Cp", "chiRho", "chiT")
LEGACY_GRAD_AD_COLUMN = 6
GRAVITY_FIELDS = ("T", "Cp", "grad_ad", "chiRho", "chiT",
                  "Cp_x", "Cp_t", "grad_ad_x", "grad_ad_t",
                  "chiRho_x", "chiRho_t", "chiT_x", "chiT_t")


def physical_spline_partials(native_q, native_t):
    return native_q / LN10, (native_t - 2 * native_q) / LN10


def legacy_state(cache, logq, logt, density_log, temperature_log):
    sample_t = np.array([temperature_log / LN10])
    sample_q = np.array([density_log / LN10 - 2 * sample_t[0] + 12])
    values = legacy.evaluate_spline(cache, logq, logt, sample_q, sample_t)[0]
    native_q = legacy.evaluate_spline(cache, logq, logt, sample_q, sample_t, (1, 0))[0]
    native_t = legacy.evaluate_spline(cache, logq, logt, sample_q, sample_t, (0, 1))[0]
    density, temperature = np.exp(density_log), np.exp(temperature_log)
    gas = np.exp(values[legacy.FIELDS["lnPgas"]])
    radiation = CRAD * temperature**4 / 3
    gas_x, gas_t = physical_spline_partials(native_q[legacy.FIELDS["lnPgas"]],
                                             native_t[legacy.FIELDS["lnPgas"]])
    state = {"T": temperature, "rho": density, "Pgas": gas,
             "P": gas + radiation, "P_x": gas * gas_x,
             "P_t": gas * gas_t + 4 * radiation}
    for field in RESPONSE_FIELDS:
        column = legacy.FIELDS[field]
        state[field] = float(values[column])
        state[field + "_x"], state[field + "_t"] = map(float,
            physical_spline_partials(native_q[column], native_t[column]))
    state["grad_ad"] = float(values[LEGACY_GRAD_AD_COLUMN])
    state["grad_ad_x"], state["grad_ad_t"] = map(float,
        physical_spline_partials(native_q[LEGACY_GRAD_AD_COLUMN],
                                 native_t[LEGACY_GRAD_AD_COLUMN]))
    return state


class CoupledModels:
    def __init__(self):
        self.tables = TableModels()

    def valid(self, density_log, temperature_log):
        return self.tables.valid(density_log, temperature_log)

    def at(self, name, density_log, temperature_log, cell=None):
        if name == "ideal":
            return ideal_state(density_log, temperature_log)
        if not self.valid(density_log, temperature_log) and cell is None:
            raise ValueError("outside split-complete native support")
        if name == "legacy":
            if cell is not None:
                raise ValueError("Legacy forced-cell evaluation is unavailable")
            state = legacy_state(self.tables.cache, self.tables.logq, self.tables.logt,
                                 density_log, temperature_log)
        elif name == "candidate":
            sample_t = np.array([temperature_log / LN10])
            sample_q = np.array([density_log / LN10 - 2 * sample_t[0] + 12])
            if cell is not None:
                for coordinate, grid, index in ((sample_q, self.tables.logq, cell[0]),
                                                (sample_t, self.tables.logt, cell[1])):
                    for boundary in (grid[index], grid[index + 1]):
                        if abs(coordinate[0] - boundary) < 1e-12:
                            coordinate[0] = boundary
            cell_indices = None if cell is None else tuple(np.array([index]) for index in cell)
            native = evaluate(self.tables.nodes, self.tables.complete,
                              self.tables.logq, self.tables.logt, sample_q, sample_t,
                              cell_indices=cell_indices, add_radiation=True)
            raw = thermodynamic_partials(native, sample_q, sample_t)
            fields = ("T", "rho", "P", "Pgas", "P_x", "P_t", "Cv",
                      "Cp", "chiRho", "chiT", "grad_ad", "Cp_x", "Cp_t",
                      "chiRho_x", "chiRho_t", "chiT_x", "chiT_t",
                      "grad_ad_x", "grad_ad_t")
            state = {field: float(raw[field][0]) for field in fields}
        else:
            raise ValueError("Unknown model " + name)
        if (not all(np.isfinite(value) for value in state.values())
                or min(state["P"], state["Cv"], state["Cp"], state["chiRho"]) <= 0):
            raise ValueError("nonfinite or nonpositive pressure/response")
        return state


def ideal_state(density_log, temperature_log, gas_constant=1e8, cv=1.5e8):
    density, temperature = np.exp(density_log), np.exp(temperature_log)
    pressure = density * temperature * gas_constant
    cp = cv + gas_constant
    state = {"T": temperature, "rho": density, "P": pressure,
             "Pgas": pressure, "P_x": pressure, "P_t": pressure,
             "Cv": cv, "Cp": cp, "chiRho": 1.0, "chiT": 1.0,
             "grad_ad": gas_constant / cp}
    for field in ("Cp", "chiRho", "chiT", "grad_ad"):
        state[field + "_x"] = 0.0
        state[field + "_t"] = 0.0
    return state
