"""Native-only FreeEOS potential preprocessing; experimental, not a runtime EOS."""

import hashlib
import math
from pathlib import Path

import numpy as np


LN10 = np.log(10.0)
CRAD = 8.0 * np.pi**5 * (1.380649e-16)**4 / (
    15.0 * (6.62607015e-27)**3 * (2.99792458e10)**3
)
LOW_ORDERS = ((0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2))
EVAL_ORDERS = tuple((density_order, total - density_order)
                    for total in range(4) for density_order in range(total + 1))


def checksum(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_table(path, validate_endpoints=True):
    records = []
    with Path(path).open() as stream:
        next(stream)
        header = np.fromstring(next(stream), sep=" ")
        if header.size != 11 or header[0] != 51:
            raise ValueError("Expected FreeEOS format 51 header")
        if not np.all(np.isfinite(header)) or header[6] <= 0 or header[10] <= 0:
            raise ValueError("Nonfinite header or nonpositive grid spacing")
        if header[3] != int(header[3]) or header[7] != int(header[7]) or min(header[3], header[7]) < 2:
            raise ValueError("Grid counts must be integer and at least two")
        if header[1] < 0 or header[2] < 0 or header[1] + header[2] > 1:
            raise ValueError("Invalid source composition header")
        if validate_endpoints:
            derived = np.array([header[4] + (header[3] - 1) * header[6],
                                header[8] + (header[7] - 1) * header[10]])
            if not np.allclose(derived, header[[5, 9]], rtol=0, atol=1e-10):
                raise ValueError("Displayed grid endpoint differs from count/minimum/step")
        for line in stream:
            columns = line.split()
            if len(columns) == 22 and columns[0][0] in "+-0123456789":
                records.append(np.asarray(columns, dtype=float))
    density_count, temperature_count = int(header[7]), int(header[3])
    if len(records) != density_count * temperature_count:
        raise ValueError("Incomplete FreeEOS source table")
    rows = np.asarray(records).reshape(density_count, temperature_count, 22)
    logq = header[8] + np.arange(density_count) * header[10]
    logt = header[4] + np.arange(temperature_count) * header[6]
    if not np.allclose(rows[:, :, 0], logt[None, :], rtol=0, atol=1e-10):
        raise ValueError("Temperature grid does not match header")
    return header, logq, logt, rows


def to_native(physical, orders=LOW_ORDERS):
    native = np.full_like(physical, np.nan)
    for density_order, temperature_order in orders:
        native[density_order, temperature_order] = sum(
            math.comb(temperature_order, added_density) * 2**added_density
            * physical[density_order + added_density, temperature_order - added_density]
            for added_density in range(temperature_order + 1))
    return native


def to_physical(native, orders=LOW_ORDERS):
    physical = np.full_like(native, np.nan)
    for density_order, temperature_order in orders:
        physical[density_order, temperature_order] = sum(
            math.comb(temperature_order, added_density) * (-2)**added_density
            * native[density_order + added_density, temperature_order - added_density]
            for added_density in range(temperature_order + 1))
    return physical


def radiation_potential(sample_q, sample_t, orders=EVAL_ORDERS):
    temperature = 10.0**np.asarray(sample_t)
    density = 10.0**(np.asarray(sample_q) + 2 * np.asarray(sample_t) - 12)
    free_energy = -CRAD * temperature**4 / (3 * density)
    highest_order = max(max(derivative_orders) for derivative_orders in orders) + 1
    result = np.full((highest_order, highest_order, *free_energy.shape), np.nan)
    for density_order, temperature_order in orders:
        result[density_order, temperature_order] = (-1)**density_order * 2**temperature_order * free_energy
    return result


def source_constraints(rows, logq, logt, split_radiation=False):
    temperature = np.broadcast_to(10.0**logt[None, :], rows.shape[:2])
    density = 10.0**(logq[:, None] + 2.0 * logt[None, :] - 12.0)
    gas_pressure = 10.0**rows[:, :, 1]
    pressure = gas_pressure + CRAD * temperature**4 / 3.0
    energy, entropy = 10.0**rows[:, :, 2], 10.0**rows[:, :, 3]
    physical = np.full((3, 3, *density.shape), np.nan)
    physical[0, 0] = energy - temperature * entropy
    physical[1, 0] = pressure / density
    physical[0, 1] = -temperature * entropy
    physical[2, 0] = physical[1, 0] * (rows[:, :, 4] - 1.0)
    physical[1, 1] = physical[1, 0] * rows[:, :, 5]
    physical[0, 2] = -temperature * (entropy + rows[:, :, 7])
    if split_radiation:
        radiation_pressure = CRAD * temperature**4 / 3
        radiation_scale = radiation_pressure / density
        physical[0, 0] += radiation_scale
        physical[1, 0] = gas_pressure / density
        physical[0, 1] += 4 * radiation_scale
        physical[2, 0] = (pressure * rows[:, :, 4] - gas_pressure) / density
        physical[1, 1] = (pressure * rows[:, :, 5] - 4 * radiation_pressure) / density
        physical[0, 2] += 16 * radiation_scale
    valid = (rows[:, :, 17] == 0) & (logt[None, :] >= 3.2 - 1e-12)
    valid &= abs(rows[:, :, 18] - np.log10(density)) < 1e-10
    valid &= np.all(np.isfinite(rows[:, :, 1:11]), axis=2)
    for orders in LOW_ORDERS:
        valid &= np.isfinite(physical[orders])
    valid &= (pressure > 0) & (energy > 0) & (entropy > 0) & (rows[:, :, 7] > 0) & (rows[:, :, 4] > 0)
    native = to_native(physical)
    native[:, :, ~valid] = np.nan
    return native, valid, {"rho": density, "T": temperature, "P": pressure,
                           "Pgas": gas_pressure, "E": energy, "S": entropy,
                           "Cv": rows[:, :, 7]}


def local_derivative(values, valid, spacing, axis, radius=2):
    if radius not in (1, 2, 3) or not np.isfinite(spacing) or spacing <= 0:
        raise ValueError("Require positive spacing and centered radius 1, 2, or 3")
    if values.shape != valid.shape or axis not in (0, 1):
        raise ValueError("Derivative values and validity must be matching 2D arrays")
    offsets = np.arange(-radius, radius + 1)
    vandermonde = offsets[None, :] ** np.arange(offsets.size)[:, None]
    target = np.zeros(offsets.size)
    target[1] = 1.0
    weights = np.linalg.solve(vandermonde.astype(float), target) / spacing
    result = np.full_like(values, np.nan, dtype=float)
    accepted = np.zeros_like(valid, dtype=bool)
    interior = [slice(None), slice(None)]
    interior[axis] = slice(radius, values.shape[axis] - radius)
    interior = tuple(interior)
    if values.shape[axis] <= 2 * radius:
        return result, accepted
    total = np.zeros_like(values[interior])
    safe = np.ones_like(valid[interior])
    for offset, weight in zip(offsets, weights):
        shifted = [slice(None), slice(None)]
        shifted[axis] = slice(radius + offset, values.shape[axis] - radius + offset)
        shifted = tuple(shifted)
        safe &= valid[shifted] & np.isfinite(values[shifted])
        if abs(weight) > 1e-14 / spacing:
            total += weight * np.where(valid[shifted], values[shifted], 0.0)
    accepted[interior] = safe
    result[interior] = np.where(safe, total, np.nan)
    return result, accepted


def construct_twists(low, valid, du, dtau, radius=2, route="mean"):
    if route not in ("mean", "first", "second"):
        raise ValueError("Twist route must be mean, first, or second")
    nodes = low.copy()
    routes = {}
    first21, mask21a = local_derivative(low[2, 0], valid, dtau, 1, radius)
    second21, mask21b = local_derivative(low[1, 1], valid, du, 0, radius)
    first12, mask12a = local_derivative(low[1, 1], valid, dtau, 1, radius)
    second12, mask12b = local_derivative(low[0, 2], valid, du, 0, radius)
    twist_valid = mask21a & mask21b & mask12a & mask12b
    for orders, first, second in (((2, 1), first21, second21), ((1, 2), first12, second12)):
        routes[orders] = (first, second, twist_valid.copy())
        nodes[orders] = {"mean": (first + second) / 2, "first": first, "second": second}[route]
    first22, mask22a = local_derivative(nodes[2, 1], twist_valid, dtau, 1, radius)
    second22, mask22b = local_derivative(nodes[1, 2], twist_valid, du, 0, radius)
    complete = twist_valid & mask22a & mask22b
    routes[2, 2] = (first22, second22, complete.copy())
    nodes[2, 2] = {"mean": (first22 + second22) / 2, "first": first22, "second": second22}[route]
    nodes[:, :, ~complete] = np.nan
    return nodes, complete, routes


def valid_cells(valid):
    return valid[:-1, :-1] & valid[1:, :-1] & valid[:-1, 1:] & valid[1:, 1:]


def query_cells(logq, logt, sample_q, sample_t):
    sample_q, sample_t = np.broadcast_arrays(sample_q, sample_t)
    inside = np.isfinite(sample_q) & np.isfinite(sample_t)
    inside &= (sample_q >= logq[0]) & (sample_q <= logq[-1])
    inside &= (sample_t >= logt[0]) & (sample_t <= logt[-1])
    density_index = np.clip(np.searchsorted(logq, sample_q, side="right") - 1, 0, len(logq) - 2)
    temperature_index = np.clip(np.searchsorted(logt, sample_t, side="right") - 1, 0, len(logt) - 2)
    return density_index, temperature_index, inside


def quintic_basis(position, spacing, order):
    coefficients = np.array([
        [1, 0, 0, -10, 15, -6], [0, 1, 0, -6, 8, -3],
        [0, 0, .5, -1.5, 1.5, -.5], [0, 0, 0, 10, -15, 6],
        [0, 0, 0, -4, 7, -3], [0, 0, 0, .5, -1, .5],
    ], dtype=float)
    result = np.empty((len(position), 2, 3))
    for endpoint in range(2):
        for derivative in range(3):
            polynomial = np.polynomial.polynomial.polyder(coefficients[3 * endpoint + derivative], order)
            result[:, endpoint, derivative] = np.polynomial.polynomial.polyval(position, polynomial) * spacing**(derivative - order)
    return result


def evaluate(nodes, complete, logq, logt, sample_q, sample_t, orders=EVAL_ORDERS, cell_indices=None, add_radiation=False):
    sample_q, sample_t = np.asarray(sample_q), np.asarray(sample_t)
    density_index, temperature_index, inside = query_cells(logq, logt, sample_q, sample_t)
    if cell_indices is not None:
        density_index, temperature_index = map(np.asarray, cell_indices)
        if density_index.shape != sample_q.shape or temperature_index.shape != sample_t.shape:
            raise ValueError("Forced cell indices must match query arrays")
        if (np.any(density_index < 0) or np.any(density_index >= len(logq) - 1)
                or np.any(temperature_index < 0) or np.any(temperature_index >= len(logt) - 1)):
            raise ValueError("Forced cell indices outside grid")
        inside &= (sample_q >= logq[density_index]) & (sample_q <= logq[density_index + 1])
        inside &= (sample_t >= logt[temperature_index]) & (sample_t <= logt[temperature_index + 1])
    if not np.all(inside & valid_cells(complete)[density_index, temperature_index]):
        raise ValueError("Query outside table or complete native-only stencil domain")
    dq, dt = logq[1] - logq[0], logt[1] - logt[0]
    density_basis = {order: quintic_basis((sample_q - logq[density_index]) / dq, dq * LN10, order)
                     for order in range(4)}
    temperature_basis = {order: quintic_basis((sample_t - logt[temperature_index]) / dt, dt * LN10, order)
                         for order in range(4)}
    result = np.full((4, 4, sample_q.size), np.nan)
    for density_order, temperature_order in orders:
        total = np.zeros(sample_q.size)
        for density_end in range(2):
            for temperature_end in range(2):
                for nodal_density in range(3):
                    for nodal_temperature in range(3):
                        total += nodes[nodal_density, nodal_temperature,
                                       density_index + density_end, temperature_index + temperature_end] * (
                            density_basis[density_order][:, density_end, nodal_density]
                            * temperature_basis[temperature_order][:, temperature_end, nodal_temperature])
        result[density_order, temperature_order] = total
    if add_radiation:
        radiation = radiation_potential(sample_q, sample_t, orders)
        for derivative_orders in orders:
            result[derivative_orders] += radiation[derivative_orders]
    return result


def thermodynamics(native, sample_q, sample_t):
    physical = to_physical(native)
    density = 10.0**(sample_q + 2.0 * sample_t - 12.0)
    temperature = 10.0**sample_t
    pressure = density * physical[1, 0]
    energy = physical[0, 0] - physical[0, 1]
    entropy = -physical[0, 1] / temperature
    cv = (physical[0, 1] - physical[0, 2]) / temperature
    with np.errstate(divide="ignore", invalid="ignore"):
        chirho = 1 + physical[2, 0] / physical[1, 0]
        chit = physical[1, 1] / physical[1, 0]
        gamma3 = pressure * chit / (density * temperature * cv)
        gamma1 = chirho + chit * gamma3
        cp = cv + pressure * chit**2 / (density * temperature * chirho)
        grada = gamma3 / gamma1
    return {"P": pressure, "Pgas": pressure - CRAD * temperature**4 / 3.0,
            "E": energy, "S": entropy, "Cv": cv, "Cp": cp, "chiRho": chirho,
            "chiT": chit, "gamma1": gamma1, "gamma3": 1 + gamma3, "gamma3_minus1": gamma3, "grad_ad": grada,
            "P_rho": physical[1, 0] + physical[2, 0],
            "P_T": density * physical[1, 1] / temperature,
            "E_rho": (physical[1, 0] - physical[1, 1]) / density,
            "E_T": cv, "S_rho": -physical[1, 1] / (density * temperature),
            "S_T": cv / temperature, "rho": density, "T": temperature}


def thermodynamic_partials(native, sample_q, sample_t):
    state = thermodynamics(native, sample_q, sample_t)
    physical = to_physical(native, EVAL_ORDERS)
    temperature, density = state["T"], state["rho"]
    fx, fxx, fxt = physical[1, 0], physical[2, 0], physical[1, 1]
    fxxx, fxxt, fxtt = physical[3, 0], physical[2, 1], physical[1, 2]
    derivatives = {}
    derivatives["P"] = (density * (fx + fxx), density * fxt)
    derivatives["Pgas"] = (derivatives["P"][0], derivatives["P"][1] - 4 * CRAD * temperature**4 / 3)
    derivatives["E"] = (fx - fxt, temperature * state["Cv"])
    derivatives["S"] = (-fxt / temperature, state["Cv"])
    derivatives["Cv"] = ((fxt - fxtt) / temperature,
                         (2 * physical[0, 2] - physical[0, 1] - physical[0, 3]) / temperature)
    with np.errstate(divide="ignore", invalid="ignore"):
        derivatives["chiRho"] = ((fxxx * fx - fxx**2) / fx**2,
                                 (fxxt * fx - fxx * fxt) / fx**2)
        derivatives["chiT"] = ((fxxt * fx - fxt * fxx) / fx**2,
                               (fxtt * fx - fxt**2) / fx**2)
        excess = state["Cp"] - state["Cv"]
        cp_density = (fxx * state["chiT"]**2 + 2 * fx * state["chiT"] * derivatives["chiT"][0]) / (temperature * state["chiRho"])
        cp_density -= excess * derivatives["chiRho"][0] / state["chiRho"]
        cp_temperature = (fxt * state["chiT"]**2 + 2 * fx * state["chiT"] * derivatives["chiT"][1]) / (temperature * state["chiRho"])
        cp_temperature -= excess * (1 + derivatives["chiRho"][1] / state["chiRho"])
        derivatives["Cp"] = (derivatives["Cv"][0] + cp_density, derivatives["Cv"][1] + cp_temperature)
        derivatives["gamma3"] = (fxxt / (temperature * state["Cv"]) - state["gamma3_minus1"] * derivatives["Cv"][0] / state["Cv"],
                                  fxtt / (temperature * state["Cv"]) - state["gamma3_minus1"] * (1 + derivatives["Cv"][1] / state["Cv"]))
        derivatives["gamma1"] = tuple(derivatives["chiRho"][axis] + derivatives["chiT"][axis] * state["gamma3_minus1"]
                                       + state["chiT"] * derivatives["gamma3"][axis] for axis in range(2))
        derivatives["grad_ad"] = tuple((derivatives["gamma3"][axis] * state["gamma1"]
                                        - state["gamma3_minus1"] * derivatives["gamma1"][axis]) / state["gamma1"]**2
                                       for axis in range(2))
    for field, (density_partial, temperature_partial) in derivatives.items():
        state[field + "_x"] = density_partial
        state[field + "_t"] = temperature_partial
    return state


def identity_residuals(state, balanced=False):
    pressure, density, temperature = state["P"], state["rho"], state["T"]
    first = density**2 * state["E_rho"]
    second = temperature * state["P_T"]
    if not balanced:
        with np.errstate(divide="ignore", invalid="ignore"):
            return {"dpe": (first + second - pressure) / pressure,
                    "dse": temperature * state["S_T"] / state["E_T"] - 1,
                    "dsp": -density**2 * state["S_rho"] / state["P_T"] - 1}
    numerators = {"dpe": first + second - pressure,
                  "dse": temperature * state["S_T"] - state["E_T"],
                  "dsp": density**2 * state["S_rho"] + state["P_T"]}
    denominators = {"dpe": np.maximum.reduce([abs(first), abs(second), abs(pressure)]),
                    "dse": np.maximum(abs(temperature * state["S_T"]), abs(state["E_T"])),
                    "dsp": np.maximum(abs(density**2 * state["S_rho"]), abs(state["P_T"]))}
    return {name: np.divide(numerator, denominators[name], out=np.zeros_like(numerator), where=denominators[name] != 0)
            for name, numerator in numerators.items()}
