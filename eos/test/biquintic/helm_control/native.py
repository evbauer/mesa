"""Forced-cell HELM traces in the actual linear interpolation coordinates."""

import math
import struct
from decimal import Decimal, localcontext
from pathlib import Path

import numpy as np


SLOTS = ((0, 0), (1, 0), (0, 1), (2, 0), (0, 2),
         (1, 1), (2, 1), (1, 2), (2, 2))
COEFFICIENTS = ((1, 0, 0, -10, 15, -6),
                (0, 1, 0, -6, 8, -3),
                (0, 0, .5, -1.5, 1.5, -.5))
REAL = np.longdouble


def basis(fraction, width, derivative, dtype=REAL):
    if width <= 0 or derivative not in range(4):
        raise ValueError("positive linear width and derivative 0..3 required")
    weights = np.empty((2, 3), dtype=dtype)
    for slot, coefficients in enumerate(COEFFICIENTS):
        derived = np.polynomial.polynomial.polyder(np.asarray(coefficients, dtype=dtype), derivative)
        weights[0, slot] = np.polynomial.polynomial.polyval(dtype(fraction), derived)
        weights[1, slot] = (-1)**(slot + derivative) * np.polynomial.polynomial.polyval(
            dtype(1) - dtype(fraction), derived)
        weights[:, slot] *= dtype(width)**(slot - derivative)
    return weights


def tensor_jet(corners, density_width, temperature_width, density_fraction,
               temperature_fraction, density_order, temperature_order, dtype=REAL):
    density_weights = basis(density_fraction, density_width, density_order, dtype)
    temperature_weights = basis(temperature_fraction, temperature_width, temperature_order, dtype)
    terms = []
    for slot, (density_slot, temperature_slot) in enumerate(SLOTS):
        terms.extend((np.asarray(corners[slot], dtype=dtype)
                      * density_weights[:, density_slot, None]
                      * temperature_weights[None, :, temperature_slot]).ravel())
    return np.sum(terms, dtype=dtype)


def read_cache(path, density_count=2701, temperature_count=1001):
    def record(stream, count, dtype):
        marker = stream.read(4)
        if len(marker) != 4:
            raise ValueError("missing Fortran record marker")
        size = struct.unpack('<i', marker)[0]
        if size != count * np.dtype(dtype).itemsize:
            raise ValueError("unexpected HELM record size")
        values = np.fromfile(stream, dtype=dtype, count=count)
        if len(values) != count or stream.read(4) != marker:
            raise ValueError("truncated or inconsistent HELM record")
        return values

    with Path(path).open('rb') as stream:
        if record(stream, 1, '<i4')[0] != density_count:
            raise ValueError("density dimension mismatch")
        if record(stream, 1, '<i4')[0] != temperature_count:
            raise ValueError("temperature dimension mismatch")
        fields = [record(stream, density_count * temperature_count, '<f8').reshape(
            temperature_count, density_count).T for _ in SLOTS]
    if not all(np.isfinite(field).all() for field in fields):
        raise ValueError("nonfinite HELM potential fields")
    return fields


class NativeHelm:
    def __init__(self, cache, density_axis, temperature_axis):
        self.density = np.loadtxt(density_axis, dtype=np.float64).astype(REAL)
        self.temperature = np.loadtxt(temperature_axis, dtype=np.float64).astype(REAL)
        if np.any(np.diff(self.density) <= 0) or np.any(np.diff(self.temperature) <= 0):
            raise ValueError("nonmonotone HELM axes")
        self.fields = read_cache(cache, len(self.density), len(self.temperature))

    def locate(self, din, temperature):
        indices = tuple(int(np.searchsorted(axis, value, side='right')) - 1
                        for axis, value in ((self.density, din), (self.temperature, temperature)))
        if any(index < 0 or index >= len(axis) - 1 for index, axis in
               zip(indices, (self.density, self.temperature))):
            raise ValueError("outside native HELM cell domain")
        return indices

    def corners(self, cell):
        density_index, temperature_index = cell
        return np.array([field[density_index:density_index + 2,
                               temperature_index:temperature_index + 2] for field in self.fields])

    def cv(self, din, temperature, ye, cell=None, dtype=REAL):
        din, temperature, ye = dtype(din), dtype(temperature), dtype(ye)
        cell = self.locate(din, temperature) if cell is None else cell
        widths, fractions = [], []
        for value, index, axis in zip((din, temperature), cell, (self.density, self.temperature)):
            if index < 0 or index >= len(axis) - 1:
                raise ValueError("invalid forced HELM cell")
            width = dtype(axis[index + 1]) - dtype(axis[index])
            fraction = (value - dtype(axis[index])) / width
            if fraction < -1e-11 or fraction > 1 + 1e-11:
                raise ValueError("forced-cell point is not in cell")
            widths.append(width)
            fractions.append(fraction)
        corners = self.corners(cell)

        def derivative(density_order, temperature_order):
            return tensor_jet(corners, *widths, *fractions, density_order, temperature_order, dtype)

        second = derivative(0, 2)
        third = derivative(0, 3)
        mixed = derivative(1, 2)
        return np.array([-ye * temperature * second,
                         -ye * temperature * din * mixed,
                         -ye * temperature * (second + temperature * third)], dtype=dtype)


def response_jump(left, right):
    left, right = np.asarray(left), np.asarray(right)
    result = {}
    for index, field in ((1, 'Cv_x'), (2, 'Cv_t')):
        scale = max(1, abs(left[0]), abs(right[0]), abs(left[index]), abs(right[index]))
        result[field] = {'left_minus_right': float(left[index] - right[index]),
                         'scale': float(scale), 'signed_scaled': float((left[index] - right[index]) / scale)}
    return result


def decimal_basis(fraction, width, derivative):
    result = []
    for endpoint in (0, 1):
        coordinate = fraction if endpoint == 0 else 1-fraction
        values = []
        for slot, source in enumerate(COEFFICIENTS):
            coefficients = [Decimal(str(value)) for value in source]
            for _ in range(derivative):
                coefficients = [Decimal(index)*value for index, value in enumerate(coefficients)][1:]
            value = Decimal(0)
            for coefficient in reversed(coefficients):
                value = value*coordinate + coefficient
            values.append(value * width**(slot-derivative) * (1 if endpoint == 0 else (-1)**(slot+derivative)))
        result.append(values)
    return result


def decimal_cv(table, din, temperature, ye, cell, precision=80):
    with localcontext() as context:
        context.prec = precision
        din, temperature, ye = [Decimal.from_float(float(value)) for value in (din, temperature, ye)]
        density_index, temperature_index = cell
        density_axis = [Decimal.from_float(float(value)) for value in table.density[density_index:density_index+2]]
        temperature_axis = [Decimal.from_float(float(value)) for value in table.temperature[temperature_index:temperature_index+2]]
        density_width = density_axis[1]-density_axis[0]
        temperature_width = temperature_axis[1]-temperature_axis[0]
        density_fraction = (din-density_axis[0])/density_width
        temperature_fraction = (temperature-temperature_axis[0])/temperature_width
        corners = table.corners(cell)
        jets = []
        for density_order, temperature_order in ((0, 2), (0, 3), (1, 2)):
            density_weights = decimal_basis(density_fraction, density_width, density_order)
            temperature_weights = decimal_basis(temperature_fraction, temperature_width, temperature_order)
            value = Decimal(0)
            for slot, (density_slot, temperature_slot) in enumerate(SLOTS):
                for density_endpoint in (0, 1):
                    for temperature_endpoint in (0, 1):
                        value += (Decimal.from_float(float(corners[slot, density_endpoint, temperature_endpoint]))
                                  * density_weights[density_endpoint][density_slot]
                                  * temperature_weights[temperature_endpoint][temperature_slot])
            jets.append(value)
        second, third, mixed = jets
        return [-ye*temperature*second, -ye*temperature*din*mixed,
                -ye*temperature*(second+temperature*third)]


def thermal_jacobian(left, right, delta_t):
    jacobians = [-((1 + delta_t) * side[0] + delta_t * side[2]) for side in (left, right)]
    scale = max(1, abs(left[0]), abs(right[0]), *(abs(value) for value in jacobians))
    return {'delta_lnT': delta_t, 'left_over_T': float(jacobians[0]),
            'right_over_T': float(jacobians[1]), 'scale_over_T': float(scale),
            'signed_scaled': float((jacobians[0] - jacobians[1]) / scale)}


def polynomial_corners(density_axis, temperature_axis, density_power, temperature_power):
    values = np.zeros((9, 2, 2), dtype=REAL)
    for slot, (density_order, temperature_order) in enumerate(SLOTS):
        if density_order > density_power or temperature_order > temperature_power:
            continue
        prefactor = (math.factorial(density_power) / math.factorial(density_power - density_order)
                     * math.factorial(temperature_power) / math.factorial(temperature_power - temperature_order))
        values[slot] = prefactor * np.asarray(density_axis, dtype=REAL)[:, None]**(density_power - density_order) * np.asarray(
            temperature_axis, dtype=REAL)[None, :]**(temperature_power - temperature_order)
    return values
