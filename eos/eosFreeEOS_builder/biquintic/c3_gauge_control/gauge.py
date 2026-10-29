"""Known a+bT gauge transformations for a frozen coarse C3 selector audit."""

import os
import sys
from pathlib import Path

for variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[variable] = '16'

import numpy as np

HERE = Path(__file__).resolve().parent
BIQUINTIC = HERE.parent
sys.path.insert(0, str(BIQUINTIC / 'c3_spline_control'))
from regularized import solve_regularized
from spline_control import evaluate_fit
from freeeos_potential import EVAL_ORDERS, LOW_ORDERS, LN10, thermodynamic_partials


def gauge_native(sample_t, energy_shift, entropy_shift_parameter, orders=EVAL_ORDERS):
    sample_t = np.asarray(sample_t)
    jets = np.full((4, 4, *sample_t.shape), np.nan)
    temperature_term = entropy_shift_parameter * 10.0**sample_t
    for first, second in orders:
        jets[first, second] = temperature_term if first == 0 else np.zeros_like(sample_t)
        if (first, second) == (0, 0):
            jets[first, second] += energy_shift
    return jets


def shifted_low(low, logt, energy_shift, entropy_shift_parameter):
    result = low.copy()
    gauge = gauge_native(logt, energy_shift, entropy_shift_parameter, LOW_ORDERS)
    for orders in LOW_ORDERS:
        result[orders] += gauge[orders][None, :]
    return result


def remove_input_gauge(low, logt, energy_shift, entropy_shift_parameter):
    return shifted_low(low, logt, -energy_shift, -entropy_shift_parameter)


def recovered_states(coefficients, logq, logt, sample_q, sample_t,
                     energy_shift=0., entropy_shift_parameter=0., remainder_fit=False):
    fitted = evaluate_fit(coefficients, logq * LN10, logt * LN10, sample_q * LN10, sample_t * LN10)
    gauge = gauge_native(sample_t, energy_shift, entropy_shift_parameter)
    shifted = fitted.copy()
    restored_reference = fitted.copy()
    for orders in EVAL_ORDERS:
        if remainder_fit:
            shifted[orders] += gauge[orders]
        else:
            restored_reference[orders] -= gauge[orders]
    raw = thermodynamic_partials(shifted, sample_q, sample_t)
    recovered = thermodynamic_partials(restored_reference, sample_q, sample_t)
    return recovered, raw, restored_reference


def analytic_gauge_roughness(logq, logt, entropy_shift_parameter, value_scale):
    grid_tau = logt * LN10
    spacing_tau = np.median(np.diff(grid_tau))
    tau_length = grid_tau[-1] - grid_tau[0]
    average_temperature_squared = (np.exp(2 * grid_tau[-1]) - np.exp(2 * grid_tau[0])) / (2 * tau_length)
    return float(spacing_tau**8 * entropy_shift_parameter**2 * average_temperature_squared / value_scale**2)
