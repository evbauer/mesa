"""Bounded offline composition controls, never a MESA selection route."""

import sys
import warnings
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from freeeos_potential import EVAL_ORDERS, LN10, to_physical


SUPPORTS = {"linear": (.6, .8), "cubic_lagrange": (.5, .6, .8, .9)}
TARGET_X = .7


def composition_weights(scheme, target=TARGET_X, warn=True):
    support = np.asarray(SUPPORTS[scheme])
    if not np.isfinite(target) or not .6 <= target <= .8:
        raise ValueError("This bounded control only accepts .60 <= X <= .80")
    weights = np.ones(len(support))
    partials = np.zeros(len(support))
    for index, coordinate in enumerate(support):
        others = np.delete(support, index)
        denominator = np.prod(coordinate - others)
        weights[index] = np.prod(target - others) / denominator
        partials[index] = sum(np.prod(target - np.delete(others, omitted))
                              for omitted in range(len(others))) / denominator
    if not np.isclose(weights.sum(), 1, rtol=0, atol=2e-14):
        raise ValueError("Composition weights do not sum to one")
    if warn and np.any(weights < 0):
        warnings.warn("Cubic Lagrange has NEGATIVE composition weights: not convex; "
                      "positive sources do not guarantee a stable mixed EOS.", RuntimeWarning,
                      stacklevel=2)
    return dict(zip(support, weights)), dict(zip(support, partials))


def combine(arrays, weights):
    return sum(weight * arrays[coordinate] for coordinate, weight in weights.items())


def gauge_jet(sample_q, sample_t, constant, linear_temperature):
    temperature = 10.0**np.asarray(sample_t)
    native = np.zeros((4, 4, *temperature.shape))
    native[0, 0] = constant + linear_temperature * temperature
    for order in (1, 2, 3):
        native[0, order] = linear_temperature * temperature
    return native


def composition_partials(native, partial_native, state):
    physical = to_physical(partial_native, EVAL_ORDERS)
    baseline = to_physical(native, EVAL_ORDERS)
    density, temperature = state["rho"], state["T"]
    result = {"P": density * physical[1, 0],
              "Pgas": density * physical[1, 0],
              "E": physical[0, 0] - physical[0, 1],
              "S": -physical[0, 1] / temperature,
              "Cv": (physical[0, 1] - physical[0, 2]) / temperature}
    with np.errstate(divide="ignore", invalid="ignore"):
        result["chiRho"] = (physical[2, 0] - (state["chiRho"] - 1) * physical[1, 0]) / baseline[1, 0]
        result["chiT"] = (physical[1, 1] - state["chiT"] * physical[1, 0]) / baseline[1, 0]
        result["gamma3"] = (physical[1, 1] / (temperature * state["Cv"])
                            - state["gamma3_minus1"] * result["Cv"] / state["Cv"])
        result["gamma1"] = (result["chiRho"] + result["chiT"] * state["gamma3_minus1"]
                            + state["chiT"] * result["gamma3"])
        result["Cp"] = result["Cv"] + (
            physical[1, 0] * state["chiT"]**2 + 2 * baseline[1, 0] * state["chiT"] * result["chiT"]
            - baseline[1, 0] * state["chiT"]**2 * result["chiRho"] / state["chiRho"]
        ) / (temperature * state["chiRho"])
        result["grad_ad"] = (result["gamma3"] - state["grad_ad"] * result["gamma1"]) / state["gamma1"]
    return result


def grid_spacing(logq, logt):
    return (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10
