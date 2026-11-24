"""Two-variable damped Newton and independent nested-bracket reference."""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from gravity_control import gravity_residual


@dataclass
class Problem:
    model: object
    x0: float
    t0: float
    pressure_target: float
    heating: float
    energy_scale: float
    x_bounds: tuple
    t_bounds: tuple


def residual(problem, density_log, temperature_log, cell=None):
    state = problem.model(density_log, temperature_log, cell)
    gravity = gravity_residual(state, density_log - problem.x0, temperature_log - problem.t0)
    value = np.array([np.log(state["P"] / problem.pressure_target),
                      (gravity["value"] + problem.heating) / problem.energy_scale], dtype=float)
    jacobian = np.array([[state["P_x"] / state["P"], state["P_t"] / state["P"]],
                         [gravity["_x"] / problem.energy_scale,
                          gravity["_t"] / problem.energy_scale]], dtype=float)
    if not np.all(np.isfinite(value)) or not np.all(np.isfinite(jacobian)):
        raise ValueError("nonfinite coupled residual or Jacobian")
    return value, jacobian, {"P": state["P"], "epsgrav": float(gravity["value"])}


def reference(problem, sample_count=11):
    """Bracket pressure in x at fixed t, then bracket the reduced heat residual."""
    pressure_cache = {}
    pressure_slopes = []
    reduced_slopes = []

    def pressure_root(temperature_log):
        if temperature_log in pressure_cache:
            return pressure_cache[temperature_log]

        def pressure_value(density_log):
            values, jacobian, _ = residual(problem, density_log, temperature_log)
            return float(values[0]), float(jacobian[0, 0])

        lower_value, lower_slope = pressure_value(problem.x_bounds[0])
        upper_value, upper_slope = pressure_value(problem.x_bounds[1])
        pressure_slopes.extend((lower_slope, upper_slope))
        if lower_value * upper_value > 0:
            raise ValueError("pressure_not_bracketed")
        density_log = brentq(lambda point: pressure_value(point)[0], *problem.x_bounds,
                             xtol=2e-14, rtol=9e-16, maxiter=100)
        pressure_cache[temperature_log] = float(density_log)
        return float(density_log)

    def reduced(temperature_log):
        density_log = pressure_root(float(temperature_log))
        values, jacobian, _ = residual(problem, density_log, temperature_log)
        if abs(jacobian[0, 0]) < 1e-14:
            raise ValueError("singular_pressure_slope")
        reduced_slopes.append(float(jacobian[1, 1] - jacobian[1, 0] * jacobian[0, 1] / jacobian[0, 0]))
        return float(values[1])

    try:
        points = np.linspace(*problem.t_bounds, sample_count)
        values = np.array([reduced(float(point)) for point in points])
        sign_changes = int(np.count_nonzero(values[:-1] * values[1:] < 0)
                           + np.count_nonzero(values == 0))
        diagnostic = {"temperature_samples": sample_count,
                      "sampled_sign_changes_or_zeros": sign_changes,
                      "sampled_pressure_slope_min": float(min(pressure_slopes)),
                      "sampled_reduced_slope_min": float(min(reduced_slopes)),
                      "sampled_reduced_slope_max": float(max(reduced_slopes)),
                      "sampled_monotone": bool(np.all(np.asarray(reduced_slopes) > 0)
                                               or np.all(np.asarray(reduced_slopes) < 0)),
                      "uniqueness_limit": "sampled only; no global uniqueness proof"}
        if values[0] * values[-1] > 0:
            return {"status": "heat_not_bracketed", "diagnostic": diagnostic}
        temperature_log = float(brentq(reduced, *problem.t_bounds,
                                       xtol=2e-14, rtol=9e-16, maxiter=100))
        density_log = pressure_root(temperature_log)
        return {"status": "bracketed", "root": [density_log, temperature_log],
                "residual_inf": float(np.max(np.abs(residual(problem, density_log, temperature_log)[0]))),
                "diagnostic": diagnostic}
    except (ValueError, OverflowError, FloatingPointError) as error:
        return {"status": "reference_failure", "reason": str(error),
                "diagnostic": {"sampled_pressure_slope_count": len(pressure_slopes),
                               "sampled_reduced_slope_count": len(reduced_slopes)}}


def solve(problem, start, tolerance=1e-10, max_iterations=25, minimum_damping=2**-12,
          condition_limit=1e12):
    current = np.asarray(start, dtype=float)
    history, rejected = [], []
    box = np.asarray((problem.x_bounds, problem.t_bounds), dtype=float)
    if np.any(current < box[:, 0]) or np.any(current > box[:, 1]):
        return {"status": "start_outside_box", "history": history, "rejected_history": rejected}
    for iteration in range(max_iterations + 1):
        try:
            value, jacobian, _ = residual(problem, *current)
        except (ValueError, OverflowError, FloatingPointError) as error:
            return {"status": "invalid_iterate", "reason": str(error),
                    "iterations": iteration, "history": history, "rejected_history": rejected}
        norm = float(np.max(np.abs(value)))
        condition = float(np.linalg.cond(jacobian))
        history.append({"x": float(current[0]), "t": float(current[1]),
                        "residual": value.tolist(), "residual_inf": norm,
                        "jacobian": jacobian.tolist(), "condition": condition})
        if norm <= tolerance:
            return {"status": "converged", "root": current.tolist(),
                    "iterations": iteration, "history": history,
                    "rejected_history": rejected, "rejected": len(rejected),
                    "damped_steps": sum(item.get("accepted_damping", 1.0) < 1.0 for item in history)}
        if iteration == max_iterations:
            break
        if not np.isfinite(condition) or condition > condition_limit:
            return {"status": "singular_or_ill_conditioned", "iterations": iteration,
                    "history": history, "rejected_history": rejected}
        try:
            direction = np.linalg.solve(jacobian, -value)
        except np.linalg.LinAlgError:
            return {"status": "singular_jacobian", "iterations": iteration,
                    "history": history, "rejected_history": rejected}
        damping = 1.0
        while damping >= minimum_damping:
            proposal = current + damping * direction
            if np.any(proposal <= box[:, 0]) or np.any(proposal >= box[:, 1]):
                reason = "outside_box"
            else:
                try:
                    proposed, _, _ = residual(problem, *proposal)
                    proposed_norm = float(np.max(np.abs(proposed)))
                    reason = "no_decrease" if proposed_norm > max(tolerance, norm * (1 - 1e-4 * damping)) else None
                except (ValueError, OverflowError, FloatingPointError) as error:
                    reason = "invalid_trial: " + str(error)
            if reason is None:
                history[-1]["accepted_damping"] = damping
                current = proposal
                break
            rejected.append({"x": float(proposal[0]), "t": float(proposal[1]),
                             "damping": damping, "reason": reason})
            damping *= .5
        else:
            return {"status": "line_search_failure", "iterations": iteration,
                    "history": history, "rejected_history": rejected, "rejected": len(rejected)}
    return {"status": "max_iterations", "iterations": max_iterations,
            "history": history, "rejected_history": rejected, "rejected": len(rejected)}
