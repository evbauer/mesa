"""Bounded scalar backward-Euler thermal solve at fixed density."""

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq


@dataclass
class Trial:
    status: str
    root: float | None = None
    iterations: int = 0
    rejected: int = 0
    bisected: int = 0
    history: list = field(default_factory=list)
    rejected_history: list = field(default_factory=list)


def residual(model, temperature_log, initial_log, heating_step):
    cv, cv_t = model(temperature_log)
    if not np.isfinite(cv) or cv <= 0 or not np.isfinite(cv_t):
        raise ValueError("nonpositive or nonfinite Cv/Cv_t")
    increment = heating_step / (np.exp(temperature_log) * cv)
    return temperature_log - initial_log - increment, 1 + increment * (1 + cv_t / cv)


def reference(model, initial_log, heating_step, lower, upper, samples):
    """Brent root with a sampled root-count/monotonicity diagnostic."""
    points = np.unique(np.r_[lower, samples, upper])
    points = points[(points >= lower) & (points <= upper)]
    if points[0] != lower or points[-1] != upper:
        raise ValueError("Reference sample does not span bracket")
    values = []
    slopes = []
    for point in points:
        value, slope = residual(model, float(point), initial_log, heating_step)
        values.append(value)
        slopes.append(slope)
    signs = np.sign(values)
    crossings = int(np.count_nonzero(signs[:-1] * signs[1:] < 0) + np.count_nonzero(signs == 0))
    diagnostic = {"sample_count": len(points), "sampled_sign_changes_or_zeros": crossings,
                  "sampled_min_jacobian": float(min(slopes)),
                  "sampled_positive_jacobian": bool(np.all(np.asarray(slopes) > 0)),
                  "uniqueness": "sampled_monotone_only" if np.all(np.asarray(slopes) > 0) else "unresolved"}
    if values[0] == 0:
        return float(lower), diagnostic
    if values[-1] == 0:
        return float(upper), diagnostic
    if values[0] * values[-1] > 0:
        diagnostic["uniqueness"] = "bracket_failure"
        return None, diagnostic
    root = brentq(lambda point: residual(model, point, initial_log, heating_step)[0],
                  lower, upper, xtol=2e-14, rtol=9e-16, maxiter=150)
    return float(root), diagnostic


def solve(model, initial_log, heating_step, lower, upper, start, tolerance=1e-11, max_iterations=40):
    """Safeguarded Newton; every trial is evaluated inside the valid bracket."""
    result = Trial("max_iterations")
    if not (lower <= start <= upper):
        result.status = "start_outside_bracket"
        return result
    try:
        lower_value = residual(model, lower, initial_log, heating_step)[0]
        upper_value = residual(model, upper, initial_log, heating_step)[0]
        if lower_value * upper_value > 0:
            result.status = "bracket_failure"
            return result
        current = float(start)
        for iteration in range(max_iterations + 1):
            value, slope = residual(model, current, initial_log, heating_step)
            result.history.append({"y": current, "residual": value, "jacobian": slope,
                                   "lower": lower, "upper": upper})
            if abs(value) <= tolerance:
                result.status, result.root, result.iterations = "converged", current, iteration
                return result
            if iteration == max_iterations:
                break
            if lower_value * value <= 0:
                upper, upper_value = current, value
            else:
                lower, lower_value = current, value
            proposal = current - value / slope if np.isfinite(slope) and slope != 0 else np.nan
            if not np.isfinite(proposal) or not lower < proposal < upper:
                if np.isfinite(proposal):
                    result.rejected_history.append({"y": float(proposal), "reason": "outside_bracket"})
                proposal = (lower + upper) / 2
                result.bisected += 1
                result.rejected += 1
            else:
                proposed_value = residual(model, proposal, initial_log, heating_step)[0]
                while abs(proposed_value) > abs(value) and abs(proposal - current) > 1e-13:
                    result.rejected_history.append({"y": float(proposal), "reason": "no_residual_decrease"})
                    proposal = (proposal + current) / 2
                    proposed_value = residual(model, proposal, initial_log, heating_step)[0]
                    result.rejected += 1
                if abs(proposed_value) > abs(value):
                    proposal = (lower + upper) / 2
                    result.bisected += 1
            current = float(proposal)
        result.iterations = max_iterations
    except (ValueError, OverflowError, FloatingPointError) as error:
        result.status = "invalid_model: " + str(error)
    return result
