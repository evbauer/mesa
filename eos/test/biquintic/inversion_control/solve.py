"""Bounded log-temperature Newton/bisection controls with analytic Jacobians."""

import numpy as np


def solve(callback, target, lower, upper, starting, safeguarded, rtol=1e-10, max_iterations=60):
    target, lower, upper, starting = [np.asarray(values, dtype=float).copy()
                                     for values in (target, lower, upper, starting)]
    count = target.size
    indices = np.arange(count)
    low_value, _, _ = callback(lower, indices)
    high_value, _, _ = callback(upper, indices)
    bracketed = np.isfinite(low_value) & np.isfinite(high_value) & (low_value <= target) & (high_value >= target)
    bracketed &= (lower < upper) & (starting >= lower) & (starting <= upper) & (target > 0)
    status = np.where(bracketed, 'running', 'invalid_bracket').astype('<U24')
    current = starting.copy()
    declared_lower, declared_upper = lower.copy(), upper.copy()
    iterations = np.zeros(count, dtype=int)
    safeguards = np.zeros(count, dtype=int)
    unsafe_cv = np.zeros(count, dtype=int)
    escapes = np.zeros(count, dtype=int)
    face_crossings = np.zeros(count, dtype=int)
    residual = np.full(count, np.nan)
    minimum_jacobian_condition = np.full(count, np.inf)
    for iteration in range(1, max_iterations + 1):
        active = np.flatnonzero(status == 'running')
        if not active.size:
            break
        energy, jacobian, cv = callback(current[active], active)
        difference = energy - target[active]
        residual[active] = difference / target[active]
        iterations[active] = iteration
        finite = np.isfinite(energy) & np.isfinite(jacobian) & np.isfinite(cv)
        positive = finite & (jacobian > 0) & (cv > 0)
        unsafe_cv[active[~positive]] += 1
        minimum_jacobian_condition[active] = np.minimum(
            minimum_jacobian_condition[active], np.where(positive, jacobian / target[active], 0))
        done = finite & (abs(residual[active]) <= rtol)
        status[active[done]] = 'converged'
        bad_value = ~np.isfinite(energy)
        status[active[bad_value]] = 'nonfinite_energy'
        pending = ~(done | bad_value)
        active, difference, jacobian, positive = (values[pending]
                                                  for values in (active, difference, jacobian, positive))
        if not active.size:
            continue
        lower[active] = np.where(difference <= 0, current[active], lower[active])
        upper[active] = np.where(difference > 0, current[active], upper[active])
        with np.errstate(divide='ignore', invalid='ignore'):
            proposal = current[active] - difference / jacobian
        if safeguarded:
            accepted = positive & np.isfinite(proposal)
            accepted &= (proposal > lower[active]) & (proposal < upper[active])
            accepted &= abs(proposal - current[active]) <= .75 * (upper[active] - lower[active])
            safeguards[active[~accepted]] += 1
            proposal = np.where(accepted, proposal, .5 * (lower[active] + upper[active]))
        else:
            status[active[~positive]] = 'unsafe_Cv'
            escaped = positive & (~np.isfinite(proposal) | (proposal < declared_lower[active])
                                  | (proposal > declared_upper[active]))
            escapes[active[escaped]] += 1
            status[active[escaped]] = 'left_declared_bracket'
        accepted = status[active] == 'running'
        if hasattr(callback, 'crossings'):
            accepted_indices = active[accepted]
            face_crossings[accepted_indices] += callback.crossings(
                current[accepted_indices], proposal[accepted], accepted_indices)
        current[active[accepted]] = proposal[accepted]
    status[status == 'running'] = 'iteration_limit'
    return {'status': status, 'lnT': current, 'relative_energy_residual': residual,
            'iterations': iterations, 'safeguards': safeguards, 'unsafe_Cv': unsafe_cv,
            'out_of_bracket_attempts': escapes, 'accepted_face_crossings': face_crossings,
            'minimum_TCv_over_target_E': minimum_jacobian_condition}
