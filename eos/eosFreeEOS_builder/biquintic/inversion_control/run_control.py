"""Fixed-density table-model self-inversion; no source accuracy or stellar-solver claim."""

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

SIDECAR = Path(__file__).resolve().parent
BIQUINTIC = SIDECAR.parent
WORKSPACE = BIQUINTIC.parents[3]
sys.path.insert(0, str(BIQUINTIC))
sys.path.insert(0, str(BIQUINTIC / 'c3_control'))

from freeeos_potential import (
    CRAD, LN10, LOW_ORDERS, construct_twists, evaluate, query_cells, read_table,
    source_constraints, thermodynamics, valid_cells,
)
from c3_potential import construct_jets, evaluate as evaluate_c3
from solve import solve


OUTPUT = WORKSPACE / 'freeeos_biquintic_progress/experiments/016_inversion_control'
TABLE = WORKSPACE / 'mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data'
SAVED = WORKSPACE / 'freeeos_biquintic_progress/experiments/001_native_prototype/matched_queries.npz'
TABLE_SHA = 'e6ca7ce27bee8197bb257a8c4886a06366e6951e58958fbed3735c877ee68502'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def own_use(logq, logt):
    return (logt >= 3.2) & ((logt >= 6.5) | (logq <= -.5) | (logq + 2 * logt - 12 <= -6.5))


class Model:
    def __init__(self, name, logq, logt, nodes=None, complete=None, split=False, c3=False):
        self.name, self.logq, self.logt = name, logq, logt
        self.nodes, self.complete, self.split, self.c3 = nodes, complete, split, c3

    def energy(self, logrho, temperature_log):
        temperature = np.exp(temperature_log)
        if self.nodes is None:
            radiation_energy = CRAD * temperature**4 / 10.0**logrho
            energy = 1.5e8 * temperature + radiation_energy
            jacobian = 1.5e8 * temperature + 4 * radiation_energy
            return energy, jacobian, jacobian / temperature
        sample_t = temperature_log / LN10
        sample_q = logrho - 2 * sample_t + 12
        evaluator = evaluate_c3 if self.c3 else evaluate
        native = evaluator(self.nodes, self.complete, self.logq, self.logt, sample_q, sample_t,
                           orders=LOW_ORDERS, add_radiation=self.split)
        state = thermodynamics(native, sample_q, sample_t)
        return state['E'], temperature * state['Cv'], state['Cv']


class Callback:
    def __init__(self, model, logrho):
        self.model, self.logrho = model, np.asarray(logrho)

    def __call__(self, temperature_log, indices):
        return self.model.energy(self.logrho[indices], temperature_log)

    def crossings(self, previous, following, indices):
        previous_t, following_t = previous / LN10, following / LN10
        previous_q = self.logrho[indices] - 2 * previous_t + 12
        following_q = self.logrho[indices] - 2 * following_t + 12
        previous_iq, previous_it, _ = query_cells(self.model.logq, self.model.logt, previous_q, previous_t)
        following_iq, following_it, _ = query_cells(self.model.logq, self.model.logt, following_q, following_t)
        return abs(following_iq - previous_iq) + abs(following_it - previous_it)


def build_models(logq, logt, rows):
    models = []
    common = np.ones((len(logq) - 1, len(logt) - 1), dtype=bool)
    for split in (False, True):
        low, valid, _ = source_constraints(rows, logq, logt, split_radiation=split)
        valid[520, :] = False
        spacing_q, spacing_t = (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10
        nodes, complete, _ = construct_twists(low, valid, spacing_q, spacing_t, radius=2, route='mean')
        models.append(Model('biquintic_split' if split else 'biquintic_total', logq, logt, nodes, complete, split))
        common &= valid_cells(complete)
        if split:
            nodes, complete, _ = construct_jets(low, valid, spacing_q, spacing_t, radius=2)
            models.append(Model('C3_split_control', logq, logt, nodes, complete, True, True))
            common &= valid_cells(complete)
    models.append(Model('analytic_ideal_gas_radiation', logq, logt))
    ratio = CRAD * (10.0**logt[None, :])**4 / (3 * 10.0**rows[..., 1])
    material = (ratio < 1e-3) & own_use(logq[:, None], logt[None, :])
    material &= (logq[:, None] >= -8) & (logq[:, None] <= 5) & (logt[None, :] <= 7.8)
    return models, common & valid_cells(material)


def bracket_path(logrho, lower, upper, logq, logt, cells):
    events = np.concatenate(([lower, upper], logt[(logt > lower) & (logt < upper)],
                             ((logrho + 12 - logq) / 2)))
    events = np.unique(events[(events >= lower) & (events <= upper)])
    probes = np.unique(np.concatenate((events, (events[:-1] + events[1:]) / 2)))
    sample_q = logrho - 2 * probes + 12
    density_index, temperature_index, inside = query_cells(logq, logt, sample_q, probes)
    if not np.all(inside & cells[density_index, temperature_index] & own_use(sample_q, probes)):
        return None
    return probes


def select_targets(models, cells, logq, logt, count):
    with np.load(SAVED) as saved:
        sample_q, sample_t = saved['logQ'].copy(), saved['logT'].copy()
        eligible = saved['common'] & saved['own_use']
    candidates = np.flatnonzero(eligible)
    generator = np.random.default_rng(20261028)
    candidates = generator.permutation(candidates)
    selected = []
    conditioning = []
    attempted = 0
    for sample_index in candidates:
        attempted += 1
        truth_t = sample_t[sample_index]
        logrho = sample_q[sample_index] + 2 * truth_t - 12
        lower, upper = truth_t - .30, truth_t + .08
        path = bracket_path(logrho, lower, upper, logq, logt, cells)
        if path is None:
            continue
        accepted = True
        conditions = []
        for model in models:
            energy, jacobian, cv = model.energy(np.full(path.size, logrho), path * LN10)
            condition = jacobian / energy
            accepted &= np.all(np.isfinite(energy) & np.isfinite(condition) & (energy > 0) & (cv > 0))
            accepted &= np.all((condition >= .05) & (condition <= 20))
            conditions.append(float(np.min(condition)))
        if accepted:
            selected.append((sample_index, logrho, truth_t, lower, upper))
            conditioning.append(conditions)
        if len(selected) == count:
            break
    if not selected:
        raise ValueError('No common conditioned material brackets')
    return np.asarray(selected), np.asarray(conditioning), attempted


def stats(values):
    finite = np.asarray(values)[np.isfinite(values)]
    if not finite.size:
        return {'finite': 0}
    return {'finite': int(finite.size), **dict(zip(['median', 'p99', 'max'],
                                                map(float, np.quantile(abs(finite), [.5, .99, 1]))))}


def run(args):
    if digest(TABLE) != TABLE_SHA:
        raise ValueError('Unexpected source-table SHA-256')
    header, logq, logt, rows = read_table(TABLE)
    models, cells = build_models(logq, logt, rows)
    selected, conditioning, attempted = select_targets(models, cells, logq, logt, args.targets)
    sample_ids, logrho, truth_t, lower_t, upper_t = selected.T
    lower, upper, truth = lower_t * LN10, upper_t * LN10, truth_t * LN10
    trial_records, results, jacobian_sanity = [], {}, {}
    for model in models:
        callback = Callback(model, logrho)
        target, _, _ = callback(truth, np.arange(truth.size))
        check_indices = np.arange(min(16, truth.size))
        _, exact_jacobian, _ = callback(truth[check_indices], check_indices)
        jacobian_sanity[model.name] = {}
        for step in (1e-5, 1e-6):
            plus = callback(truth[check_indices] + step, check_indices)[0]
            minus = callback(truth[check_indices] - step, check_indices)[0]
            jacobian_sanity[model.name][str(step)] = stats((plus - minus) / (2 * step * exact_jacobian) - 1)
        model_results = {}
        for start_name, starting in [('lower_cross_face', lower), ('upper_cross_face', upper)]:
            initial_crossings = callback.crossings(starting, truth, np.arange(truth.size))
            for safeguarded in (False, True):
                key = start_name + ('_safeguarded' if safeguarded else '_full_Newton')
                solved = solve(callback, target, lower, upper, starting, safeguarded,
                               rtol=args.rtol, max_iterations=args.max_iterations)
                success = solved['status'] == 'converged'
                temperature_error = np.expm1(solved['lnT'] - truth)
                unique_status, counts = np.unique(solved['status'], return_counts=True)
                model_results[key] = {
                    'trials': int(success.size), 'converged': int(success.sum()),
                    'status_counts': dict(zip(map(str, unique_status), map(int, counts))),
                    'converged_relative_T_error': stats(temperature_error[success]),
                    'converged_relative_energy_residual': stats(solved['relative_energy_residual'][success]),
                    'iterations_all': stats(solved['iterations']), 'iterations_converged': stats(solved['iterations'][success]),
                    'safeguards_total': int(solved['safeguards'].sum()), 'unsafe_Cv_total': int(solved['unsafe_Cv'].sum()),
                    'out_of_bracket_attempts_total': int(solved['out_of_bracket_attempts'].sum()),
                    'accepted_face_crossings_total': int(solved['accepted_face_crossings'].sum()),
                    'initial_to_truth_face_crossings': stats(initial_crossings),
                }
                for target_index in range(target.size):
                    trial_records.append({
                        'model': model.name, 'method': key, 'saved_sample_id': int(sample_ids[target_index]),
                        'logRho': logrho[target_index], 'truth_logT': truth_t[target_index],
                        'lower_logT': lower_t[target_index], 'upper_logT': upper_t[target_index],
                        'target_E_self_model': target[target_index],
                        'initial_to_truth_face_crossings': int(initial_crossings[target_index]),
                        'recovered_logT': solved['lnT'][target_index] / LN10,
                        'relative_T_error': temperature_error[target_index],
                        **{name: values[target_index].item() for name, values in solved.items() if name != 'lnT'},
                    })
        results[model.name] = model_results
    output = {
        'claim': 'Fixed-density self-inversion of each table model; NOT source accuracy or a stellar-solver test.',
        'table_sha256': TABLE_SHA, 'saved_queries_sha256': digest(SAVED), 'command': sys.argv,
        'evaluator_sha256': {'freeeos_potential.py': digest(BIQUINTIC / 'freeeos_potential.py'),
                             'c3_potential.py': digest(BIQUINTIC / 'c3_control/c3_potential.py')},
        'coordinates': 'saved off-node truth; variable lnT with fixed log10rho, q=log10rho-2*log10T+12',
        'jacobian': 'analytic dE/dlnT=T*Cv from the same model; no finite-difference Jacobian',
        'jacobian_sanity_16_truth_points': jacobian_sanity,
        'selected_coordinates': int(selected.shape[0]), 'model_specific_targets': int(len(models) * selected.shape[0]),
        'candidates_attempted': attempted, 'model_names': [model.name for model in models],
        'conditioning_min_TCv_over_E_at_checked_path_points': np.min(conditioning, axis=0).tolist(),
        'selection': {'material_Prad_over_Pgas_max': 1e-3, 'q_range': [-8, 5], 'logT_max': 7.8,
                      'truth_to_lower_dex': .30, 'truth_to_upper_dex': .08,
                      'validity': 'all crossed cells and faces; common native stencils, own-use, material-node masks; iq520 excluded',
                      'Cv_screening': 'positive finite at all faces and segment midpoints; not an analytic positivity proof'},
        'solver': {'relative_energy_tolerance': args.rtol, 'max_iterations': args.max_iterations,
                   'safeguard': 'bracket sign update; accept positive-Jacobian Newton inside bracket with step <=75% of bracket width, else midpoint bisection',
                   'full_Newton': 'no damping/bisection; terminate if step leaves the declared valid bracket, without querying unverified cells'},
        'results': results,
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    products = [OUTPUT / name for name in ['results.json', 'trials.csv', 'targets.npz', 'SUMMARY.md']]
    if not args.overwrite and any(path.exists() for path in products):
        raise ValueError('Products exist; use --overwrite to repeat this sidecar')
    products[0].write_text(json.dumps(output, indent=2, allow_nan=False) + '\n')
    with products[1].open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(trial_records[0]))
        writer.writeheader()
        writer.writerows(trial_records)
    np.savez(products[2], saved_sample_id=sample_ids.astype(int), logRho=logrho, truth_logT=truth_t,
             lower_logT=lower_t, upper_logT=upper_t, conditioning=conditioning)
    lines = ['# Bounded fixed-density inversion control', '', output['claim'], '',
             f'{selected.shape[0]} saved coordinates; {len(models) * selected.shape[0]} model-specific self targets.',
             'Analytic Jacobian T*Cv; brackets cross q/T faces; no source reference or stellar solve.', '',
             '| Model | Start/method | Converged | Iteration median / max | Safeguards | Unsafe Cv | Bracket escapes | T error max (converged) |',
             '| --- | --- | ---: | --- | ---: | ---: | ---: | ---: |']
    for model_name, model_results in results.items():
        for method, result in model_results.items():
            iterations = result['iterations_all']
            error = result['converged_relative_T_error'].get('max', float('nan'))
            lines.append(f"| {model_name} | {method} | {result['converged']}/{result['trials']} | "
                         f"{iterations['median']:g} / {iterations['max']:g} | {result['safeguards_total']} | "
                         f"{result['unsafe_Cv_total']} | {result['out_of_bracket_attempts_total']} | {error:.3e} |")
    lines.extend(['', 'Full-step escape means leaving this declared valid bracket, not proven global divergence.',
                  'Targets are generated independently for each same model; recovered-temperature errors are self-inversion only.',
                  'Face/midpoint screening does not prove monotonicity between every checked point; run-time Cv guards remain active.'])
    products[3].write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--targets', type=int, default=100)
    parser.add_argument('--rtol', type=float, default=1e-10)
    parser.add_argument('--max-iterations', type=int, default=60)
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    if not 1 <= args.targets <= 125 or not 0 < args.rtol < 1 or not 1 <= args.max_iterations <= 100:
        parser.error('Require 1..125 coordinates (<=500 self targets), positive tolerance, and 1..100 iterations')
    run(args)


if __name__ == '__main__':
    main()
