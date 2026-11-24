"""Exploratory printed-decimal rounding sensitivity of the material prototype."""

import argparse
from decimal import Decimal
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from freeeos_potential import (
    CRAD, LN10, checksum, construct_twists, evaluate, local_derivative,
    query_cells, read_table, source_constraints, thermodynamic_partials, valid_cells,
)
from run_audit import own_use
from validate_table import material_nodes


FIELDS = {1: "logPgas", 2: "logE", 3: "logS", 4: "chiRho", 5: "chiT",
          6: "Cp", 7: "Cv", 8: "dE_dRho", 9: "dS_dT", 10: "dS_dRho"}
GROUPS = {"energy_entropy": (2, 3), "six_constraints": (1, 2, 3, 4, 5, 7),
          "constraints_and_redundant": tuple(FIELDS)}
VALUES = ("P", "Pgas", "E", "S", "Cv", "Cp", "chiRho", "chiT",
          "gamma1", "gamma3", "grad_ad")
RESPONSES = ("Cv", "Cp", "chiRho", "chiT", "gamma1", "gamma3", "grad_ad")
METRICS = VALUES + tuple(field + suffix for field in RESPONSES for suffix in ("_x", "_t"))
ROUTES = (((0, 0), 0, (1, 0)), ((0, 0), 1, (0, 1)),
          ((1, 0), 0, (2, 0)), ((1, 0), 1, (1, 1)),
          ((0, 1), 0, (1, 1)), ((0, 1), 1, (0, 2)))


def distribution(values):
    values = np.abs(np.asarray(values))
    if not values.size or not np.all(np.isfinite(values)):
        raise ValueError("Expected finite, nonempty diagnostic values")
    return dict(zip(("median", "p90", "p99", "max"),
                    map(float, np.quantile(values, (.5, .9, .99, 1)))))


def printed_interval(token):
    center = Decimal(token)
    if not center.is_finite():
        raise ValueError("Nonfinite printed value")
    half = Decimal(5).scaleb(center.as_tuple().exponent - 1)
    lower_decimal, upper_decimal = center - half, center + half
    lower, upper = float(lower_decimal), float(upper_decimal)
    if Decimal.from_float(lower) < lower_decimal:
        lower = np.nextafter(lower, np.inf)
    if Decimal.from_float(upper) > upper_decimal:
        upper = np.nextafter(upper, -np.inf)
    if not lower <= float(center) <= upper:
        raise ValueError("Printed interval does not contain parsed source value")
    return float(half), lower, upper


def read_intervals(path, rows, support):
    indices = np.flatnonzero(support)
    half = np.empty((len(indices), len(FIELDS)))
    lower, upper = np.empty_like(half), np.empty_like(half)
    records_seen, selected_seen = 0, 0
    examples = {}
    with Path(path).open() as stream:
        next(stream)
        next(stream)
        for line in stream:
            tokens = line.split()
            if len(tokens) != 22 or tokens[0][0] not in "+-0123456789":
                continue
            if selected_seen < len(indices) and records_seen == indices[selected_seen]:
                for position, column in enumerate(FIELDS):
                    half[selected_seen, position], lower[selected_seen, position], upper[selected_seen, position] = printed_interval(tokens[column])
                    if float(tokens[column]) != rows.reshape(-1, 22)[records_seen, column]:
                        raise ValueError("String records are not aligned with core table reader")
                    if selected_seen == 0:
                        examples[FIELDS[column]] = {"token": tokens[column], "half_last_decimal": float(half[0, position]),
                                                   "lower_float_inward": float(lower[0, position]), "upper_float_inward": float(upper[0, position])}
                selected_seen += 1
            records_seen += 1
    if records_seen != rows.shape[0] * rows.shape[1] or selected_seen != len(indices):
        raise ValueError("Incomplete source interval scan")
    return indices, half, lower, upper, examples


def perturb(rows, indices, lower, upper, columns, seed):
    changed = rows.copy()
    flattened = changed.reshape(-1, 22)
    for position, column in enumerate(FIELDS):
        if column in columns:
            generator = np.random.default_rng(np.random.SeedSequence([seed, column]))
            candidate = generator.uniform(lower[:, position], upper[:, position])
            flattened[indices, column] = np.clip(candidate, lower[:, position], upper[:, position])
            if np.any(flattened[indices, column] < lower[:, position]) or np.any(flattened[indices, column] > upper[:, position]):
                raise AssertionError("Perturbation outside printed interval")
    preserved = [column for column in range(22) if column not in columns]
    np.testing.assert_array_equal(changed[:, :, preserved], rows[:, :, preserved])
    return changed


def face_queries(record, logq, logt):
    result = {}
    for direction, axis in (("u", 0), ("tau", 1)):
        metrics = record["faces"][direction]["response_log_coordinate_jacobian_jump_scaled_by_response_or_partial"]
        worst = metrics["Cp_t" if axis == 0 else "Cv_t"]["worst"]
        sample_q, sample_t = np.array([worst["logQ"]]), np.array([worst["logT"]])
        density_index, temperature_index, inside = query_cells(logq, logt, sample_q, sample_t)
        if not np.all(inside):
            raise ValueError("Recorded face outside source grid")
        boundary = int(np.argmin(abs((logq if axis == 0 else logt) - (sample_q[0] if axis == 0 else sample_t[0]))))
        if axis == 0:
            if not np.isclose(sample_q[0], logq[boundary], rtol=0, atol=1e-12):
                raise ValueError("Recorded u face is not a grid face")
            sample_q[0] = logq[boundary]
            right = np.array([boundary]), temperature_index
            left = np.array([boundary - 1]), temperature_index
        else:
            if not np.isclose(sample_t[0], logt[boundary], rtol=0, atol=1e-12):
                raise ValueError("Recorded tau face is not a grid face")
            sample_t[0] = logt[boundary]
            right = density_index, np.array([boundary])
            left = density_index, np.array([boundary - 1])
        result[direction] = {"q": sample_q, "t": sample_t, "left": left, "right": right,
                             "recorded_worst": worst, "target": "Cp_t" if axis == 0 else "Cv_t"}
    return result


def selected_queries(samples, complete, rows, logq, logt, count):
    density_index, temperature_index, inside = query_cells(logq, logt, samples["logQ"], samples["logT"])
    accepted = valid_cells(complete) & valid_cells(material_nodes(rows, logq, logt))
    selected = inside & accepted[density_index, temperature_index] & own_use(samples["logQ"], samples["logT"])
    candidates = np.flatnonzero(selected)
    if len(candidates) < count:
        raise ValueError("Insufficient clean saved material queries")
    choice = np.sort(np.random.default_rng(6006).choice(candidates, count, replace=False))
    return choice, samples["logQ"][choice], samples["logT"][choice]


def nodal_selection(logq, logt, sample_q, sample_t, faces):
    selected = np.zeros((len(logq), len(logt)), dtype=bool)
    density_index, temperature_index, _ = query_cells(logq, logt, sample_q, sample_t)
    cells = [(density_index, temperature_index)]
    cells.extend(face[side] for face in faces.values() for side in ("left", "right"))
    for density_index, temperature_index in cells:
        for density_end in (0, 1):
            for temperature_end in (0, 1):
                selected[density_index + density_end, temperature_index + temperature_end] = True
    return selected


def stencil_support(selected, valid):
    support = np.zeros_like(selected)
    for density_offset in range(-4, 5):
        for temperature_offset in range(-4, 5):
            destination = (slice(max(0, density_offset), min(selected.shape[0], selected.shape[0] + density_offset)),
                           slice(max(0, temperature_offset), min(selected.shape[1], selected.shape[1] + temperature_offset)))
            source = (slice(max(0, -density_offset), min(selected.shape[0], selected.shape[0] - density_offset)),
                      slice(max(0, -temperature_offset), min(selected.shape[1], selected.shape[1] - temperature_offset)))
            support[destination] |= selected[source]
    return support & valid


def states_at(nodes, complete, logq, logt, sample_q, sample_t, faces, split):
    queries = thermodynamic_partials(evaluate(nodes, complete, logq, logt, sample_q, sample_t, add_radiation=split), sample_q, sample_t)
    face_states = {direction: {side: thermodynamic_partials(
        evaluate(nodes, complete, logq, logt, face["q"], face["t"], cell_indices=face[side], add_radiation=split), face["q"], face["t"])
        for side in ("left", "right")} for direction, face in faces.items()}
    for state in [queries] + [state for sides in face_states.values() for state in sides.values()]:
        if not all(np.all(np.isfinite(state[field])) for field in METRICS):
            raise ValueError("Nonfinite perturbed thermodynamic output")
    return queries, face_states


def defects(low, valid, selected, spacings):
    result = {}
    for base, axis, target in ROUTES:
        derivative, accepted = local_derivative(low[base], valid, spacings[axis], axis, radius=2)
        chosen = selected & accepted
        name = f"d{base}_axis{axis}_vs_{target}"
        result[name] = {"values": derivative[chosen] - low[target][chosen], "chosen": chosen,
                        "scale": np.maximum(np.maximum(abs(derivative[chosen]), abs(low[target][chosen])), np.finfo(float).tiny)}
    return result


def redundant_closure(rows, source, selected):
    density, temperature = source["rho"][selected], source["T"][selected]
    pressure_scale = source["P"][selected] / density
    mixed = pressure_scale * rows[:, :, 5][selected]
    pairs = {"mixed_energy": (mixed, pressure_scale - density * rows[:, :, 8][selected]),
             "mixed_entropy": (mixed, -density * temperature * rows[:, :, 10][selected]),
             "Cv_entropy": (rows[:, :, 7][selected], temperature * rows[:, :, 9][selected])}
    return {name: {"values": first - second,
                   "scale": np.maximum(np.maximum(abs(first), abs(second)), np.finfo(float).tiny)}
            for name, (first, second) in pairs.items()}


def variability(values, baseline, scale, seeds, sample_q, sample_t):
    changes = values - baseline
    scaled = changes / scale
    seed_index, query_index = np.unravel_index(np.argmax(abs(changes)), changes.shape)
    scaled_seed, scaled_query = np.unravel_index(np.argmax(abs(scaled)), changes.shape)
    return {"absolute_change": distribution(changes), "field_scaled_change": distribution(scaled),
            "worst_field_scaled_change": {"seed": seeds[scaled_seed], "logQ": float(sample_q[scaled_query]),
                "logT": float(sample_t[scaled_query]), "baseline": float(baseline[scaled_query]),
                "perturbed": float(values[scaled_seed, scaled_query]), "signed_change": float(changes[scaled_seed, scaled_query]),
                "scale": float(scale[scaled_query]), "signed_scaled_change": float(scaled[scaled_seed, scaled_query])},
            "worst_absolute_change": {"seed": seeds[seed_index], "logQ": float(sample_q[query_index]),
                "logT": float(sample_t[query_index]), "baseline": float(baseline[query_index]),
                "perturbed": float(values[seed_index, query_index]), "signed_change": float(changes[seed_index, query_index])}}


def compare_face(realizations, baseline, field, seeds):
    first, second = float(baseline["left"][field][0]), float(baseline["right"][field][0])
    scale = max(abs(first), abs(second), abs(float(baseline["left"][field.split("_x")[0].split("_t")[0]][0])),
                abs(float(baseline["right"][field.split("_x")[0].split("_t")[0]][0])), np.finfo(float).tiny)
    left = np.array([sides["left"][field][0] for sides in realizations])
    right = np.array([sides["right"][field][0] for sides in realizations])
    jump, changes = first - second, left - right - (first - second)
    maximum_noise = float(np.max(abs(changes)))
    worst = int(np.argmax(abs(changes)))
    ratio = maximum_noise / abs(jump) if jump != 0 else None
    return {"baseline_left": first, "baseline_right": second, "baseline_signed_jump": jump,
            "baseline_field_scaled_jump": jump / scale, "field_scale": scale,
            "left_absolute_change": distribution(left - first), "right_absolute_change": distribution(right - second),
            "jump_absolute_change": distribution(changes), "jump_field_scaled_change": distribution(changes / scale),
            "worst_jump_change": {"seed": seeds[worst], "left": float(left[worst]), "right": float(right[worst]),
                                  "signed_jump": float(left[worst] - right[worst]), "signed_jump_change": float(changes[worst])},
            "max_jump_noise_over_abs_baseline_jump": ratio, "perturbed_jump_range": [float(np.min(left - right)), float(np.max(left - right))],
            "sampled_rounding_negligible_vs_jump": ratio is not None and ratio < .01}


def main():
    workspace = Path(__file__).resolve().parents[5]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=workspace / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data")
    parser.add_argument("--queries", type=Path, default=workspace / "freeeos_report/sample_residuals.csv")
    parser.add_argument("--faces", type=Path, default=workspace / "freeeos_biquintic_progress/experiments/002_table_validation/results.json")
    parser.add_argument("--output", type=Path, default=workspace / "freeeos_biquintic_progress/experiments/006_rounding_sensitivity")
    parser.add_argument("--seeds", type=int, default=8)
    parser.add_argument("--queries-count", type=int, default=1000)
    parser.add_argument("--split-radiation", action="store_true")
    arguments = parser.parse_args()
    allowed = (workspace / "freeeos_biquintic_progress/experiments/006_rounding_sensitivity").resolve()
    if arguments.output.resolve() != allowed:
        parser.error("Outputs are restricted to experiments/006_rounding_sensitivity")
    if arguments.seeds < 2 or arguments.queries_count < 1:
        parser.error("Require at least two seeds and one query")
    core_path = Path(__file__).resolve().parents[1] / "freeeos_potential.py"
    core_hash = checksum(core_path)
    header, logq, logt, rows = read_table(arguments.source)
    record = json.loads(arguments.faces.read_text())
    if checksum(arguments.source) != record["provenance"]["source_sha256"]:
        raise ValueError("Face record source differs from source table")
    spacings = (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10
    low, valid, source = source_constraints(rows, logq, logt, split_radiation=arguments.split_radiation)
    nodes, complete, _ = construct_twists(low, valid, *spacings, radius=2, route="mean")
    faces = face_queries(record, logq, logt)
    clean_cells = valid_cells(complete) & valid_cells(material_nodes(rows, logq, logt))
    for face in faces.values():
        if not own_use(face["q"], face["t"]).all() or not all(clean_cells[face[side]].all() for side in ("left", "right")):
            raise ValueError("Requested face is not a clean material own-use face")
    samples = np.genfromtxt(arguments.queries, delimiter=",", names=True)
    query_indices, sample_q, sample_t = selected_queries(samples, complete, rows, logq, logt, arguments.queries_count)
    selected = nodal_selection(logq, logt, sample_q, sample_t, faces)
    support = stencil_support(selected, valid)
    indices, half, lower, upper, examples = read_intervals(arguments.source, rows, support)
    baseline, baseline_faces = states_at(nodes, complete, logq, logt, sample_q, sample_t, faces, arguments.split_radiation)
    baseline_defects = defects(low, valid, selected, spacings)
    baseline_closure = redundant_closure(rows, source, selected)
    seeds = list(range(2026102800, 2026102800 + arguments.seeds))
    results = {"status": "exploratory printed-table rounding, not solver/source accuracy or a rigorous interval bound",
        "provenance": {"source_sha256": checksum(arguments.source), "queries_sha256": checksum(arguments.queries),
            "face_record_sha256": checksum(arguments.faces), "face_record_preprocessing_sha256": record["provenance"]["preprocessing_sha256"],
            "preprocessing_sha256": core_hash, "script_sha256": checksum(__file__), "numpy": np.__version__,
            "python": sys.version, "composition": {"X": float(header[1]), "Z": float(header[2])}, "radiation_constant": float(CRAD), "command": sys.argv},
        "model": {"interval": "printed decimal value +/- half of 10**Decimal(token).exponent; logarithmic fields stay in log10 units",
            "sampling": "independent uniform draws between inward-rounded printed-decimal endpoints per source row/column; binary64 rounding, with endpoint clipping as a safety guard",
            "limitations": "exploratory independent rounding, not rigorous bounds, source precision proof, or correlated physical perturbations; no coordinate/provenance perturbation",
            "seeds": seeds, "groups": {name: [FIELDS[column] for column in columns] for name, columns in GROUPS.items()},
            "unconsumed_but_perturbed": ["Cp", "dE_dRho", "dS_dT", "dS_dRho"], "scheme": "radius2_mean", "split_radiation": arguments.split_radiation,
            "perturbed_valid_nodes": int(support.sum()), "support": "valid native nodes in conservative radius4 square around requested cell corners, covering radius2 two-stage twists",
            "support_nodes_outside_material_query_cut": int((support & ~material_nodes(rows, logq, logt)).sum()),
            "interval_examples": examples, "half_interval_distributions_stored_units": {field: distribution(half[:, position]) for position, field in enumerate(FIELDS.values())},
            "scaling": "max(|baseline field|,|baseline derivative|), no unit floor; face scales include both sides; absolute changes also reported",
            "negligible_threshold": "max sampled change of intercell jump <1% of baseline jump; exploratory comparison, not an accuracy tolerance"},
        "queries": {"count": len(sample_q), "saved_csv_indices": query_indices.tolist(), "selection_seed": 6006,
            "domain": "saved native complete material cells, q[-8,5], logT[3.2,7.8], corner Prad/Pgas<1e-3, own-use cuts"},
        "f_conditioning": {"cancellation_factor": distribution((abs(source["E"][support]) + abs(source["T"][support] * source["S"][support])) / np.maximum(abs(source["E"][support] - source["T"][support] * source["S"][support]), np.finfo(float).tiny)),
            "linearized_f_half_interval_sum_erg_per_g": distribution(LN10 * (source["E"][support] * half[:, 1] + source["T"][support] * source["S"][support] * half[:, 2]))},
        "groups": {}}
    for direction, face in faces.items():
        if not arguments.split_radiation:
            np.testing.assert_allclose([baseline_faces[direction][side][face["target"]][0] for side in ("left", "right")],
                                       [face["recorded_worst"][side] for side in ("left", "right")], rtol=1e-10, atol=1e-5)
    previous_outputs, previous_faces = {}, {}
    for name, columns in GROUPS.items():
        query_runs, face_runs = [], {direction: [] for direction in faces}
        defect_runs = {route: [] for route in baseline_defects}
        closure_runs = {route: [] for route in baseline_closure}
        for seed in seeds:
            changed = perturb(rows, indices, lower, upper, columns, seed)
            changed_low, changed_valid, changed_source = source_constraints(changed, logq, logt, split_radiation=arguments.split_radiation)
            np.testing.assert_array_equal(changed_valid, valid)
            changed_nodes, changed_complete, _ = construct_twists(changed_low, changed_valid, *spacings, radius=2, route="mean")
            np.testing.assert_array_equal(changed_complete, complete)
            query_state, face_states = states_at(changed_nodes, changed_complete, logq, logt, sample_q, sample_t, faces, arguments.split_radiation)
            query_runs.append(np.stack([query_state[field] for field in METRICS]))
            for direction in faces:
                face_runs[direction].append(face_states[direction])
            for route, data in defects(changed_low, changed_valid, selected, spacings).items():
                np.testing.assert_array_equal(data["chosen"], baseline_defects[route]["chosen"])
                defect_runs[route].append(data["values"])
            for route, data in redundant_closure(changed, changed_source, selected).items():
                closure_runs[route].append(data["values"])
        query_values = np.stack(query_runs)
        face_values = np.array([[[[state[side][field][0] for field in METRICS] for side in ("left", "right")]
                                for state in face_runs[direction]] for direction in faces])
        if name == "constraints_and_redundant":
            np.testing.assert_array_equal(query_values, previous_outputs["six_constraints"])
            np.testing.assert_array_equal(face_values, previous_faces["six_constraints"])
        previous_outputs[name] = query_values
        previous_faces[name] = face_values
        group = {"masks_preserved_all_seeds": True, "query_variability": {}, "faces": {}, "nodal_fd_integrability": {},
                 "redundant_printed_closure": {}}
        for position, field in enumerate(METRICS):
            base_field = field[:-2] if field.endswith(("_x", "_t")) else field
            scale = np.maximum(np.maximum(abs(baseline[field]), abs(baseline[base_field])), np.finfo(float).tiny)
            group["query_variability"][field] = variability(query_values[:, position], baseline[field], scale, seeds, sample_q, sample_t)
        for direction, face in faces.items():
            metrics = {field: compare_face(face_runs[direction], baseline_faces[direction], field, seeds) for field in METRICS}
            target = metrics[face["target"]]
            group["faces"][direction] = {"logQ": float(face["q"][0]), "logT": float(face["t"][0]),
                "recorded_target": face["target"], "recorded_worst": face["recorded_worst"], "metrics": metrics,
                "assessment": "sampled print rounding does not explain the recorded response jump" if target["sampled_rounding_negligible_vs_jump"] else "sampled rounding is not negligible; attribution inconclusive"}
        for route, data in baseline_defects.items():
            nodal_q, nodal_t = np.broadcast_arrays(logq[:, None], logt[None, :])
            noise = variability(np.stack(defect_runs[route]), data["values"], data["scale"], seeds, nodal_q[data["chosen"]], nodal_t[data["chosen"]])
            baseline_max = float(np.max(abs(data["values"])))
            worst = int(np.argmax(abs(data["values"])))
            group["nodal_fd_integrability"][route] = {"baseline_absolute_defect": distribution(data["values"]),
                "baseline_scaled_defect": distribution(data["values"] / data["scale"]), "rounding_variability": noise,
                "baseline_worst": {"logQ": float(nodal_q[data["chosen"]][worst]), "logT": float(nodal_t[data["chosen"]][worst]),
                                   "signed_defect": float(data["values"][worst]), "scale": float(data["scale"][worst])},
                "max_noise_over_max_baseline_defect": noise["absolute_change"]["max"] / baseline_max if baseline_max else None,
                "limitation": "centered radius2 finite-grid derivative defect includes truncation; not a direct-source integrability error"}
        for route, data in baseline_closure.items():
            group["redundant_printed_closure"][route] = {"baseline_scaled_closure_defect": distribution(data["values"] / data["scale"]),
                "rounding_variability": variability(np.stack(closure_runs[route]), data["values"], data["scale"], seeds,
                                                   nodal_q[selected], nodal_t[selected]),
                "limitation": "independent printed derivative closure, not direct-source accuracy; redundant columns do not constrain the interpolant"}
        results["groups"][name] = group
        print(name, {direction: group["faces"][direction]["metrics"][faces[direction]["target"]]["max_jump_noise_over_abs_baseline_jump"] for direction in faces}, flush=True)
    results["redundant_columns_no_effect_asserted"] = True
    results["assessment"] = {"all_target_faces_sampled_rounding_below_one_percent_of_jump": all(
        group["faces"][direction]["metrics"][faces[direction]["target"]]["sampled_rounding_negligible_vs_jump"]
        for group in results["groups"].values() for direction in faces),
        "scope": "sampled independent rounding cannot account for target jumps if the observed gap passes; not exclusion of all possible correlated interval perturbations"}
    if checksum(core_path) != core_hash:
        raise RuntimeError("Core changed during experiment; discard this mixed snapshot")
    arguments.output.mkdir(parents=True, exist_ok=True)
    (arguments.output / "results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")
    np.savez_compressed(arguments.output / "samples.npz", logQ=sample_q, logT=sample_t, csv_indices=query_indices,
                        seeds=seeds, metrics=np.array(METRICS), face_directions=np.array(list(faces)),
                        **previous_outputs, **{name + "_faces": values for name, values in previous_faces.items()})
    lines = ["# Printed rounding sensitivity", "", "Independent uniform half-last-decimal perturbations; exploratory, not a rigorous bound or source-accuracy test.", "",
             f"{len(sample_q)} saved clean material queries; {len(seeds)} seeds per group; radius2_mean; split radiation={arguments.split_radiation}.", "", "## Recorded material faces"]
    for name, group in results["groups"].items():
        for direction, face in group["faces"].items():
            comparison = face["metrics"][face["recorded_target"]]
            lines.append(f"- {name}, {direction} ({face['logQ']:.12g}, {face['logT']:.12g}), {face['recorded_target']}: baseline jump {comparison['baseline_signed_jump']:.6g}; maximum sampled jump noise {comparison['jump_absolute_change']['max']:.6g}; ratio {comparison['max_jump_noise_over_abs_baseline_jump']:.3g}. {face['assessment']}.")
    lines.extend(["", "## Saved-query conditioning", "", "Maximum field-scaled variability over the 1,000 queries and all seeds (six selected constraints):"])
    for field in ("P", "Cv", "chiRho", "chiT", "Cv_t", "Cp_t"):
        metric = results["groups"]["six_constraints"]["query_variability"][field]
        worst = metric["worst_field_scaled_change"]
        lines.append(f"- {field}: maximum field-scaled change {metric['field_scaled_change']['max']:.6g} at ({worst['logQ']:.12g}, {worst['logT']:.12g}), seed {worst['seed']}; maximum absolute change {metric['absolute_change']['max']:.6g} (absolute extrema are recorded separately).")
    conclusion = "Target face jumps exceed sampled print-rounding changes by more than 100-fold." if results["assessment"]["all_target_faces_sampled_rounding_below_one_percent_of_jump"] else "The sampled changes are not negligible at every target face; attribution is inconclusive."
    lines.extend(["", "## Nodal defects", ""])
    for route, data in results["groups"]["six_constraints"]["nodal_fd_integrability"].items():
        lines.append(f"- {route}: maximum sampled defect change / maximum baseline defect = {data['max_noise_over_max_baseline_defect']:.6g}. This compares different extrema, not a pointwise bound.")
    lines.extend(["", "The additional redundant/Cp columns produce bitwise-identical query and face outputs to six-constraint perturbations on shared seeds, while independent printed-derivative closure is audited separately. Masks/coordinates/provenance remain unchanged.",
                  "", conclusion + " This does not distinguish C2 interpolation/truncation, missing twists, or source inter-node consistency. Nodal finite-difference defects also include truncation. No solver/default conclusion follows.",
                  "", "Full distributions, worst coordinates/seeds, interval examples, fingerprints, and scaling definitions are in results.json; samples.npz stores the query realizations."])
    (arguments.output / "SUMMARY.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
