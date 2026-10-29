"""Native-only shifted quartic twist/coverage CONTROL; not a runtime fallback."""

import os

os.environ["OPENBLAS_NUM_THREADS"] = "8"

import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from active_faces import classifications
from freeeos_potential import (
    LN10, checksum, construct_twists, evaluate, query_cells, read_table,
    source_constraints, thermodynamic_partials, valid_cells,
)
from validate_table import material_nodes, table_state


VALUES = ("P", "Pgas", "E", "S", "Cv", "Cp", "chiRho", "chiT", "gamma1", "gamma3", "grad_ad")
JACOBIANS = tuple(field + suffix for field in ("Cv", "Cp", "chiRho", "chiT", "gamma1", "gamma3", "grad_ad") for suffix in ("_x", "_t"))


def distribution(values):
    values = np.abs(np.asarray(values))
    finite = values[np.isfinite(values)]
    result = {"count": int(values.size), "nonfinite_count": int(values.size - finite.size)}
    if finite.size:
        result.update(dict(zip(("median", "p99", "max"), map(float, np.quantile(finite, (.5, .99, 1))))))
    else:
        result.update({key: None for key in ("median", "p99", "max")})
    return result


def state_distribution(values, sample_q, sample_t):
    result = distribution(values)
    finite = np.isfinite(values)
    for name, sign in (("signed_minimum", 1), ("signed_maximum", -1)):
        if finite.any():
            index = int(np.argmin(np.where(finite, sign * values, np.inf)))
            result[name] = {"value": float(values[index]), "logQ": float(sample_q[index]), "logT": float(sample_t[index])}
    return result


def original_bounds(shape):
    density_index, temperature_index = np.indices(shape)
    return np.stack((density_index, density_index, temperature_index, temperature_index), axis=-1)


def union_bounds(first, second):
    return np.stack((np.minimum(first[..., 0], second[..., 0]), np.maximum(first[..., 1], second[..., 1]),
                     np.minimum(first[..., 2], second[..., 2]), np.maximum(first[..., 3], second[..., 3])), axis=-1)


def polynomial_derivative(values, valid, spacing, axis, shifted=True, donor_bounds=None):
    if values.shape != valid.shape or values.ndim != 2 or axis not in (0, 1) or not np.isfinite(spacing) or spacing <= 0:
        raise ValueError("Require matching 2D arrays, axis0/1 and positive spacing")
    result = np.full(values.shape, np.nan)
    accepted = np.zeros_like(valid)
    start = np.full(values.shape, -1, dtype=np.int16)
    position_record = np.full(values.shape, -1, dtype=np.int8)
    support = np.full((*values.shape, 4), -1, dtype=np.int16)
    bounds = original_bounds(values.shape) if donor_bounds is None else donor_bounds
    length = values.shape[axis]
    if length < 5:
        return result, accepted, {"start": start, "position": position_record, "support": support}
    node_indices = np.indices(values.shape)[axis]
    for position in ((2, 1, 3, 0, 4) if shifted else (2,)):
        offsets = np.arange(5) - position
        vandermonde = offsets[None, :] ** np.arange(5)[:, None]
        target = np.zeros(5)
        target[1] = 1.
        weights = np.linalg.solve(vandermonde.astype(float), target) / spacing
        interior = [slice(None), slice(None)]
        interior[axis] = slice(position, length - 4 + position)
        interior = tuple(interior)
        safe = ~accepted[interior]
        total = np.zeros_like(values[interior])
        combined = None
        for offset, weight in zip(offsets, weights):
            donors = [slice(None), slice(None)]
            donors[axis] = slice(position + offset, length - 4 + position + offset)
            donors = tuple(donors)
            safe &= valid[donors] & np.isfinite(values[donors])
            if abs(weight) > 1e-14 / spacing:
                total += weight * np.where(valid[donors], values[donors], 0.)
            combined = bounds[donors].copy() if combined is None else union_bounds(combined, bounds[donors])
        result[interior] = np.where(safe, total, result[interior])
        accepted[interior] |= safe
        start[interior] = np.where(safe, node_indices[interior] - position, start[interior])
        position_record[interior] = np.where(safe, position, position_record[interior])
        support[interior] = np.where(safe[..., None], combined, support[interior])
    return result, accepted, {"start": start, "position": position_record, "support": support}


def build_jets(low, valid, du, dtau, shifted=True):
    nodes = low.copy()
    routes, stencils = {}, {}
    first21, mask21a, stencils["21_first"] = polynomial_derivative(low[2, 0], valid, dtau, 1, shifted)
    second21, mask21b, stencils["21_second"] = polynomial_derivative(low[1, 1], valid, du, 0, shifted)
    first12, mask12a, stencils["12_first"] = polynomial_derivative(low[1, 1], valid, dtau, 1, shifted)
    second12, mask12b, stencils["12_second"] = polynomial_derivative(low[0, 2], valid, du, 0, shifted)
    twist_valid = mask21a & mask21b & mask12a & mask12b
    nodes[2, 1], nodes[1, 2] = (first21 + second21) / 2, (first12 + second12) / 2
    twist_bounds = stencils["21_first"]["support"]
    for name in ("21_second", "12_first", "12_second"):
        twist_bounds = union_bounds(twist_bounds, stencils[name]["support"])
    twist_bounds[~twist_valid] = -1
    first22, mask22a, stencils["22_first"] = polynomial_derivative(nodes[2, 1], twist_valid, dtau, 1, shifted, twist_bounds)
    second22, mask22b, stencils["22_second"] = polynomial_derivative(nodes[1, 2], twist_valid, du, 0, shifted, twist_bounds)
    complete = twist_valid & mask22a & mask22b
    nodes[2, 2] = (first22 + second22) / 2
    support = union_bounds(twist_bounds, union_bounds(stencils["22_first"]["support"], stencils["22_second"]["support"]))
    support[~complete] = -1
    nodes[:, :, ~complete] = np.nan
    for orders, first, second, mask in (((2, 1), first21, second21, twist_valid),
                                      ((1, 2), first12, second12, twist_valid), ((2, 2), first22, second22, complete)):
        routes[str(orders)] = (first, second, mask)
    return nodes, complete, {"stencils": stencils, "twist_valid": twist_valid, "support": support, "routes": routes}


def reach_report(complete, metadata, logq, logt):
    density, temperature = np.indices(complete.shape)
    support = metadata["support"]
    density_reach = np.maximum(density - support[..., 0], support[..., 1] - density)[complete]
    temperature_reach = np.maximum(temperature - support[..., 2], support[..., 3] - temperature)[complete]
    if (density_reach > 8).any() or (temperature_reach > 8).any():
        raise AssertionError("Two-stage donor reach exceeded eight native spacings")
    if np.min(logt[support[..., 2][complete]]) < 3.2 - 1e-12:
        raise AssertionError("A donor below the source floor entered the control")
    return {"max_density_node_reach": int(np.max(density_reach)), "max_temperature_node_reach": int(np.max(temperature_reach)),
            "max_logQ_reach": float(np.max(density_reach) * (logq[1] - logq[0])),
            "max_logT_reach": float(np.max(temperature_reach) * (logt[1] - logt[0])),
            "minimum_donor_logT": float(np.min(logt[support[..., 2][complete]])),
            "support_interpretation": "bounding boxes of numerical and eligibility donors; not every point inside the box is used; exact windows/masks are archived",
            "stencils": {name: {"accepted": int(np.sum(data["position"] >= 0)), "shifted": int(np.sum((data["position"] >= 0) & (data["position"] != 2))),
                                "query_position_counts": {str(position): int(np.sum(data["position"] == position)) for position in range(5)}}
                         for name, data in metadata["stencils"].items()}}


def state_at(nodes, complete, logq, logt, sample_q, sample_t, split, cells=None):
    native = evaluate(nodes, complete, logq, logt, sample_q, sample_t, cell_indices=cells, add_radiation=split)
    return thermodynamic_partials(native, sample_q, sample_t)


def positivity(state, sample_q, sample_t):
    result = {}
    for field in ("P", "Pgas", "E", "S", "Cv", "Cp", "P_rho", "chiRho"):
        failed = (state[field] <= 0) | ~np.isfinite(state[field])
        result[field] = {"count": int(failed.sum()), "nonfinite_count": int((~np.isfinite(state[field])).sum()),
            "total": len(sample_q), "examples": [
            {"logQ": float(sample_q[index]), "logT": float(sample_t[index]),
             "value": float(state[field][index]) if np.isfinite(state[field][index]) else None}
            for index in np.flatnonzero(failed)[:8]]}
        if failed.any():
            result[field]["extrema"] = state_distribution(state[field], sample_q, sample_t)
    return result


def comparison(first, second, sample_q, sample_t, fields=VALUES):
    result = {}
    for field in fields:
        scale = abs(second[field])
        if field.endswith(("_x", "_t")):
            scale = np.maximum(scale, abs(second[field[:-2]]))
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            difference = (first[field] - second[field]) / np.maximum(scale, np.finfo(float).tiny)
        metric = distribution(difference)
        finite = np.isfinite(difference)
        if finite.any():
            index = int(np.argmax(np.where(finite, abs(difference), -1)))
            metric["worst"] = {"logQ": float(sample_q[index]), "logT": float(sample_t[index]), "signed_scaled_difference": float(difference[index]),
                               "candidate": float(first[field][index]), "reference": float(second[field][index])}
        result[field] = metric
    return result


def face_comparison(left, right, sample_q, sample_t):
    metrics = {}
    for field in VALUES + JACOBIANS:
        scale = np.maximum(abs(left[field]), abs(right[field]))
        if field in JACOBIANS:
            scale = np.maximum.reduce([scale, abs(left[field[:-2]]), abs(right[field[:-2]])])
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            gap = (left[field] - right[field]) / np.maximum(scale, np.finfo(float).tiny)
        metric = distribution(gap)
        finite = np.isfinite(gap)
        if finite.any():
            index = int(np.argmax(np.where(finite, abs(gap), -1)))
            metric["worst"] = {"logQ": float(sample_q[index]), "logT": float(sample_t[index]), "signed_gap": float(gap[index])}
        metrics[field] = metric
    return metrics


def edge_faces(nodes, complete, core_nodes, core_complete, logq, logt, rows, split, count):
    cells = valid_cells(complete) & valid_cells(material_nodes(rows, logq, logt))
    cells &= logt[None, :-1] < 3.36 + 1e-12
    generator = np.random.default_rng(20261028020)
    result = {}
    for axis, direction in ((0, "u"), (1, "tau")):
        adjacent = cells[:-1, :] & cells[1:, :] if axis == 0 else cells[:, :-1] & cells[:, 1:]
        density, temperature = np.nonzero(adjacent)
        choice = generator.choice(len(density), min(count, len(density)), replace=False)
        density, temperature = density[choice], temperature[choice]
        fraction = generator.uniform(.1, .9, len(choice))
        if axis == 0:
            sample_q = logq[density + 1]
            sample_t = logt[temperature] + fraction * (logt[1] - logt[0])
            right_cells = density + 1, temperature
        else:
            sample_q = logq[density] + fraction * (logq[1] - logq[0])
            sample_t = logt[temperature + 1]
            right_cells = density, temperature + 1
        left_cells = density, temperature
        left = state_at(nodes, complete, logq, logt, sample_q, sample_t, split, left_cells)
        right = state_at(nodes, complete, logq, logt, sample_q, sample_t, split, right_cells)
        overlap = valid_cells(core_complete)[left_cells] & valid_cells(core_complete)[right_cells]
        groups = {}
        for name, chosen in (("new_coverage", ~overlap), ("overlap", overlap)):
            group_left = {field: values[chosen] for field, values in left.items()}
            group_right = {field: values[chosen] for field, values in right.items()}
            groups[name] = {"count": int(chosen.sum()),
                "scaled_jumps": face_comparison(group_left, group_right, sample_q[chosen], sample_t[chosen])}
        result[direction] = {"count": len(sample_q), "overlap_core_faces": int(overlap.sum()),
            "scaled_jumps": face_comparison(left, right, sample_q, sample_t), "groups": groups,
            "coordinates": {"minimum_logT": float(sample_t.min()), "maximum_logT": float(sample_t.max())},
            "positivity_left": positivity(left, sample_q, sample_t), "positivity_right": positivity(right, sample_q, sample_t),
            "caveat": "shared nine jets retain C2 potential, not C3; new coverage does not remove response-Jacobian jumps"}
    return result


def main():
    workspace = Path(__file__).resolve().parents[5]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--split-radiation", action="store_true")
    parser.add_argument("--overlap-sample", type=int, default=5000)
    parser.add_argument("--faces", type=int, default=1000)
    arguments = parser.parse_args()
    if arguments.overlap_sample < 1 or arguments.faces < 1:
        parser.error("Sample counts must be positive")
    output = workspace / "freeeos_biquintic_progress/experiments/020_edge_control"
    if arguments.split_radiation:
        output /= "split_radiation"
    source_path = workspace / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
    archive_path = workspace / "freeeos_report/plots/offgrid_sample_data.npz"
    rejected_path = workspace / "freeeos_biquintic_progress/experiments/019_coverage_gap/rejected_queries.npz"
    coverage_path = rejected_path.with_name("results.json")
    core_path = Path(__file__).resolve().parents[1] / "freeeos_potential.py"
    fingerprints = {str(path): checksum(path) for path in (source_path, archive_path, rejected_path, coverage_path, core_path,
        core_path.with_name("validate_table.py"), core_path.with_name("active_faces.py"), Path(__file__))}
    inherited = json.loads(coverage_path.read_text())
    if fingerprints[str(source_path)] != inherited["provenance"]["source_sha256"] or fingerprints[str(archive_path)] != inherited["provenance"]["archive_sha256"]:
        raise ValueError("Archived coverage provenance differs")
    started = time.perf_counter()
    _, logq, logt, rows = read_table(source_path)
    low, valid, _ = source_constraints(rows, logq, logt, split_radiation=arguments.split_radiation)
    du, dtau = (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10
    core_nodes, core_complete, _ = construct_twists(low, valid, du, dtau)
    centered_nodes, centered_complete, _ = build_jets(low, valid, du, dtau, shifted=False)
    np.testing.assert_array_equal(centered_complete, core_complete)
    np.testing.assert_array_equal(centered_nodes[:, :, core_complete], core_nodes[:, :, core_complete])
    nodes, complete, metadata = build_jets(low, valid, du, dtau)
    if np.any(core_complete & ~complete):
        raise AssertionError("Shifted preference lost centered coverage")
    np.testing.assert_array_equal(nodes[:, :, core_complete], core_nodes[:, :, core_complete])
    with np.load(archive_path) as archive:
        sample_q, sample_t = archive["logQ"].ravel(), archive["logT"].ravel()
    density, temperature, inside = query_cells(logq, logt, sample_q, sample_t)
    baseline = inside & valid_cells(core_complete)[density, temperature]
    control = inside & valid_cells(complete)[density, temperature]
    clean_source_cells = inside & valid_cells(valid)[density, temperature]
    if np.any(control & ~clean_source_cells):
        raise AssertionError("Control accepted a non-native/invalid source cell")
    with np.load(rejected_path) as rejected:
        rejected_indices = rejected["archive_flat_indices"]
        probe = rejected["probe"]
        np.testing.assert_allclose(rejected["logQ"], sample_q[rejected_indices], rtol=0, atol=0)
        np.testing.assert_allclose(rejected["logT"], sample_t[rejected_indices], rtol=0, atol=0)
    np.testing.assert_array_equal(rejected_indices, np.flatnonzero(~baseline))
    recovered = control[rejected_indices]
    masks = classifications(probe[:, 2], probe[:, 5])
    masks["source_clean"] = clean_source_cells[rejected_indices]
    rejected_t = sample_t[rejected_indices]
    bands = {"all": np.ones(len(recovered), bool), "below_3.2": rejected_t < 3.2,
             "3.2_to_3.28": (rejected_t >= 3.2) & (rejected_t < 3.28), "atleast_3.28": rejected_t >= 3.28}
    results = {"purpose": "shifted/one-sided native quartic missing-twist coverage CONTROL, NOT runtime fallback/backend recommendation",
        "provenance": {"input_sha256": fingerprints, "command": sys.argv, "python": sys.version, "numpy": np.__version__, "inherited_fractions": "experiment019 default-MESA probe; not rerun or assigned to candidate runtime"},
        "method": {"degree": 4, "contiguous_donors": 5, "query_window_position_preference": [2, 1, 3, 0, 4],
            "routes": "same redundant 21/12 mean routes and two-stage 22 mean route as core; shift only when centered five valid donors unavailable",
            "source_floor": 3.2, "flags_density_finite_positivity_exclusions_unchanged": True, "split_radiation": arguments.split_radiation,
            "normalization": "value/holdout departures: abs(reference); Jacobian comparisons: max(abs(reference derivative), abs(reference value)); face gaps: max(abs(left), abs(right)), also including both values for Jacobians; x/t denote physical natural-log rho/T partials",
            "positivity_interpretation": "nonpositive/nonfinite counts are sampled diagnostics; E/S signs are printed-convention checks, not gauge-invariant stability conditions",
            "uncertainty": "one-sided first-derivative L1 weight norm can be ~7.1x centered; two-stage estimates amplify noise/truncation; not missing exact higher source derivatives"},
        "counts": {"archived_queries": len(sample_q), "core_accepted": int(baseline.sum()), "control_accepted": int(control.sum()),
                   "newly_accepted": int((control & ~baseline).sum()), "source_clean_rejected": int((clean_source_cells & ~baseline).sum()),
                   "source_clean_still_rejected": int((clean_source_cells & ~control).sum()), "core_complete_nodes": int(core_complete.sum()), "control_complete_nodes": int(complete.sum())},
        "minimum_complete_logT": {"core": float(logt[np.nonzero(core_complete)[1]].min()), "control": float(logt[np.nonzero(complete)[1]].min())},
        "stencil_support": reach_report(complete, metadata, logq, logt), "coverage_by_band": {}}
    for band, selected in bands.items():
        results["coverage_by_band"][band] = {category: {"previously_rejected": int((selected & category_mask).sum()),
            "recovered": int((selected & category_mask & recovered).sum()), "remaining_rejected": int((selected & category_mask & ~recovered).sum())}
            for category, category_mask in masks.items()}
    if control[sample_t < 3.2 - 1e-12].any():
        raise AssertionError("Control crossed excluded low-temperature floor")
    results["recovered_query_diagnostics"] = {}
    for category in ("all_successful", "active", "pure_FreeEOS", "blend"):
        selected = recovered & masks[category]
        query_indices = rejected_indices[selected]
        state = state_at(nodes, complete, logq, logt, sample_q[query_indices], sample_t[query_indices], arguments.split_radiation)
        results["recovered_query_diagnostics"][category] = {"count": len(query_indices), "positivity": positivity(state, sample_q[query_indices], sample_t[query_indices]),
            "values": {field: state_distribution(state[field], sample_q[query_indices], sample_t[query_indices]) for field in VALUES},
            "jacobians": {field: state_distribution(state[field], sample_q[query_indices], sample_t[query_indices]) for field in JACOBIANS},
            "value_caveat": "no corresponding accepted core prediction and no arbitrary-query source oracle; distributions are not errors"}
    overlap = np.flatnonzero(baseline & control)
    choice = np.random.default_rng(20261028020).choice(overlap, min(arguments.overlap_sample, len(overlap)), replace=False)
    core_state = state_at(core_nodes, core_complete, logq, logt, sample_q[choice], sample_t[choice], arguments.split_radiation)
    control_state = state_at(nodes, complete, logq, logt, sample_q[choice], sample_t[choice], arguments.split_radiation)
    for field in VALUES + JACOBIANS:
        np.testing.assert_array_equal(control_state[field], core_state[field])
    results["overlap"] = {"sample": len(choice), "all_values_jacobians_bitwise_equal": True,
                          "differences": comparison(control_state, core_state, sample_q[choice], sample_t[choice], VALUES + JACOBIANS)}
    coarse_q, coarse_t, coarse_rows = logq[::2], logt[::2], rows[::2, ::2]
    coarse_low, coarse_valid, _ = source_constraints(coarse_rows, coarse_q, coarse_t, split_radiation=arguments.split_radiation)
    coarse_nodes, coarse_complete, coarse_metadata = build_jets(coarse_low, coarse_valid, (coarse_q[1] - coarse_q[0]) * LN10, (coarse_t[1] - coarse_t[0]) * LN10)
    coarse_core_nodes, coarse_core_complete, _ = construct_twists(coarse_low, coarse_valid, (coarse_q[1] - coarse_q[0]) * LN10, (coarse_t[1] - coarse_t[0]) * LN10)
    held_density, held_temperature = np.nonzero(valid & material_nodes(rows, logq, logt) & (logt[None, :] >= 3.2 - 1e-12) & (logt[None, :] <= 3.6 + 1e-12))
    selected = (held_density % 2 == 1) & (held_temperature % 2 == 1)
    held_density, held_temperature = held_density[selected], held_temperature[selected]
    held_q, held_t = logq[held_density], logt[held_temperature]
    coarse_density, coarse_temperature, held_inside = query_cells(coarse_q, coarse_t, held_q, held_t)
    held_control = held_inside & valid_cells(coarse_complete)[coarse_density, coarse_temperature]
    held_core = held_inside & valid_cells(coarse_core_complete)[coarse_density, coarse_temperature]
    results["stride2_lowedge_holdouts"] = {"source_native_material_withheld": len(held_q), "accepted_control": int(held_control.sum()), "accepted_core": int(held_core.sum()),
        "stencil_support": reach_report(coarse_complete, coarse_metadata, coarse_q, coarse_t),
        "minimum_coarse_valid_logT": float(coarse_t[np.nonzero(coarse_valid)[1]].min()),
        "limitation": "globally even stride2 coarse-only training and same odd/odd held-out native rows, clean logT>=3.2; minimum coarse source floor3.22 excludes heldout3.2; archived fidelity not full-resolution source accuracy", "groups": {}}
    reference = table_state(rows[held_density, held_temperature])
    for label, chosen in (("new_coverage", held_control & ~held_core), ("overlap", held_control & held_core)):
        state = state_at(coarse_nodes, coarse_complete, coarse_q, coarse_t, held_q[chosen], held_t[chosen], arguments.split_radiation)
        group = {"count": int(chosen.sum()), "archived_relative_departure": comparison(state, {field: values[chosen] for field, values in reference.items()}, held_q[chosen], held_t[chosen], tuple(reference)),
                 "positivity": positivity(state, held_q[chosen], held_t[chosen]),
                 "minimum_logT": float(held_t[chosen].min()), "maximum_logT": float(held_t[chosen].max())}
        if label == "overlap":
            core_state = state_at(coarse_core_nodes, coarse_core_complete, coarse_q, coarse_t, held_q[chosen], held_t[chosen], arguments.split_radiation)
            group["core_archived_relative_departure"] = comparison(core_state, {field: values[chosen] for field, values in reference.items()}, held_q[chosen], held_t[chosen], tuple(reference))
        results["stride2_lowedge_holdouts"]["groups"][label] = group
    results["lowedge_material_faces"] = edge_faces(nodes, complete, core_nodes, core_complete, logq, logt, rows, arguments.split_radiation, arguments.faces)
    results["elapsed_seconds"] = time.perf_counter() - started
    for path, expected in fingerprints.items():
        if checksum(path) != expected:
            raise RuntimeError("Read-only input changed during the experiment")
    output.mkdir(parents=True, exist_ok=True)
    (output / "results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")
    stencil_arrays = {name + "_" + field: values for name, data in metadata["stencils"].items() for field, values in data.items()}
    np.savez_compressed(output / "stencils.npz", source_valid=valid, core_complete=core_complete, control_complete=complete,
                        twist_valid=metadata["twist_valid"], source_support_bounds=metadata["support"], **stencil_arrays)
    np.savez_compressed(output / "query_masks.npz", core_accepted=baseline, control_accepted=control,
                        source_cells_clean=clean_source_cells, archive_rejected_indices=rejected_indices, recovered_rejected=recovered)
    lines = ["# Native shifted-stencil coverage CONTROL", "", "NOT a recommended fallback or runtime backend. No low-temperature bug repair; donors/query cells remain native and clean with logT>=3.2.", "",
        f"Archived queries: {len(sample_q)}. Centered accepted {baseline.sum()}; control accepted {control.sum()}; newly accepted {(control & ~baseline).sum()}.",
        f"Complete-node floor: centered {results['minimum_complete_logT']['core']:.4g}, control {results['minimum_complete_logT']['control']:.4g}.", "",
        "| Band / category (MESA class where applicable) | Previously rejected | Recovered | Still rejected |", "|---|---:|---:|---:|"]
    for band, data in results["coverage_by_band"].items():
        for category in ("all_successful", "source_clean", "active", "pure_FreeEOS", "blend"):
            row = data[category]
            lines.append(f"| {band} / {category} | {row['previously_rejected']} | {row['recovered']} | {row['remaining_rejected']} |")
    lines.extend(["", f"Centered windows win whenever safe; optional shifted contiguous quartic windows never bridge invalid donors. Existing complete nine jets and {results['overlap']['sample']:,} sampled overlap values/Jacobians match core bitwise.",
        "", f"Maximum two-stage donor reach: q={results['stencil_support']['max_density_node_reach']} nodes, T={results['stencil_support']['max_temperature_node_reach']} nodes; minimum donor logT={results['stencil_support']['minimum_donor_logT']}.",
        "", "## Quality and uncertainty", ""])
    diagnostics = results["recovered_query_diagnostics"]
    for category in ("active", "all_successful"):
        failed = ", ".join(f"{field}={data['count']}" for field, data in diagnostics[category]["positivity"].items() if data["count"])
        lines.append(f"Recovered {category} ({diagnostics[category]['count']} queries): nonpositive/nonfinite counts {failed or 'all zero for tested fields'}. E/S signs are printed-convention checks, not gauge-invariant stability conditions.")
    worst_cv = diagnostics["all_successful"]["values"]["Cv"]["signed_minimum"]
    lines.extend(["", f"Minimum recovered Cv={worst_cv['value']:.6g} at logQ={worst_cv['logQ']:.12g}, logT={worst_cv['logT']:.12g}. Broad archived-query and inherited active classifications are reported separately, not interchangeable; own_use was not computed in this main experiment.",
        "", "### Stride2 native-row holdout", "", "Coarse-only globally even nodes, odd/odd withheld material rows: source floor3.22 (no donor below3.2, no extrapolation to3.2). Printed-row comparisons, NOT source accuracy. New and overlap populations differ; their tails cannot isolate stencil causation.",
        "", "| Population / field | Nodes | Median relative departure | p99 | Max |", "|---|---:|---:|---:|---:|"])
    for name, group in results["stride2_lowedge_holdouts"]["groups"].items():
        for field in ("Pgas", "Cv", "Cp", "chiRho", "chiT", "gamma3"):
            metric = group["archived_relative_departure"][field]
            lines.append(f"| {name} / {field} | {group['count']} | {metric['median']:.6g} | {metric['p99']:.6g} | {metric['max']:.6g} |")
    worst_holdout = results["stride2_lowedge_holdouts"]["groups"]["new_coverage"]["archived_relative_departure"]["Cv"]["worst"]
    holdout_failures = sum(metric["count"] for group in results["stride2_lowedge_holdouts"]["groups"].values() for metric in group["positivity"].values())
    lines.extend(["", f"Worst new holdout Cv: logQ={worst_holdout['logQ']:.12g}, logT={worst_holdout['logT']:.12g}; candidate={worst_holdout['candidate']:.6g}, printed={worst_holdout['reference']:.6g}. Tested holdout positivity violation count={holdout_failures} (summed across fields). Gamma3 is the printed Gamma3, not Gamma3-minus1.",
        "", "### Exact low-edge shared faces", "", f"Up to {arguments.faces:,} sampled faces/orientation, material cells through logT~3.38 (selected cell starts<=3.36); face samples are NOT fraction-classified active runtime states. Forced adjacent cells at identical coordinates, no extrapolated finite steps. Physical x/t partials use natural-log rho/T. Jacobian gaps normalize by max(abs(left derivative), abs(right derivative), abs(left value), abs(right value)).", "",
        "| Orientation / population / partial | Faces | Median scaled jump | p99 | Max |", "|---|---:|---:|---:|---:|"])
    for direction, data in results["lowedge_material_faces"].items():
        for name, group in data["groups"].items():
            for field in ("Cv_t", "chiT_t", "Cp_t"):
                metric = group["scaled_jumps"][field]
                lines.append(f"| {direction} / {name} / {field} | {group['count']} | {metric['median']:.6g} | {metric['p99']:.6g} | {metric['max']:.6g} |")
    lines.extend(["", "Value continuity is consistent with C2 to roundoff; response Jacobian tails remain. The worst u/tau Cv_t coordinates and all value/Jacobian distributions are in results.json. Recovered off-node values have no source oracle; raw distributions are magnitudes with signed extrema/coordinates, not errors.",
        "", "A quartic one-sided first derivative can amplify input noise about 7.1 times more than centered differences; two stages compound truncation/rounding. Estimated third/fourth twists remain uncertain. Shared jets are C2, not C3; response Jacobian jumps need not improve.",
        "", "Coverage below logT3.2 stays deliberately rejected even when archived MESA has a positive FreeEOS fraction. Source flags/density mismatches/short native runs remain exclusions. Fractions are inherited classifications, not a rerun or a claim for a candidate runtime path.",
        "", "No gauge/source accuracy, solver robustness, boundary blending or runtime fallback safety is established. Exact window starts/query positions and recursive support bounds are archived in stencils.npz; query decisions are in query_masks.npz."])
    (output / "SUMMARY.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({"counts": results["counts"], "active_coverage": {band: data["active"] for band, data in results["coverage_by_band"].items()}}, indent=2))


if __name__ == "__main__":
    main()
