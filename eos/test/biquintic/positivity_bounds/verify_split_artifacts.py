"""Independent split artifact recount, root reconstruction and core witness replay."""

import sys

sys.dont_write_bytecode = True

import json
from pathlib import Path

import numpy as np

from bounds import DTYPE, FIELDS, STATUS, bernstein_bounds, cell_power, normalize_closed_cell, sign_powers
from bound_split import OUTPUT, TOTAL_OUTPUT, radiation_cv_numerator, radiation_guard
from run_bounds import BIQUINTIC, TABLE, WORKSPACE, collect_jets
from freeeos_potential import LN10, checksum, construct_twists, evaluate, read_table, source_constraints, thermodynamics, valid_cells


def main():
    result = json.loads((OUTPUT / "results.json").read_text())
    archive = np.load(OUTPUT / "cell_bounds.npz", allow_pickle=False)
    baseline = np.load(TOTAL_OUTPUT / "cell_bounds.npz", allow_pickle=False)
    witnesses = json.loads((OUTPUT / "negative_witnesses.json").read_text())
    unresolved = json.loads((OUTPUT / "unresolved_cells.json").read_text())
    assert result["calibration_passed"]
    assert result["floating_proof_claim"] is False
    assert result["source_unchanged"] and result["baseline024_unchanged"]
    for name, digest in result["artifacts_sha256"].items():
        assert checksum(OUTPUT / name) == digest, name
    for name, digest in result["baseline024_sha256"].items():
        assert checksum(TOTAL_OUTPUT / name) == digest, name
    for name, digest in result["source_sha256"].items():
        assert checksum(WORKSPACE / name) == digest, name
    for name, digest in result["sidecar_sha256"].items():
        assert checksum(Path(__file__).with_name(name)) == digest, name
    for name in ("density_cell_index", "temperature_cell_index", "full_one_branch_own_use", "material_corner_screen", "maximum_corner_Prad_over_Pgas"):
        np.testing.assert_array_equal(archive[name], baseline[name])
    density_indices, temperature_indices = archive["density_cell_index"], archive["temperature_cell_index"]
    count = len(density_indices)
    cells = list(zip(density_indices.tolist(), temperature_indices.tolist()))
    assert count == len(set(cells)) == result["domain"]["complete_cells"]
    full_own, material = archive["full_one_branch_own_use"], archive["material_corner_screen"]
    masks = {"all_complete": np.ones(count, dtype=bool), "full_one_branch_own_use": full_own,
             "material_full_one_branch_own_use": material & full_own,
             "material_not_full_one_branch_own_use": material & ~full_own,
             "not_full_one_branch_own_use": ~full_own}
    for field in FIELDS:
        status, representation_status = archive[f"{field}__status"], archive[f"{field}__representation_status"]
        assert np.all(archive[f"{field}__final_lower_bound"][status == 1] > 0)
        assert np.all(archive[f"{field}__representation_lower_bound"][representation_status == 1] > 0)
        assert np.all(archive[f"{field}__visited_nodes"] <= result["budgets"]["max_nodes_per_mode_cell_field"])
        assert np.all(archive[f"{field}__representation_visited_nodes"] <= result["budgets"]["max_nodes_per_mode_cell_field"])
        assert np.all(archive[f"{field}__maximum_depth"] <= result["budgets"]["max_depth"])
        assert np.all(archive[f"{field}__unresolved_leaves"][status == 1] == 0)
        for entries, expected_status in ((witnesses[field], STATUS["negative_witness"]), (unresolved[field], STATUS["unresolved"])):
            entry_cells = [tuple(entry["cell"]) for entry in entries]
            assert len(entry_cells) == len(set(entry_cells))
            assert set(entry_cells) == {cells[index] for index in np.flatnonzero(status == expected_status)}
        for name, mask in masks.items():
            summary = result["counts"][name]
            assert summary["cell_count"] == int(mask.sum())
            field_summary = summary["fields"][field]
            for label, status_value in STATUS.items():
                assert field_summary[label] == int(np.count_nonzero(mask & (status == status_value)))
                assert field_summary["representation"][label] == int(np.count_nonzero(mask & (representation_status == status_value)))
            assert field_summary["restoration_guard_only_unresolved"] == int(np.count_nonzero(mask & (status == 0) & (representation_status == 1)))
            assert field_summary["total024_unresolved_now_positive"] == int(np.count_nonzero(mask & (baseline[f"{field}__status"] == 0) & (status == 1)))
            assert field_summary["total024_unresolved_now_representation_positive"] == int(np.count_nonzero(mask & (baseline[f"{field}__status"] == 0) & (representation_status == 1)))
            assert field_summary["root_positive"] == int(np.count_nonzero(mask & (archive[f"{field}__root_lower_bound"] > 0)))
    for name, mask in masks.items():
        for prefix, label in (("status", "both_restored_core_guarded_positive"), ("representation_status", "both_representation_positive")):
            positive = mask & (archive[f"Cv_numerator__{prefix}"] == 1) & (archive[f"P_rho__{prefix}"] == 1)
            assert result["counts"][name][label] == int(positive.sum())
    _, logq, logt, rows = read_table(TABLE)
    low, valid, _ = source_constraints(rows, logq, logt, split_radiation=True)
    du, dtau = (logq[1]-logq[0])*LN10, (logt[1]-logt[0])*LN10
    nodes, complete, _ = construct_twists(low, valid, du, dtau, radius=2, route="mean")
    expected_density, expected_temperature = np.nonzero(valid_cells(complete))
    np.testing.assert_array_equal(density_indices, expected_density)
    np.testing.assert_array_equal(temperature_indices, expected_temperature)
    selected = np.linspace(0, count-1, min(64, count), dtype=int)
    density_batch, temperature_batch = density_indices[selected], temperature_indices[selected]
    power, magnitude = cell_power(collect_jets(nodes, density_batch, temperature_batch), du, dtau)
    powers, scales = sign_powers(power, magnitude, du, dtau)
    density_stretch = (logq[density_batch+1].astype(DTYPE)-logq[density_batch].astype(DTYPE))/DTYPE(logq[1]-logq[0])
    temperature_stretch = (logt[temperature_batch+1].astype(DTYPE)-logt[temperature_batch].astype(DTYPE))/DTYPE(logt[1]-logt[0])
    radiation_min = radiation_cv_numerator(logq[density_batch+1], logt[temperature_batch])
    radiation_max = radiation_cv_numerator(logq[density_batch], logt[temperature_batch+1])
    for field in FIELDS:
        normalized_power = normalize_closed_cell(powers[field], density_stretch, temperature_stretch)
        normalized_scale = normalize_closed_cell(scales[field], density_stretch, temperature_stretch)
        coefficients, polynomial_margin, *_ = bernstein_bounds(normalized_power, normalized_scale)
        np.testing.assert_array_equal(coefficients.min(axis=(1, 2)), archive[f"{field}__root_polynomial_minimum"][selected])
        np.testing.assert_array_equal(polynomial_margin, archive[f"{field}__polynomial_margin"][selected])
        analytic_lower = radiation_min if field == "Cv_numerator" else DTYPE(0)
        analytic_guard = radiation_guard(field, radiation_max, False)
        restore_guard = radiation_guard(field, radiation_max, True)-analytic_guard
        root_lower = coefficients.min(axis=(1, 2))+analytic_lower-polynomial_margin-analytic_guard-restore_guard
        np.testing.assert_array_equal(root_lower, archive[f"{field}__root_lower_bound"][selected])
        entries = witnesses[field]
        if entries:
            sample_q = np.array([entry["witness"]["logQ"] for entry in entries])
            sample_t = np.array([entry["witness"]["logT"] for entry in entries])
            witness_density = np.array([entry["cell"][0] for entry in entries])
            witness_temperature = np.array([entry["cell"][1] for entry in entries])
            native = evaluate(nodes, complete, logq, logt, sample_q, sample_t, add_radiation=True,
                              cell_indices=(witness_density, witness_temperature))
            state = thermodynamics(native, sample_q, sample_t)
            values = state["Cv" if field == "Cv_numerator" else "P_rho"]
            assert np.all(np.isfinite(values) & (values < 0))
            np.testing.assert_array_equal(values, [entry["witness"]["core_value"] for entry in entries])
            for entry in entries:
                assert entry["witness"]["core_add_radiation"] is True
                assert not entry["full_one_branch_own_use"]
    print(json.dumps({"passed": True, "recounted_cells": count, "recomputed_root_cells": len(selected),
                      "replayed_add_radiation_true_witnesses": {field: len(entries) for field, entries in witnesses.items()},
                      "unique_negative_witness_cells": len({tuple(entry["cell"]) for entries in witnesses.values() for entry in entries}),
                      "baseline024_and_protected_source_hashes_unchanged": True}, indent=2))


if __name__ == "__main__":
    main()
