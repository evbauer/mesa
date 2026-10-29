"""Recount artifacts and independently re-evaluate all witnesses plus 64 roots."""

import sys

sys.dont_write_bytecode = True

import json

import numpy as np

from bounds import (
    DTYPE, FIELDS, STATUS, bernstein_bounds, cell_power, domain_masks,
    normalize_closed_cell, sign_powers,
)
from run_bounds import (
    BIQUINTIC, OUTPUT, TABLE, WORKSPACE, collect_jets,
)
from freeeos_potential import (
    CRAD, LN10, checksum, construct_twists, evaluate, read_table,
    source_constraints, thermodynamics, valid_cells,
)


def main():
    result = json.loads((OUTPUT / "results.json").read_text())
    archive = np.load(OUTPUT / "cell_bounds.npz", allow_pickle=False)
    witnesses = json.loads((OUTPUT / "negative_witnesses.json").read_text())
    unresolved = json.loads((OUTPUT / "unresolved_cells.json").read_text())
    for name, digest in result["artifacts_sha256"].items():
        assert checksum(OUTPUT / name) == digest, name
    for name, digest in {**result["source_sha256"], **result["sidecar_sha256"]}.items():
        assert checksum(WORKSPACE / name) == digest, name
    assert result["calibration_passed"]
    assert result["floating_proof_claim"] is False
    assert result["source_unchanged"]
    assert result["Hermite_reproduction"]["passed"]
    assert result["Hermite_reproduction"]["jet_count"] == 36*result["Hermite_reproduction"]["cell_count"]
    density_indices, temperature_indices = archive["density_cell_index"], archive["temperature_cell_index"]
    count = len(density_indices)
    assert count == result["domain"]["complete_cells"]
    cells = list(zip(density_indices.tolist(), temperature_indices.tolist()))
    assert len(set(cells)) == count
    full_own, material = archive["full_one_branch_own_use"], archive["material_corner_screen"]
    subsets = {"all_complete": np.ones(count, dtype=bool), "full_one_branch_own_use": full_own,
               "material_full_one_branch_own_use": material & full_own,
               "material_not_full_one_branch_own_use": material & ~full_own,
               "not_full_one_branch_own_use": ~full_own}
    index_lookup = {cell: index for index, cell in enumerate(cells)}
    for field in FIELDS:
        status = archive[f"{field}__status"]
        root_minimum = archive[f"{field}__root_minimum"]
        margin = archive[f"{field}__roundoff_margin"]
        positive = status == STATUS["positive"]
        assert np.all(archive[f"{field}__final_lower_bound"][positive] > 0)
        assert np.all(status[root_minimum > margin] == STATUS["positive"])
        assert np.all(archive[f"{field}__unresolved_leaves"][positive] == 0)
        assert np.all(archive[f"{field}__visited_nodes"] >= 1)
        assert np.all(archive[f"{field}__visited_nodes"] <= result["budgets"]["max_nodes_per_cell_field"])
        assert np.all(archive[f"{field}__maximum_depth"] <= result["budgets"]["max_depth"])
        for filename, entries, expected_status in (("negative_witnesses", witnesses[field], STATUS["negative_witness"]),
                                                  ("unresolved_cells", unresolved[field], STATUS["unresolved"])):
            entry_cells = [tuple(entry["cell"]) for entry in entries]
            assert len(entry_cells) == len(set(entry_cells)), filename
            assert set(entry_cells) == {cells[index] for index in np.flatnonzero(status == expected_status)}, filename
            for entry in entries:
                index = index_lookup[tuple(entry["cell"])]
                assert entry["full_one_branch_own_use"] == bool(full_own[index])
                assert entry["material_corner_screen"] == bool(material[index])
        for name, mask in subsets.items():
            summary = result["counts"][name]
            assert summary["cell_count"] == int(mask.sum())
            field_summary = summary["fields"][field]
            for label, status_value in STATUS.items():
                assert field_summary[label] == int(np.count_nonzero(mask & (status == status_value)))
            assert sum(field_summary[label] for label in STATUS) == int(mask.sum())
            root_positive = int(np.count_nonzero(mask & (root_minimum > margin)))
            assert field_summary["root_positive"] == root_positive
            assert field_summary["adaptive_positive"] == field_summary["positive"]-root_positive
            for metric in ("visited_nodes", "unresolved_leaves"):
                assert field_summary[metric] == int(archive[f"{field}__{metric}"][mask].sum())
    for name, mask in subsets.items():
        positive = mask & (archive["Cv_numerator__status"] == STATUS["positive"]) & (archive["P_rho__status"] == STATUS["positive"])
        assert result["counts"][name]["both_numerically_certified_positive"] == int(positive.sum())
    _, logq, logt, rows = read_table(TABLE)
    low, valid, state = source_constraints(rows, logq, logt, split_radiation=False)
    du, dtau = (logq[1]-logq[0])*LN10, (logt[1]-logt[0])*LN10
    nodes, complete, _ = construct_twists(low, valid, du, dtau, radius=2, route="mean")
    expected_density, expected_temperature = np.nonzero(valid_cells(complete))
    np.testing.assert_array_equal(density_indices, expected_density)
    np.testing.assert_array_equal(temperature_indices, expected_temperature)
    ratio = CRAD * state["T"]**4 / (3*state["Pgas"])
    expected_own, expected_material, _, ratio_max = domain_masks(logq, logt, ratio)
    np.testing.assert_array_equal(full_own, expected_own[density_indices, temperature_indices])
    np.testing.assert_array_equal(material, expected_material[density_indices, temperature_indices])
    np.testing.assert_array_equal(archive["maximum_corner_Prad_over_Pgas"], ratio_max[density_indices, temperature_indices])
    selected = np.linspace(0, count-1, min(64, count), dtype=int)
    density_batch, temperature_batch = density_indices[selected], temperature_indices[selected]
    power, magnitude = cell_power(collect_jets(nodes, density_batch, temperature_batch), du, dtau)
    powers, scales = sign_powers(power, magnitude, du, dtau)
    density_stretch = (logq[density_batch+1].astype(DTYPE)-logq[density_batch].astype(DTYPE))/DTYPE(logq[1]-logq[0])
    temperature_stretch = (logt[temperature_batch+1].astype(DTYPE)-logt[temperature_batch].astype(DTYPE))/DTYPE(logt[1]-logt[0])
    for field in FIELDS:
        normalized_power = normalize_closed_cell(powers[field], density_stretch, temperature_stretch)
        normalized_scale = normalize_closed_cell(scales[field], density_stretch, temperature_stretch)
        coefficients, margin, *_ = bernstein_bounds(normalized_power, normalized_scale)
        np.testing.assert_array_equal(coefficients.min(axis=(1, 2)), archive[f"{field}__root_minimum"][selected])
        np.testing.assert_array_equal(coefficients.max(axis=(1, 2)), archive[f"{field}__root_maximum"][selected])
        np.testing.assert_array_equal(margin, archive[f"{field}__roundoff_margin"][selected])
        entries = witnesses[field]
        if entries:
            sample_q = np.array([entry["witness"]["logQ"] for entry in entries])
            sample_t = np.array([entry["witness"]["logT"] for entry in entries])
            witness_density = np.array([entry["cell"][0] for entry in entries])
            witness_temperature = np.array([entry["cell"][1] for entry in entries])
            native = evaluate(nodes, complete, logq, logt, sample_q, sample_t,
                              cell_indices=(witness_density, witness_temperature))
            state = thermodynamics(native, sample_q, sample_t)
            values = state["Cv" if field == "Cv_numerator" else "P_rho"]
            assert np.all(np.isfinite(values) & (values < 0))
            np.testing.assert_array_equal(values, [entry["witness"]["core_value"] for entry in entries])
            for entry, density_index, temperature_index, sample_density, sample_temperature in zip(
                    entries, witness_density, witness_temperature, sample_q, sample_t):
                assert logq[density_index] <= sample_density <= logq[density_index+1]
                assert logt[temperature_index] <= sample_temperature <= logt[temperature_index+1]
                assert entry["witness"]["actual_density_phase"] >= 0
                assert entry["witness"]["actual_density_phase"] <= 1
                assert entry["witness"]["actual_temperature_phase"] >= 0
                assert entry["witness"]["actual_temperature_phase"] <= 1
    print(json.dumps({"passed": True, "recounted_cells": count,
                      "recomputed_root_cells": len(selected),
                      "core_reverified_negative_witnesses": {field: len(entries) for field, entries in witnesses.items()},
                      "unique_negative_witness_cells": len({tuple(entry["cell"]) for entries in witnesses.values() for entry in entries}),
                      "all_source_core_sidecar_and_data_hashes_match": True,
                      "core_source_path": str(BIQUINTIC / "freeeos_potential.py")}, indent=2))


if __name__ == "__main__":
    main()
