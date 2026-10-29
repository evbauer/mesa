"""Reuse archived fits at three lead-measured active faces; no runtime rerun."""

import json
from pathlib import Path

import numpy as np

from validate_multislice import CONFIGURATIONS
from freeeos_potential import LN10, checksum, construct_twists, evaluate, query_cells, read_table, source_constraints, thermodynamic_partials
from smooth_patch import evaluate_fit


def face_scale(left, right, field):
    return max(1., abs(float(left[field][0])), abs(float(right[field][0])),
               abs(float(left[field[:-2]][0])), abs(float(right[field[:-2]][0])))


def main():
    workspace = Path(__file__).resolve().parents[5]
    output = workspace / "freeeos_biquintic_progress/experiments/015_smooth_validation"
    record_path = workspace / "freeeos_biquintic_progress/experiments/014_active_domain/face_results.json"
    lead = json.loads(record_path.read_text())
    results = json.loads((output / "results.json").read_text())
    patch = results["patches"]["active_chiT"]
    source_path = workspace / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
    if checksum(source_path) != lead["provenance"]["source_sha256"]:
        raise ValueError("Lead face record and source differ")
    _, logq, logt, rows = read_table(source_path)
    low, valid, _ = source_constraints(rows, logq, logt)
    nodes, complete, _ = construct_twists(low, valid, (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10)
    archive_path = output / "active_chiT/fits_and_queries.npz"
    record = {"scope": "three lead runtime-selected faces; archived C4 fits only, no runtime or active-fraction rerun; total-potential Python biquintic control",
              "provenance": {"lead_faces_sha256": checksum(record_path), "fit_archive_sha256": checksum(archive_path),
                             "source_sha256": checksum(source_path), "script_sha256": checksum(__file__)},
              "scaling": "lead convention max(1,abs(both responses),abs(both physical natural-log response derivatives))",
              "active_distribution_context": {axis: lead["faces"][axis]["active"]["jumps"]["Cv_t"]["p99"] for axis in ("u", "tau")},
              "faces": {}}
    with np.load(archive_path) as archived:
        patch_q, patch_t = archived["patch_logQ"], archived["patch_logT"]
        for direction, category in (("u", "pure_FreeEOS"), ("u", "blend"), ("tau", "pure_FreeEOS")):
            worst = lead["faces"][direction][category]["jumps"]["Cv_t"]["worst"]
            sample_q, sample_t = np.array([worst["logQ"]]), np.array([worst["logT"]])
            if not (patch_q[0] < sample_q[0] < patch_q[-1] and patch_t[0] < sample_t[0] < patch_t[-1]):
                raise AssertionError("Selected active patch misses a required lead face")
            density_index, temperature_index, _ = query_cells(logq, logt, sample_q, sample_t)
            if direction == "u":
                boundary = int(np.argmin(abs(logq - sample_q[0])))
                sample_q[0] = logq[boundary]
                left_indices, right_indices = (np.array([boundary - 1]), temperature_index), (np.array([boundary]), temperature_index)
            else:
                boundary = int(np.argmin(abs(logt - sample_t[0])))
                sample_t[0] = logt[boundary]
                left_indices, right_indices = (density_index, np.array([boundary - 1])), (density_index, np.array([boundary]))
            left = thermodynamic_partials(evaluate(nodes, complete, logq, logt, sample_q, sample_t, cell_indices=left_indices), sample_q, sample_t)
            right = thermodynamic_partials(evaluate(nodes, complete, logq, logt, sample_q, sample_t, cell_indices=right_indices), sample_q, sample_t)
            face = {"lead_selected_face": worst, "covered_by_selected_patch": True,
                    "python_biquintic_scaled_jumps": {field: float((left[field][0] - right[field][0]) / face_scale(left, right, field)) for field in ("Cv_t", "chiT_t")},
                    "candidates": {}}
            for candidate_name, _, _ in CONFIGURATIONS:
                face["candidates"][candidate_name] = {}
                for resolution, grid_q, grid_t in (("full", patch_q, patch_t), ("coarse", patch_q[::2], patch_t[::2])):
                    coefficients = archived[candidate_name + "_" + resolution + "_coefficients"]
                    state = thermodynamic_partials(evaluate_fit(coefficients, grid_q * LN10, grid_t * LN10, sample_q * LN10, sample_t * LN10), sample_q, sample_t)
                    steps, gaps = {}, []
                    for fraction in (1e-3, 1e-4, 1e-5):
                        step = fraction * (grid_q[1] - grid_q[0] if direction == "u" else grid_t[1] - grid_t[0]) * LN10
                        delta_q, delta_t = (step / LN10, 0.) if direction == "u" else (0., step / LN10)
                        first = thermodynamic_partials(evaluate_fit(coefficients, grid_q * LN10, grid_t * LN10, (sample_q - delta_q) * LN10, (sample_t - delta_t) * LN10), sample_q - delta_q, sample_t - delta_t)
                        second = thermodynamic_partials(evaluate_fit(coefficients, grid_q * LN10, grid_t * LN10, (sample_q + delta_q) * LN10, (sample_t + delta_t) * LN10), sample_q + delta_q, sample_t + delta_t)
                        scaled = {field: float((first[field][0] - second[field][0]) / face_scale(first, second, field)) for field in ("Cv_t", "chiT_t")}
                        steps[str(fraction)] = scaled
                        gaps.append(scaled)
                    face["candidates"][candidate_name][resolution] = {"finite_one_sided_scaled_gaps": steps,
                        "linear_zero_step_intercept_scaled": {field: (10 * gaps[-1][field] - gaps[-2][field]) / 9 for field in gaps[-1]},
                        "at_face": {field: float(state[field][0]) for field in ("Cv", "chiT", "Cv_t", "chiT_t")},
                        "value_departure_vs_continuous_biquintic_value": {field: float(state[field][0] / left[field][0] - 1) for field in ("Cv", "chiT")}}
            record["faces"][direction + "_" + category] = face
    (output / "active_exact_faces.json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    lines = ["# Multiregion interpretation and active-face follow-up", "", "Table-only evidence: a smoother potential is not necessarily a more faithful EOS.", "",
             "## Holdout tradeoffs", "", "All 32 full/coarse fits solved; no sampled nodal/off-node/held-out positivity violations. Each patch has 5,000 strictly off-node queries; 746 native odd/odd holdouts total. Five analytic tests pass.", "",
             "| Region | Lowest C4 coarse Cv p99 | Weight | Global coarse biquintic Cv p99 |", "|---|---:|---|---:|"]
    for name, data in results["patches"].items():
        values = [(candidate["coarse"]["archived_withheld_departure_NOT_accuracy"]["all"]["relative_departure"]["Cv"]["p99"], label)
                  for label, candidate in data["candidates"].items() if candidate["coarse"]["status"] == "solved"]
        best, label = min(values)
        baseline = data["coarse_biquintic_global_donor_control"]["archived_withheld_departure"]["all"]["relative_departure"]["Cv"]["p99"]
        lines.append(f"| {name} | {best:.6g} | {label} | {baseline:.6g} |")
    lines.extend(["", "Increasing response weight improves full-resolution nodal Cv reconciliation but can worsen coarse holdouts, especially in H/He ionization. Dense-tail and active-chiT patches favor thermodynamic relative weights over the native-weight control; no single weight wins all regions or all fields.",
        "", "Boundary holdouts matter: patch-only coarse biquintic has only 18–36 complete interior holdouts, while global clean donors permit all 170–204. Do not compare different query-set p99 values as a donor-method effect; both controls on the identical common subset are archived.",
        "", "## Exact active faces", "", "The active patch q=[0.35,1.37], logT=[6.5,6.98] covers the pure-FreeEOS u face, blended u face, and pure-FreeEOS tau face from experiment 014. Runtime fractions are inherited observations, not recomputed for or applied to the fitted potential.", "",
        "| Face | Lead fraction | Biquintic Cv_t / chiT_t jumps | thermo_0.1 full Cv_t / chiT_t zero-step intercept |", "|---|---:|---|---|"])
    for name, face in record["faces"].items():
        baseline = face["python_biquintic_scaled_jumps"]
        fitted = face["candidates"]["thermo_0.1"]["full"]["linear_zero_step_intercept_scaled"]
        lines.append(f"| {name} | {face['lead_selected_face']['fraction_FreeEOS']:.6g} | {baseline['Cv_t']:.6g} / {baseline['chiT_t']:.6g} | {fitted['Cv_t']:.3g} / {fitted['chiT_t']:.3g} |")
    lines.extend(["", "These finite-step intercepts approach zero, consistent with the C4 basis; they are not exact arithmetic continuity certificates or evidence of correct derivative values. The fitted Cv/chiT values and derivative levels shift, and their departures at these faces are recorded.",
        "", "Active p99 Cv_t jumps from experiment 014 are u=" + str(record["active_distribution_context"]["u"]) + ", tau=" + str(record["active_distribution_context"]["tau"]) + ". Dense all-own-use extremes that are Skye-only must not be substituted for this runtime-active population.",
        "", "## Actionable experimental gates", "", "- Keep weight sweeps and boundary/interior holdouts together: exact nodal agreement is deliberately relaxed by the C4 soft fit.",
        "- Use thermo_0.01/0.1 as candidates for further table-only active-region comparisons, not a runtime/default choice. Check chiT/Cp/Gamma and worst rows, not only Cv p99.",
        "- Fourth-potential-derivative relative intercepts reach a roundoff-sensitive floor in small ionization derivatives; response-Jacobian intercepts are much smaller. Report absolute/relative conditioning, not an unsupported uniform machine-precision C4 claim.",
        "- Patch ends are not blended or globally C4; gauges, source accuracy, mixture response, actual active fractions and Newton inversion remain separate lead gates.",
        "", "No source/runtime code edits, source builds, commits or pushes. Numerical fingerprints and the exact-face sampling/scaling are in results.json and active_exact_faces.json."])
    (output / "INTERPRETATION.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({name: face["candidates"]["thermo_0.1"]["full"]["linear_zero_step_intercept_scaled"] for name, face in record["faces"].items()}, indent=2))


if __name__ == "__main__":
    main()
