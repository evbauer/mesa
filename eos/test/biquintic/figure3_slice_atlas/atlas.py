"""Rerank archived 060 pure FreeEOS faces and generate exact slice atlas data."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
CODE = HERE.parent
ROOT = HERE.parents[4]
SOURCE = ROOT / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
CACHE = ROOT / "mesa/data/eosDT_data/cache/mesa-FreeEOS_02z70x.bin"
OLD = ROOT / "freeeos_biquintic_progress/experiments/060_pure_freeeos_faces"
NEW = ROOT / "freeeos_biquintic_progress/experiments/071_figure3_slice_atlas"
sys.path.insert(0, str(CODE))
sys.path.insert(0, str(ROOT / "freeeos_report"))
import audit_freeeos as baseline
from freeeos_potential import (CRAD, LN10, checksum, construct_twists, evaluate, query_cells,
                               read_table, source_constraints, thermodynamic_partials, valid_cells)
from active_slice_extended.run_extended import legacy_state, potential_state


def load_construction():
    _, logq, logt, rows = read_table(SOURCE)
    low, valid, _ = source_constraints(rows, logq, logt, split_radiation=False, minimum_logT=3.2)
    nodes, complete, _ = construct_twists(low, valid, np.diff(logq)[0] * LN10,
                                          np.diff(logt)[0] * LN10, radius=2, route="mean")
    return logq, logt, nodes, complete


def rerank():
    NEW.mkdir(parents=True, exist_ok=True)
    old_report = json.loads((OLD / "results.json").read_text())
    old_presentation = json.loads((OLD / "figure_presentation.json").read_text())
    for name, expected in ((SOURCE, old_report["provenance"]["source_sha256"]),
                           (OLD / "face_manifest.npz", old_report["provenance"]["manifest_sha256"]),
                           (OLD / "face_results.npz", old_presentation["results_sha256"])):
        if checksum(name) != expected:
            raise ValueError(f"Source fingerprint mismatch: {name}")
    logq, logt, nodes, complete = load_construction()
    manifest = np.load(OLD / "face_manifest.npz")
    old = np.load(OLD / "face_results.npz")
    output = {}
    checks = {}
    for direction, axis in (("u", 0), ("tau", 1)):
        sample_q, sample_t = manifest[f"{direction}_logQ"], manifest[f"{direction}_logT"]
        density_index = manifest[f"{direction}_density_index"]
        temperature_index = manifest[f"{direction}_temperature_index"]
        selected = old[f"{direction}_pure"] & (old[f"{direction}_support_code"] == 0)
        positions = np.flatnonzero(selected)
        left_cv = np.full(len(sample_q), np.nan)
        right_cv = np.full(len(sample_q), np.nan)
        left_cvt = np.full(len(sample_q), np.nan)
        right_cvt = np.full(len(sample_q), np.nan)
        jump = np.full(len(sample_q), np.nan)
        for batch in np.array_split(positions, max(1, (len(positions) + 4095) // 4096)):
            first = density_index[batch], temperature_index[batch]
            second = density_index[batch] + (axis == 0), temperature_index[batch] + (axis == 1)
            left = thermodynamic_partials(evaluate(nodes, complete, logq, logt, sample_q[batch], sample_t[batch],
                                                   cell_indices=first), sample_q[batch], sample_t[batch])
            right = thermodynamic_partials(evaluate(nodes, complete, logq, logt, sample_q[batch], sample_t[batch],
                                                    cell_indices=second), sample_q[batch], sample_t[batch])
            for target, values in ((left_cv, left["Cv"]), (right_cv, right["Cv"]),
                                   (left_cvt, left["Cv_t"]), (right_cvt, right["Cv_t"])):
                target[batch] = values
            denominator = np.maximum.reduce([abs(left["Cv"]), abs(right["Cv"]),
                                             abs(left["Cv_t"]), abs(right["Cv_t"])])
            differences = left["Cv_t"] - right["Cv_t"]
            if not np.all(np.isfinite(denominator)) or not np.all(np.isfinite(differences)):
                raise ValueError(f"Nonfinite {direction} face input")
            jump[batch] = np.divide(differences, denominator, out=np.zeros(len(batch)), where=denominator != 0)
        old_scale = np.maximum.reduce([np.ones(len(positions)), abs(left_cv[positions]), abs(right_cv[positions]),
                                       abs(left_cvt[positions]), abs(right_cvt[positions])])
        reproduced = (left_cvt[positions] - right_cvt[positions]) / old_scale
        old_error = np.max(abs(reproduced - old[f"{direction}_Cv_t_jump"][positions]))
        if old_error > 2e-12:
            raise ValueError(f"060 {direction} baseline reproduction failed: {old_error}")
        output.update({f"{direction}_{name}": values for name, values in
                       (("Cv_left", left_cv), ("Cv_right", right_cv), ("Cv_t_left", left_cvt),
                        ("Cv_t_right", right_cvt), ("jump", jump))})
        checks[direction] = {"selected_count": len(positions), "old_jump_max_abs_error": float(old_error),
                             "new_jump_abs_max": float(np.max(abs(jump[positions]))),
                             "new_jump_abs_p99": float(np.quantile(abs(jump[positions]), .99))}
    np.savez_compressed(NEW / "reranked_faces.npz", **output)
    selected = select_cases(manifest, old, output)
    original = selected[0]
    original["logQ"] = float(logq[np.argmin(abs(logq - original["logQ"]))])
    original["logRho"] = original["logQ"] + 2 * original["logT"] - 12
    q_face = np.array([original["logQ"]])
    t_face = np.array([original["logT"]])
    q_index = int(np.argmin(abs(logq - original["logQ"])))
    t_index = int(np.searchsorted(logt, original["logT"]) - 1)
    first = thermodynamic_partials(evaluate(nodes, complete, logq, logt, q_face, t_face,
                                            cell_indices=(np.array([q_index - 1]), np.array([t_index]))), q_face, t_face)
    second = thermodynamic_partials(evaluate(nodes, complete, logq, logt, q_face, t_face,
                                             cell_indices=(np.array([q_index]), np.array([t_index]))), q_face, t_face)
    original["Cv_left"], original["Cv_right"] = float(first["Cv"][0]), float(second["Cv"][0])
    original["Cv_t_left"], original["Cv_t_right"] = float(first["Cv_t"][0]), float(second["Cv_t"][0])
    denominator = max(*(abs(original[key]) for key in ("Cv_left", "Cv_right", "Cv_t_left", "Cv_t_right")))
    original["jump_signed"] = (original["Cv_t_left"] - original["Cv_t_right"]) / denominator if denominator else 0
    original["jump_abs"] = abs(original["jump_signed"])
    original["rho_cgs"] = 10**original["logRho"]
    original["T_kelvin"] = 10**original["logT"]
    original["density_index"] = q_index - 1
    original["temperature_index"] = t_index
    _, _, _, source_rows = read_table(SOURCE)
    face_log_pgas = (source_rows[q_index, t_index, 1] + source_rows[q_index, t_index + 1, 1]) / 2
    original["source_logPrad_over_Pgas_proxy"] = float(np.log10(CRAD / 3) + 4 * original["logT"] - face_log_pgas)
    report = {"construction": "original 060 total potential; floor logT=3.2; radius-two mean twists; no analytic radiation split",
              "normalization": "abs(Cv_t_left-Cv_t_right)/max(abs(Cv_left),abs(Cv_right),abs(Cv_t_left),abs(Cv_t_right)); zero for all-zero finite inputs",
              "selected_population": "060 fully FreeEOS probe routing and supported source faces at tangential midpoints",
              "source_sha256": checksum(SOURCE), "cache_sha256": checksum(CACHE),
              "manifest_sha256": checksum(OLD / "face_manifest.npz"),
              "old_face_results_sha256": checksum(OLD / "face_results.npz"),
              "new_face_results_sha256": checksum(NEW / "reranked_faces.npz"),
              "checks": checks, "cases": selected}
    (NEW / "selection.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"checks": checks, "cases": [(case["id"], case["orientation"], case["logQ"],
                                                    case["logT"], case["jump_abs"]) for case in selected]}, indent=2))


def select_cases(manifest, old, results):
    cases = [{"id": "original", "orientation": "u", "source": "066 saved exact slice",
              "logQ": .77, "logT": 6.756476983450105, "logRho": .77 + 2 * 6.756476983450105 - 12,
              "anchor_is_060_midpoint": False}]
    conditions = (
        ("u_peak", "u", lambda q, t, r: np.ones(len(q), dtype=bool), "largest supported pure u jump"),
        ("u_cool", "u", lambda q, t, r: (t >= 4.5) & (t < 5.5) & (q > -2), "largest cool material u jump in stated window"),
        ("tau_hot", "tau", lambda q, t, r: t >= 6.5, "largest hot tau jump"),
        ("tau_warm", "tau", lambda q, t, r: (t >= 5.5) & (t < 6.5) & (q > -2), "largest warm material tau jump"),
        ("tau_cool", "tau", lambda q, t, r: (t >= 3.5) & (t < 4.5) & (q > -2), "largest cool material tau jump"),
        ("tau_radiation", "tau", lambda q, t, r: (t >= 3.5) & (t < 4.5) & (r > 1), "largest tau jump where source Prad/Pgas proxy exceeds ten"),
    )
    for name, direction, predicate, rule in conditions:
        q, t = manifest[f"{direction}_logQ"], manifest[f"{direction}_logT"]
        radiation = old[f"{direction}_source_logPrad_over_Pgas_proxy"]
        jump = results[f"{direction}_jump"]
        eligible = old[f"{direction}_pure"] & (old[f"{direction}_support_code"] == 0)
        eligible &= predicate(q, t, radiation) & np.isfinite(jump)
        positions = np.flatnonzero(eligible)
        if len(positions) == 0:
            raise ValueError(f"No eligible {name} case")
        index = int(positions[np.argmax(abs(jump[positions]))])
        cases.append(case_record(name, direction, index, rule, manifest, old, results, len(positions)))
    direction = "u"
    q, t = manifest[f"{direction}_logQ"], manifest[f"{direction}_logT"]
    jump = results[f"{direction}_jump"]
    eligible = old[f"{direction}_pure"] & (old[f"{direction}_support_code"] == 0)
    eligible &= (t >= 4.5) & (t < 5.5) & (q > -2) & (q < 0) & np.isfinite(jump)
    positions = np.flatnonzero(eligible)
    target = np.quantile(abs(jump[positions]), .25)
    distance = abs(q[positions] + 1) / .5 + abs(t[positions] - 5) / .5
    score = abs(np.log10(abs(jump[positions]) / target)) + .1 * distance
    index = int(positions[np.argmin(score)])
    cases.append(case_record("quiet", direction, index,
                             "closest to lower-quartile absolute jump, mild distance penalty from (logQ,logT)=(-1,5)",
                             manifest, old, results, len(positions)))
    return cases


def case_record(name, direction, index, rule, manifest, old, results, eligible_count):
    sample_q = float(manifest[f"{direction}_logQ"][index])
    sample_t = float(manifest[f"{direction}_logT"][index])
    jump = float(results[f"{direction}_jump"][index])
    return {"id": name, "orientation": direction, "source": "060 face manifest and independently reranked original construction",
            "selection_rule": rule, "eligible_count": eligible_count, "face_index": index,
            "density_index": int(manifest[f"{direction}_density_index"][index]),
            "temperature_index": int(manifest[f"{direction}_temperature_index"][index]),
            "logQ": sample_q, "logT": sample_t, "logRho": sample_q + 2 * sample_t - 12,
            "rho_cgs": 10**(sample_q + 2 * sample_t - 12), "T_kelvin": 10**sample_t,
            "jump_signed": jump, "jump_abs": abs(jump),
            "Cv_left": float(results[f"{direction}_Cv_left"][index]),
            "Cv_right": float(results[f"{direction}_Cv_right"][index]),
            "Cv_t_left": float(results[f"{direction}_Cv_t_left"][index]),
            "Cv_t_right": float(results[f"{direction}_Cv_t_right"][index]),
            "source_logPrad_over_Pgas_proxy": float(old[f"{direction}_source_logPrad_over_Pgas_proxy"][index]),
            "anchor_pure_probe": bool(old[f"{direction}_pure"][index]),
            "anchor_support_code": int(old[f"{direction}_support_code"][index])}


def make_curves():
    selection = json.loads((NEW / "selection.json").read_text())
    logq, logt, nodes, complete = load_construction()
    _, cache = baseline.read_cache(CACHE, len(logq), len(logt))
    supported_cells = valid_cells(complete)
    probe_lines = []
    metadata = {}
    for case in selection["cases"]:
        name = case["id"]
        if name == "original":
            continue
        anchor_q = case["logQ"]
        anchor_t = case["logT"]
        q_min = max(float(logq[0]), anchor_q - .12)
        q_max = min(float(logq[-1]), anchor_q + .12)
        query_q = np.linspace(q_min, q_max, 401)
        query_t = np.full(len(query_q), anchor_t)
        probe_lines.extend(f"{density:.17g} {anchor_t:.17g}\n" for density in query_q + 2 * anchor_t - 12)
        legacy = legacy_state(cache, logq, logt, query_q, query_t)
        cell_q, cell_t, inside = query_cells(logq, logt, query_q, query_t)
        supported = inside & supported_cells[cell_q, cell_t]
        curves = {"logQ": query_q, "logT": query_t, "logRho": query_q + 2 * anchor_t - 12,
                  "supported": supported}
        for field in ("Cv", "Cv_x", "Cv_t", "dse_Cv"):
            curves[f"legacy_{field}"] = legacy[field]
        for label, forced_t in (("lower", np.searchsorted(logt, anchor_t, side="left") - 1),
                                ("upper", np.searchsorted(logt, anchor_t, side="right") - 1)) if case["orientation"] == "tau" else (("single", None),):
            if forced_t is None:
                chosen = supported
                indices = (cell_q[chosen], cell_t[chosen])
            else:
                chosen = inside & supported_cells[cell_q, forced_t]
                indices = (cell_q[chosen], np.full(int(chosen.sum()), forced_t))
            state = potential_state(nodes, complete, logq, logt, query_q[chosen], query_t[chosen], indices)
            curves[f"{label}_supported"] = chosen
            for field in ("Cv", "Cv_x", "Cv_t", "dse_Cv"):
                values = np.full(len(query_q), np.nan)
                values[chosen] = state[field]
                curves[f"{label}_{field}"] = values
        trace_q_parts = []
        trace_segment_parts = []
        trace_field_parts = {f"legacy_{field}": [] for field in ("Cv", "Cv_x", "Cv_t", "dse_Cv")}
        trace_field_parts.update({f"{label}_{field}": [] for label in
                                  (("lower", "upper") if case["orientation"] == "tau" else ("single",))
                                  for field in ("Cv", "Cv_x", "Cv_t", "dse_Cv")})
        for index in range(len(logq) - 1):
            lower = max(q_min, float(logq[index]))
            upper = min(q_max, float(logq[index + 1]))
            if upper - lower < 1e-12:
                continue
            mid_t = int(np.searchsorted(logt, anchor_t, side="right") - 1)
            cell_t_options = ((mid_t - 1, mid_t) if case["orientation"] == "tau" else (mid_t,))
            if not all(supported_cells[index, local_t] for local_t in cell_t_options):
                continue
            internal = query_q[(query_q > lower + 1e-12) & (query_q < upper - 1e-12)]
            trace_q = np.r_[lower, internal, upper]
            trace_t = np.full(len(trace_q), anchor_t)
            trace_q_parts.append(trace_q)
            trace_segment_parts.append(np.full(len(trace_q), len(trace_q_parts) - 1))
            legacy_trace = legacy_state(cache, logq, logt, trace_q, trace_t)
            for field in ("Cv", "Cv_x", "Cv_t", "dse_Cv"):
                trace_field_parts[f"legacy_{field}"].append(legacy_trace[field])
            for label, local_t in zip((("lower", "upper") if case["orientation"] == "tau" else ("single",)), cell_t_options):
                forced = (np.full(len(trace_q), index), np.full(len(trace_q), local_t))
                state = potential_state(nodes, complete, logq, logt, trace_q, trace_t, forced)
                for field in ("Cv", "Cv_x", "Cv_t", "dse_Cv"):
                    trace_field_parts[f"{label}_{field}"].append(state[field])
        curves["trace_logQ"] = np.concatenate(trace_q_parts)
        curves["trace_segment"] = np.concatenate(trace_segment_parts)
        curves.update({f"trace_{key}": np.concatenate(parts) for key, parts in trace_field_parts.items()})
        np.savez_compressed(NEW / f"{name}_curves.npz", **curves)
        metadata[name] = {"count": len(query_q), "logQ_range": [q_min, q_max],
                          "source_support_count": int(supported.sum()), "source_gap_count": int((~supported).sum()),
                          "exact_cell_trace_segments": len(trace_q_parts),
                          "curves_sha256": checksum(NEW / f"{name}_curves.npz")}
    (NEW / "routing_queries.txt").write_text("".join(probe_lines))
    (NEW / "curve_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("rank", "curves"))
    arguments = parser.parse_args()
    {"rank": rerank, "curves": make_curves}[arguments.stage]()


if __name__ == "__main__":
    main()
