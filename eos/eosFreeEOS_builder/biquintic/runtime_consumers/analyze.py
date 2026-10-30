"""Audit frozen forcing and summarize completed actual-API consumer repetitions."""

from collections import Counter, defaultdict
from itertools import zip_longest
import json
from pathlib import Path

import numpy as np

from run_consumers import EXPERIMENTS, OUTPUT, population, write_json
from runtime_environment import sha256


def load_records(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def reference_root(record):
    return record.get("reference_root") if record["consumer"] == "scalar" else record["reference"].get("root")


def metric(values):
    values = np.asarray([value for value in values if value is not None], dtype=float)
    if not len(values):
        return {"count": 0}
    return dict(count=len(values), median=float(np.median(values)), maximum=float(max(values)), minimum=float(min(values)))


def main():
    output = OUTPUT / "full"
    records = load_records(output / "trials.jsonl")
    original = {"scalar": {record["trial_id"]: record for record in load_records(EXPERIMENTS / "039_thermal_consumer/trials.jsonl")},
                "coupled": {record["trial_id"]: record for record in load_records(EXPERIMENTS / "042_coupled_consumer/trials.jsonl")}}
    groups = defaultdict(list)
    mismatches = []
    sample_count_changes = []
    for record in records:
        prior = original[record["consumer"]][record["trial_id"]]
        keys = (("initial_y", "solver_start_y", "H_erg_per_g_s", "dt_seconds", "bracket", "tolerance_absolute_R", "target_y")
                if record["consumer"] == "scalar" else
                ("initial", "nominal_target", "pressure_target", "heating", "energy_scale", "box"))
        for key in keys:
            if record[key] != prior[key]:
                mismatches.append(dict(consumer=record["consumer"], route=record["route"], trial_id=record["trial_id"], field=key))
        if record["consumer"] == "coupled" and record["solver_start"] != [prior["newton"]["history"][0][key] for key in ("x", "t")]:
            mismatches.append(dict(consumer="coupled", trial_id=record["trial_id"], field="solver_start"))
        if record["consumer"] == "scalar":
            before, after = prior["reference_diagnostic"]["sample_count"], record["reference_diagnostic"].get("sample_count")
            if before != after:
                sample_count_changes.append(dict(trial_id=record["trial_id"], route=record["route"], before=before, after=after))
        groups[record["consumer"], record["route"], record["model"]].append(record)
    grouped = []
    for (consumer, route, model), members in groups.items():
        root_errors = [abs(record["root_error"]) if consumer == "scalar" and record.get("root_error") is not None
                       else record.get("root_error_inf") for record in members]
        shifts = [abs(record["reference_shift"]) if consumer == "scalar" and record.get("reference_shift") is not None
                  else record.get("reference_shift_inf") for record in members]
        conditions = [item["condition"] for record in members for item in record["newton"]["history"] if "condition" in item]
        grouped.append(dict(consumer=consumer, route=route, model=model, count=len(members),
                            statuses=dict(Counter(record["status"] for record in members)),
                            iterations=dict(Counter(record["newton"].get("iterations", 0) for record in members)),
                            iteration_difference_from_original=dict(Counter(record["iteration_difference"] for record in members)),
                            rejected=sum(record["newton"].get("rejected", 0) for record in members),
                            bisected=sum(record["newton"].get("bisected", 0) for record in members),
                            damped_steps=sum(record["newton"].get("damped_steps", 0) for record in members),
                            own_reference_error=metric(root_errors), reference_shift_from_original=metric(shifts),
                            maximum_condition=max(conditions) if conditions else None,
                            selected_face_crossed=sum(record["crossings"]["selected_face_crossed"] for record in members)))
    tolerance_groups = defaultdict(list)
    for record in records:
        if record["consumer"] == "scalar":
            tolerance_groups[record["route"], record["model"], record["tolerance_absolute_R"]].append(record)
    tolerance_summary = [dict(route=route, model=model, tolerance=tolerance, count=len(members),
                              statuses=dict(Counter(record["status"] for record in members)),
                              own_reference_error=metric([abs(record["root_error"]) for record in members if record.get("root_error") is not None]),
                              changed_iterations=[dict(trial_id=record["trial_id"], case_id=record["case_id"],
                                                       difference=record["iteration_difference"], iterations=record["newton"]["iterations"])
                                                  for record in members if record["iteration_difference"] != 0])
                         for (route, model, tolerance), members in tolerance_groups.items()]

    paired = defaultdict(dict)
    for record in records:
        if record["family"] == "common":
            key = (record["consumer"], record["route"], record["case_id"], record["sign"], record["scale"],
                   record.get("density_over_temperature_increment"))
            paired[key][record["model"]] = record
    pair_groups = defaultdict(list)
    for key, pair in paired.items():
        candidate, legacy = pair["candidate"], pair["legacy"]
        pair_groups[key[0], key[1], candidate["category"]].append(dict(
            case_id=candidate["case_id"], iteration_difference=candidate["newton"]["iterations"] - legacy["newton"]["iterations"],
            reference_root_separation=float(np.max(np.abs(np.asarray(reference_root(candidate)) - reference_root(legacy))))))
    pairing = [dict(consumer=consumer, route=route, category=category, pairs=len(members),
                    iteration_differences=dict(Counter(record["iteration_difference"] for record in members)),
                    own_model_reference_separation=metric([record["reference_root_separation"] for record in members]))
               for (consumer, route, category), members in pair_groups.items()]
    route_pairs = defaultdict(dict)
    for record in records:
        route_pairs[record["consumer"], record["model"], record["trial_id"]][record["route"]] = record
    route_differences = []
    for key, pair in route_pairs.items():
        component, final = pair["component"], pair["final"]
        route_differences.append(dict(consumer=key[0], model=key[1], trial_id=key[2],
                                     iteration_difference=final["newton"]["iterations"] - component["newton"]["iterations"],
                                     reference_separation=float(np.max(np.abs(np.asarray(reference_root(final)) - reference_root(component)))),
                                     status_equal=final["status"] == component["status"]))

    actual_populations = Counter()
    errors, invalid_fractions, callback_violations, retention_violations, nonfinite_success_rows = [], [], [], [], []
    state_lookup = {}
    with (output / "api/queries.txt").open() as queries, (output / "api/results.txt").open() as results:
        for index, (query, result) in enumerate(zip_longest(queries, results), 1):
            if query is None or result is None:
                raise AssertionError("Saved query/result stream lengths differ")
            scenario, route, logdensity, logtemperature = query.split()
            raw = np.fromstring(result, sep=" ")
            if len(raw) != 125 or raw[0] != index:
                raise ValueError("Saved API protocol mismatch")
            if raw[1] != 0:
                errors.append(dict(row=index, ierr=int(raw[1]), scenario=scenario, route=route))
                label = "error"
            else:
                if not np.all(np.isfinite(raw[4:])):
                    nonfinite_success_rows.append(index)
                if raw[3] != 1:
                    retention_violations.append(index)
                fractions = raw[23:31]
                label = population(fractions)
                if abs(sum(fractions) - 1) > 1e-12 or min(fractions) < -1e-12 or max(fractions) > 1 + 1e-12:
                    invalid_fractions.append(index)
                expected_callbacks = 0 if scenario == "off" else int(route == "component" or fractions[2] > 0)
                if raw[2] != expected_callbacks:
                    callback_violations.append(index)
            actual_populations[f"{route}/{scenario}/{label}"] += 1
            state_lookup[scenario, route, float(logdensity), float(logtemperature)] = label
    path_populations = Counter()
    blend_states = []
    paths = load_records(output / "states.jsonl")
    for state in paths:
        key = ("candidate" if state["model"] == "candidate" else "off", state["route"], state["logRho"], state["logT"])
        label = state_lookup.get(key, "not_queried_at_exact_serialized_coordinate")
        path_populations[f"{state['consumer']}/{state['route']}/{state['model']}/{state['stage']}/{label}"] += 1
        if label != "pure_FreeEOS":
            blend_states.append(dict(state, population=label))
    write_json(output / "nonpure_state_records.json", blend_states)
    write_json(output / "component_final_differences.json", route_differences)
    derivative_groups = defaultdict(list)
    for record in json.loads((output / "jacobian_checks.json").read_text()):
        derivative_groups[record["consumer"], record["route"], record["kind"]].append(record)
    fd_summary = [dict(consumer=consumer, route=route, kind=kind,
                       statuses=dict(Counter(record["status"] for record in members)),
                       max_scaled_error=max(record["max_scaled_error"] for record in members if record["status"] == "checked"))
                  for (consumer, route, kind), members in derivative_groups.items()]
    failures = [record for record in records if record["status"] != "converged"]
    expected_queries = json.loads((output / "api/client.json").read_text())["actual_queries"]
    if sum(actual_populations.values()) != expected_queries:
        raise AssertionError("Saved protocol count differs from client completion record")
    result = dict(trials=len(records), forcing_mismatch_count=len(mismatches), forcing_mismatches=mismatches,
                  scalar_reference_sample_count_changes=sample_count_changes, groups=grouped, common_pairs=pairing,
                  scalar_by_tolerance=tolerance_summary,
                  actual_query_populations=dict(actual_populations), actual_query_count=sum(actual_populations.values()),
                  actual_api_errors=errors, invalid_fraction_rows=invalid_fractions, callback_count_violations=callback_violations,
                  retention_flag_violations=retention_violations, nonfinite_success_rows=nonfinite_success_rows,
                  callback_audit="Every successful candidate component=1; final positive FreeEOS fraction=1, exactly zero fraction=0; off=0",
                  retention_audit="Require driver full27-by20 callback retention/finite flag=1 and every serialized value/partial finite; unexported internal rows are checked by frozen driver",
                  saved_state_populations=dict(path_populations), fd_summary=fd_summary,
                  failed_trials=[{key: record[key] for key in ("consumer", "route", "trial_id", "case_id", "model", "family", "status")}
                                 for record in failures],
                  trial_file_sha256=sha256(output / "trials.jsonl"),
                  source_or_stellar_accuracy_claim=False, fallback_count=0,
                  note="Reference/model/root differences are separate; finite route samples do not establish continuous-path routing.")
    write_json(output / "analysis.json", result)
    print(json.dumps({key: result[key] for key in ("trials", "forcing_mismatch_count", "actual_query_count", "actual_api_errors", "failed_trials")}, indent=2))
    if errors or invalid_fractions or callback_violations or retention_violations or nonfinite_success_rows or mismatches:
        raise AssertionError("Independent stored-protocol/forcing audit failed; inspect analysis.json")


if __name__ == "__main__":
    main()
