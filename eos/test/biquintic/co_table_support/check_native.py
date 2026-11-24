import csv
import json
import math
from pathlib import Path
import subprocess
import sys


SOURCE = Path(__file__).resolve().parent
WORKSPACE = SOURCE.parents[4]
sys.dont_write_bytecode = True
sys.path.insert(0, str(SOURCE.parent / "runtime_bridge"))
from runtime_environment import MESA, MASSFILE, installed_inventory, private_environment, sha256


EXPERIMENT = WORKSPACE / "freeeos_biquintic_progress/experiments/075_co_table_support"
ARCHIVE = WORKSPACE / "freeeos_biquintic_progress/experiments/073_mesa_test_suite"
FAMILY = WORKSPACE / "freeeos_biquintic_progress/experiments/062_family_cache/family"
DRIVER = EXPERIMENT / "native_build/driver/native_driver"
DUAL_DRIVER = EXPERIMENT / "native_build/driver/dual_handle_driver"
OLD_DRIVER = WORKSPACE / "freeeos_biquintic_progress/experiments/065_native_integration/build/driver/native_driver"
CARBON = ARCHIVE / "z1_default_routing/mass_frac_carbon.txt"


def run_case(name, queries, enabled=True, ordinary=False, unknown=False, table_directory=FAMILY,
             massfile=MASSFILE, old=False, expect_process_failure=False):
    directory = EXPERIMENT / "native_api" / name
    directory.mkdir(parents=True, exist_ok=True)
    environment = private_environment(directory)
    controls = ["&eos", " use_WDM_AV=.false.",
                f" FreeEOS_use_HELM_interpolation={'.true.' if enabled else '.false.'}",
                f" FreeEOS_biquintic_table_dir='{table_directory}'"]
    if ordinary:
        controls.extend(f" suffix_for_FreeEOS_Z({layer})=''" for layer in (13, 14, 15))
    if unknown:
        controls.append(" suffix_for_FreeEOS_Z(13)='_XX'")
    inlist = directory / "inlist_eos"
    inlist.write_text("\n".join(controls + ["/"]) + "\n")
    inputs = "\n".join(" ".join([route] + [f"{number:.17g}" for number in state])
                       for route, *state in queries) + "\n"
    (directory / "queries.txt").write_text(inputs)
    command = [str(OLD_DRIVER if old else DRIVER), str(MESA), str(inlist), str(massfile)]
    (directory / "command.json").write_text(json.dumps(command, indent=2) + "\n")
    result = subprocess.run(command, input=inputs, text=True, capture_output=True,
                            env=environment, timeout=180)
    (directory / "stdout.txt").write_text(result.stdout)
    (directory / "stderr.txt").write_text(result.stderr)
    (directory / "returncode.txt").write_text(str(result.returncode) + "\n")
    if expect_process_failure and result.returncode and "missing or bad eos data" in result.stderr:
        return []
    if result.returncode:
        raise RuntimeError(f"{name}: driver exit {result.returncode}; see logs")
    rows = [[float(value) for value in line.split()] for line in result.stdout.splitlines()
            if line[:1].isdigit()]
    if len(rows) != len(queries):
        raise RuntimeError(f"{name}: incomplete output")
    return rows


def statuses(rows):
    return [int(row[1]) for row in rows]


def check_dual_handles(configured_pressure, ordinary_pressure):
    directory = EXPERIMENT / "native_api/dual_handles"
    directory.mkdir(parents=True, exist_ok=True)
    environment = private_environment(directory)
    command = [str(DUAL_DRIVER), str(MESA),
               str(EXPERIMENT / "native_api/configured/inlist_eos"),
               str(EXPERIMENT / "native_api/ordinary_override/inlist_eos"), str(MASSFILE)]
    (directory / "command.json").write_text(json.dumps(command, indent=2) + "\n")
    result = subprocess.run(command, text=True, capture_output=True, env=environment, timeout=180)
    (directory / "stdout.txt").write_text(result.stdout)
    (directory / "stderr.txt").write_text(result.stderr)
    (directory / "returncode.txt").write_text(str(result.returncode) + "\n")
    if result.returncode:
        raise RuntimeError("dual handle driver failed")
    lines = [line for line in result.stdout.splitlines()
             if line.startswith("COUNTS") or line[:1].isdigit()]
    if len(lines) != 7:
        raise AssertionError("dual handle output incomplete")
    for iteration, line in enumerate(lines[:6], 1):
        number, status, pressure = line.split()
        expected = configured_pressure if iteration % 2 else ordinary_pressure
        if int(number) != iteration or int(status) != 0 or abs(float(pressure)-expected) > 1e-12:
            raise AssertionError("dual handle cache contamination")
    if lines[-1] != "COUNTS 6 0 4":
        raise AssertionError("dual handle load accounting")
    return {"statuses": [0] * 6, "successes": 6, "failures": 0, "loads": 4}


def main():
    before = installed_inventory()
    metrics = {"driver_sha256": sha256(DRIVER),
               "family_manifest_sha256": sha256(FAMILY / "family.manifest")}
    queries = [("component", -7.0, 5.0, 0.0, metals)
               for metals in (0.02, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0)]
    selected = run_case("configured", queries)
    ordinary = run_case("ordinary_override", queries, ordinary=True)
    metrics["configured_statuses"] = statuses(selected)
    metrics["ordinary_statuses"] = statuses(ordinary)
    if any(status != 0 for status in metrics["configured_statuses"] + metrics["ordinary_statuses"]):
        raise AssertionError("High-Z native component query failed; inspect preserved logs")
    pressure_differences = [abs(selected[index][5] - ordinary[index][5])
                            for index in range(len(queries))]
    metrics["lnPgas_configured_ordinary_differences"] = pressure_differences
    if pressure_differences[0] != 0 or not all(difference > 1e-8 for difference in pressure_differences[1:]):
        raise AssertionError("Configuration does not change only selected high-Z stencils")
    metrics["dual_handles"] = check_dual_handles(selected[3][5], ordinary[3][5])
    step = 1e-5
    derivative_queries = [("component", -7.0, 5.0, hydrogen, metals)
                          for hydrogen, metals in ((0.05, 0.85), (0.05+step, 0.85),
                                                   (0.05-step, 0.85), (0.05, 0.85+step),
                                                   (0.05, 0.85-step))]
    derivative_rows = run_case("high_z_abundance_partials", derivative_queries)
    if statuses(derivative_rows) != [0] * 5:
        raise AssertionError("High-Z native derivative queries failed")
    center, hydrogen_plus, hydrogen_minus, metals_plus, metals_minus = derivative_rows
    residuals = []
    for field_offset in range(2):
        finite_hydrogen = (hydrogen_plus[5+field_offset] - hydrogen_minus[5+field_offset]) / (2*step)
        finite_metals = (metals_plus[5+field_offset] - metals_minus[5+field_offset]) / (2*step)
        analytic_hydrogen = center[86+field_offset] - center[88+field_offset]
        analytic_metals = center[90+field_offset] - center[88+field_offset]
        for actual, expected in ((analytic_hydrogen, finite_hydrogen),
                                 (analytic_metals, finite_metals)):
            residual = abs(actual-expected) / max(1.0, abs(expected))
            residuals.append(residual)
            if residual > 2e-4:
                raise AssertionError("High-Z native abundance derivative mismatch")
    metrics["high_z_abundance_partial_max_scaled_difference"] = max(residuals)
    off = run_case("switch_off", queries, enabled=False)
    old = run_case("previous_private_switch_off", queries, enabled=False, old=True)
    metrics["switch_off_statuses"] = statuses(off)
    if statuses(off) != statuses(old):
        raise AssertionError("Switch-off status changed")
    for current, previous in zip(off, old):
        for value, expected in zip(current[5:], previous[5:]):
            if math.isnan(value) and math.isnan(expected):
                continue
            if value != expected:
                raise AssertionError("Switch-off output changed")
    metrics["switch_off_exact_match"] = True
    missing = run_case("missing_variant", queries[2:3],
                       table_directory=EXPERIMENT / "fixtures/missing")
    unknown = run_case("unsupported_suffix", queries[2:3], unknown=True,
                       expect_process_failure=True)
    metrics["missing_variant_status"] = statuses(missing)
    metrics["unsupported_suffix_legacy_loader_fatal"] = not unknown
    if statuses(missing) != [-7950] or unknown:
        raise AssertionError("Strict native failure did not propagate")
    original = json.loads((ARCHIVE / "z1_default_routing/summary.json").read_text())
    original_queries = [("component", item["logRho"], item["logT"], 0.0, 1.0)
                        for item in original["ns_c_default_fraction_queries"]]
    with (ARCHIVE / "tightened_analysis/printed_failures.tsv").open() as stream:
        tightened_states = [item for item in csv.DictReader(stream, delimiter="\t")
                            if item["name"] == "ns_c"]
    tightened_queries = [("component", float(item["logRho"]), float(item["logT"]), 0.0, 1.0)
                         for item in tightened_states]
    original_rows = run_case("ns_c_original", original_queries, massfile=CARBON)
    tightened_rows = run_case("ns_c_tightened", tightened_queries, massfile=CARBON)
    original_ordinary = run_case("ns_c_original_ordinary", original_queries,
                                 ordinary=True, massfile=CARBON)
    tightened_ordinary = run_case("ns_c_tightened_ordinary", tightened_queries,
                                  ordinary=True, massfile=CARBON)
    original_final = run_case("ns_c_original_final",
                              [("final", *query[1:]) for query in original_queries], massfile=CARBON)
    tightened_final = run_case("ns_c_tightened_final",
                               [("final", *query[1:]) for query in tightened_queries], massfile=CARBON)
    metrics["ns_c_original_statuses"] = statuses(original_rows)
    metrics["ns_c_tightened_statuses"] = statuses(tightened_rows)
    metrics["ns_c_original_ordinary_statuses"] = statuses(original_ordinary)
    metrics["ns_c_tightened_ordinary_statuses"] = statuses(tightened_ordinary)
    metrics["ns_c_original_final_statuses"] = statuses(original_final)
    metrics["ns_c_tightened_final_statuses"] = statuses(tightened_final)
    metrics["ns_c_original_inputs"] = [query[1:3] for query in original_queries]
    metrics["ns_c_tightened_inputs"] = [query[1:3] for query in tightened_queries]
    if installed_inventory() != before:
        raise AssertionError("Installed EOS library or cache inventory changed")
    metrics["installed_inventory_unchanged"] = True
    (EXPERIMENT / "native_api/results.json").write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
