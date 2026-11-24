"""Native switch/API checks; this is not a stellar-evolution test."""

import argparse
import io
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

SOURCE = Path(__file__).resolve().parent
sys.path.insert(0, str(SOURCE.parent / "runtime_bridge"))
from runtime_environment import (MESA, MASSFILE, WORKSPACE, installed_inventory,
                                 private_environment, sha256)

OUTPUT = WORKSPACE / "freeeos_biquintic_progress/experiments/065_native_integration/api"
BUILD = OUTPUT.parent / "build"
FAMILY = WORKSPACE / "freeeos_biquintic_progress/experiments/062_family_cache/family"
DRIVER = BUILD / "driver/native_driver"
OLD_DRIVER = WORKSPACE / "freeeos_biquintic_progress/experiments/050_runtime_api/build/driver/baseline_driver"
CRAD = 7.565733250280004e-15


def run_case(name, queries, enabled=True, missing=False, old=False, conflict=False):
    directory = OUTPUT / name
    directory.mkdir(parents=True, exist_ok=False)
    environment = private_environment(directory)
    cache_directory = directory / "runtime_cache/eosDT_cache"
    expected_links = {path.name: path.resolve() for path in cache_directory.glob("*.bin")}
    inlist = directory / "inlist_eos"
    controls = "&eos\n use_WDM_AV=.false.\n"
    if not old:
        controls += f" FreeEOS_use_HELM_interpolation={'.true.' if enabled else '.false.'}\n"
        controls += f" FreeEOS_biquintic_table_dir='{FAMILY if not missing else directory / 'missing'}'\n"
        if conflict:
            controls += " use_other_FreeEOS_results=.true.\n"
    controls += "/\n"
    inlist.write_text(controls)
    records = []
    for route, density, temperature, hydrogen, metals in queries:
        values = (density, temperature) if old else (density, temperature, hydrogen, metals)
        records.append(route + " " + " ".join(f"{value:.17g}" for value in values))
    inputs = "\n".join(records) + "\n"
    (directory / "queries.txt").write_text(inputs)
    command = [str(OLD_DRIVER if old else DRIVER), str(MESA), str(inlist), str(MASSFILE)]
    (directory / "command.json").write_text(json.dumps(command, indent=2) + "\n")
    result = subprocess.run(command, input=inputs, text=True, capture_output=True,
                            env=environment, timeout=180)
    (directory / "stdout.txt").write_text(result.stdout)
    (directory / "stderr.txt").write_text(result.stderr)
    (directory / "returncode.txt").write_text(str(result.returncode) + "\n")
    if result.returncode:
        raise RuntimeError(f"{name} process exit {result.returncode}; preserved logs")
    rows = np.loadtxt(io.StringIO(result.stdout), ndmin=2)
    if len(rows) != len(queries):
        raise AssertionError(f"{name}: incomplete output")
    for name, target in expected_links.items():
        link = cache_directory / name
        if not link.is_symlink() or link.resolve() != target:
            raise RuntimeError(f"Existing private cache link replaced: {link}")
    created = sorted(path.name for path in cache_directory.glob("*.bin") if path.name not in expected_links)
    (directory / "new_private_legacy_caches.json").write_text(json.dumps(created, indent=2) + "\n")
    return rows


def main():
    global OUTPUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    OUTPUT = parser.parse_args().output.resolve()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    before = installed_inventory()
    metrics = {"status": "running", "driver_sha256": sha256(DRIVER),
               "manifest_sha256": sha256(FAMILY / "family.manifest")}
    metrics_path = OUTPUT / "results.json"
    metrics_path.write_text(json.dumps(metrics, indent=2) + "\n")
    solar = [(route, density, temperature, .7, .02)
             for density, temperature in [(-3.017, 5.013), (1.37, 6.503), (4., 6.)]
             for route in ("component", "final")]
    off = run_case("switch_off", solar, enabled=False)
    old = run_case("installed_legacy", solar, old=True)
    np.testing.assert_array_equal(off[:, 1], old[:, 1])
    successful = off[:, 1] == 0
    np.testing.assert_allclose(off[successful, 5:], old[successful, 4:], rtol=5e-14, atol=1e-14)
    assert np.all(off[:, 2:5] == 0)
    metrics["switch_off_legacy_max_scaled_difference"] = float(np.max(
        abs(off[successful, 5:] - old[successful, 4:]) / np.maximum(1, abs(old[successful, 4:]))))
    metrics["switch_off_statuses"] = off[:, 1].tolist()
    points = [(-3.017, 5.013), (-1.017, 5.513), (1.37, 6.503)]
    mixture_queries = [("component", density, temperature, hydrogen, metals)
                       for density, temperature in points
                       for hydrogen, metals in [(.7, .02), (.683, .021), (.52, .035)]]
    rows = run_case("composition", mixture_queries)
    assert np.all(rows[:, 1] == 0), rows[:, :5]
    assert np.all(np.isfinite(rows))
    residuals = []
    for query, row in zip(mixture_queries, rows):
        _, density, temperature, _, _ = query
        rho, temp = 10**density, 10**temperature
        values, density_partials, thermal_partials = row[5:32], row[32:59], row[59:86]
        gas, energy, entropy = np.exp(values[:3])
        radiation = CRAD * temp**4 / 3
        pressure = gas + radiation
        pressure_t = gas * thermal_partials[0] + 4 * radiation
        residuals.append([(rho * energy * density_partials[1] + pressure_t - pressure) / pressure,
                          temp * entropy * thermal_partials[2] / (energy * thermal_partials[1]) - 1,
                          (rho * temp * entropy * density_partials[2] + pressure_t) / pressure_t,
                          values[10] * temp / (energy * thermal_partials[1]) - 1])
    metrics["maximum_identity_residual"] = float(np.max(abs(np.asarray(residuals))))
    assert metrics["maximum_identity_residual"] < 2e-12, residuals
    metrics["composition_success_count"] = int(rows[-1, 2])
    metrics["composition_table_loads"] = int(rows[-1, 4])
    step = 1e-5
    derivative_queries = []
    for route in ("component", "final"):
        for hydrogen, metals in [(.683, .021), (.52, .035)]:
            for delta_x, delta_z in [(0, 0), (step, 0), (-step, 0), (0, step), (0, -step)]:
                derivative_queries.append((route, *points[0], hydrogen + delta_x, metals + delta_z))
    derivative_rows = run_case("abundance_partials", derivative_queries)
    assert np.all(derivative_rows[:, 1] == 0)
    differences = []
    for offset in range(0, len(derivative_rows), 5):
        center, xplus, xminus, zplus, zminus = derivative_rows[offset:offset + 5]
        expected_x = (xplus[5:7] - xminus[5:7]) / (2 * step)
        expected_z = (zplus[5:7] - zminus[5:7]) / (2 * step)
        analytic_x = center[86:88] - center[88:90]
        analytic_z = center[90:92] - center[88:90]
        differences.extend(abs(analytic_x - expected_x) / np.maximum(1, abs(expected_x)))
        differences.extend(abs(analytic_z - expected_z) / np.maximum(1, abs(expected_z)))
        np.testing.assert_allclose(analytic_x, expected_x, rtol=2e-6, atol=2e-8)
        np.testing.assert_allclose(analytic_z, expected_z, rtol=2e-6, atol=2e-8)
    metrics["composition_partial_max_scaled_difference"] = float(max(differences))
    pure = [(route, *points[0], .683, .021) for route in ("component", "final")]
    missing = run_case("missing_family", pure, missing=True)
    assert np.all(missing[:, 1] == -7950), missing[:, :5]
    cold = run_case("unsupported_cells", [(route, -10., 3.19, .7, .02)
                                         for route in ("component", "final")])
    assert np.all(cold[:, 1] == -7950), cold[:, :5]
    conflict = run_case("conflicting_hooks", pure, conflict=True)
    assert np.all(conflict[:, 1] == -7903), conflict[:, :5]
    assert installed_inventory() == before
    metrics.update(status="passed", missing_status=missing[:, 1].tolist(),
                   cold_status=cold[:, 1].tolist(), conflict_status=conflict[:, 1].tolist(),
                   installed_inventory_unchanged=True)
    metrics_path.write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
