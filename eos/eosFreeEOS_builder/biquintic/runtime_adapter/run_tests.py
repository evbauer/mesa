"""Bounded standalone adapter validation; writes only the TASK047 experiment."""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import socket
import subprocess
import sys
import time

import numpy as np


WORKSPACE = Path(__file__).resolve().parents[5]
MESA = WORKSPACE / "mesa"
EXPERIMENT = WORKSPACE / "freeeos_biquintic_progress/experiments/047_runtime_adapter"
PATTERN = MESA / "eos/eosFreeEOS_builder/mass_frac_solar.txt"
PATTERN_SHA = "bd62922804e073a2a410480994b8a28e85b0bc53e469729d6ed5136ea9da10a6"
SOURCE_SHA = "e6ca7ce27bee8197bb257a8c4886a06366e6951e58958fbed3735c877ee68502"
PROCESS_SHA = "40b513e48b512e44ea2c020cd29c000b60128661260f38c0c0b7365c2ae10a7d"
CACHE_SHA = "e2c72c6f02cb900649a6a629bfd4fc53c2e8423e5bbe6afd6a189c45cac3a204"
CRAD = 7.565733250280004e-15
LN10 = math.log(10.0)
THERMO = np.array([0, 1, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])
ORDERS = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2), (3, 0), (2, 1), (1, 2), (0, 3)]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def identity_text(cache, expected_digest=None):
    with Path(cache).open() as stream:
        header = [next(stream).rstrip() for _ in range(6)]
    assert header[0] == "FREEEOS_BIQUINTIC_V1"
    return "\n".join([expected_digest or digest(cache), *header[1:], f"3.2 {PATTERN_SHA}"]) + "\n"


def synthetic_cache(path, invalid=False, rough=False):
    lines = ["FREEEOS_BIQUINTIC_V1", "5 4", f"0.7 0.02 {CRAD:.17e} 1", SOURCE_SHA,
             f"{PROCESS_SHA} 2 mean", "-12.0 0.3 11.5 0.06"]
    for temp_index in range(4):
        natural_temp = 11.5 + 0.06 * temp_index
        for density_index in range(5):
            native = -12.0 + 0.3 * density_index
            values = []
            for temp_order in range(3):
                for density_order in range(3):
                    if density_order == 0:
                        factor = native + 0.5 * natural_temp - 12 * LN10 - 10 + 0.5 * temp_order
                    elif density_order == 1:
                        factor = 1.0
                    else:
                        factor = 0.0
                    value = 2e7 * math.exp(natural_temp) * factor
                    if rough and density_order == 2 and temp_order == 0:
                        value += 2e7 * math.exp(natural_temp) * 0.02 * (-1)**density_index
                    values.append(value)
            valid = not (invalid and density_index == 2 and temp_index == 1)
            if not valid:
                values = [0.0] * 9
            lines.append(f"{density_index+1} {temp_index+1} {int(valid)} " +
                         " ".join(f"{value:.17e}" for value in values))
    path.write_text("\n".join(lines) + "\n")


def analytic_partials(logdensity, logtemperature, radiation=CRAD, manufactured=False):
    temperature = math.exp(logtemperature)
    gas = 2e7
    radiation_f = -radiation * math.exp(4 * logtemperature - logdensity) / 3
    potential = []
    for density_order, temp_order in ORDERS:
        if density_order == 0:
            value = gas * temperature * (logdensity - 1.5 * logtemperature - 10 - 1.5 * temp_order)
        elif density_order == 1:
            value = gas * temperature
        else:
            value = 0.0
        value += radiation_f * (-1)**density_order * 4**temp_order
        if manufactured:
            value += 5e15 * math.exp(2 * logdensity / 3) * (2 / 3)**density_order * (temp_order == 0)
        potential.append(value)
    return potential


def analytic_values(logdensity, logtemperature, radiation=CRAD, manufactured=False):
    density, temperature = math.exp(logdensity), math.exp(logtemperature)
    gas = 2e7
    radiation_pressure = radiation * temperature**4 / 3
    degeneracy_energy = 5e15 * density**(2 / 3) if manufactured else 0.0
    gas_pressure = density * gas * temperature + (2 / 3) * density * degeneracy_energy
    pressure = gas_pressure + radiation_pressure
    energy = 1.5 * gas * temperature + 3 * radiation_pressure / density + degeneracy_energy
    entropy = gas * (11.5 + 1.5 * logtemperature - logdensity) + 4 * radiation_pressure / (density * temperature)
    cv = 1.5 * gas + 12 * radiation_pressure / (density * temperature)
    chi_density = (density * gas * temperature + (10 / 9) * density * degeneracy_energy) / pressure
    chi_temp = (density * gas * temperature + 4 * radiation_pressure) / pressure
    gamma_minus1 = pressure * chi_temp / (density * temperature * cv)
    gamma1 = chi_density + chi_temp * gamma_minus1
    cp = cv + pressure * chi_temp**2 / (density * temperature * chi_density)
    values = np.zeros(27)
    values[:3] = [math.log(gas_pressure), math.log(energy), math.log(entropy)]
    values[3:6] = np.arange(4, 7) / 10
    values[6:16] = [gamma_minus1 / gamma1, chi_density, chi_temp, cp, cv,
                   -3 * radiation_pressure / density**2 + (2 / 3) * degeneracy_energy / density,
                   cv / temperature, -gas / density - 4 * radiation_pressure / (density**2 * temperature),
                   gamma1, 1 + gamma_minus1]
    values[21] = 1
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--driver", type=Path, default=EXPERIMENT / "build/runtime_driver")
    parser.add_argument("--output", type=Path, default=EXPERIMENT / "validation")
    arguments = parser.parse_args()
    expected_node = os.environ.get("SLURMD_NODENAME")
    if not expected_node or socket.gethostname().split(".")[0] != expected_node:
        raise RuntimeError("Compute-node placement is required")
    output = arguments.output.resolve()
    if not output.is_relative_to(EXPERIMENT.resolve()):
        raise ValueError("Output must be in the TASK047 experiment")
    output.mkdir(parents=True, exist_ok=True)
    environment = dict(os.environ, OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", MKL_NUM_THREADS="2",
                       NUMEXPR_NUM_THREADS="2", PYTHONDONTWRITEBYTECODE="1")
    checks, commands = [], []
    started = time.time()

    def check(condition, name, **details):
        checks.append(dict(name=name, passed=bool(condition), **details))
        print(f"{'PASS' if condition else 'FAIL'} {name}", flush=True)

    def invoke(mode, stdin="", cache=None, identity=None, pattern=PATTERN, suffix=None):
        command = [str(arguments.driver.resolve()), mode]
        if mode == "hash":
            command.append(str(cache))
        elif mode != "map":
            command += [str(MESA), str(cache), str(identity), str(pattern)]
        completed = subprocess.run(command, input=stdin, text=True, capture_output=True, env=environment, timeout=90)
        commands.append(dict(command=command, stdin=stdin, returncode=completed.returncode))
        label = suffix or f"{len(commands):03d}_{mode}"
        (output / f"{label}.stdout").write_text(completed.stdout)
        (output / f"{label}.stderr").write_text(completed.stderr)
        if completed.returncode:
            raise RuntimeError(f"{label}: {completed.stderr}")
        return completed.stdout.splitlines()

    def query(cache, identity, states, mode="query"):
        lines = invoke(mode, "".join(f"{alteration} {logrho:.17e} {logtemp:.17e}\n"
                                    for alteration, logrho, logtemp in states), cache, identity)
        data, unchanged = [], []
        for line in lines:
            if line.startswith("INIT"):
                assert line == "INIT 0", line
            elif line.startswith("QUERY"):
                unchanged.append(line.split()[-1] == "T")
            else:
                data.append(np.fromstring(line, sep=" "))
        return np.array(data), unchanged

    def map_states(states):
        records = []
        for logdensity, logtemperature, radiation, manufactured in states:
            values = [math.exp(logdensity), math.exp(logtemperature), radiation,
                      *analytic_partials(logdensity, logtemperature, radiation, manufactured)]
            records.append(" ".join(f"{value:.17e}" for value in values))
        return np.array([np.fromstring(line, sep=" ") for line in invoke("map", "\n".join(records) + "\n")])

    for length in [0, 3, 55, 56, 63, 64, 65, 119, 120, 127, 128, 65535, 65536, 65537, 131137]:
        path = output / f"hash_{length}.bin"
        path.write_bytes(bytes((position * 131 + 17) % 256 for position in range(length)))
        result = invoke("hash", cache=path)[0].split()
        check(result == ["0", digest(path)], f"sha256_{length}_bytes")

    analytic_states = [(math.log(1e-7), math.log(1e5), radiation, manufactured)
                       for radiation, manufactured in [(0.0, False), (CRAD, False), (CRAD, True)]]
    analytic_output = map_states(analytic_states)
    hybrid = np.array([[100 * slot + species for slot in range(1, 28)] for species in range(1, 21)]).ravel()
    for state_index, state in enumerate(analytic_states):
        actual = analytic_output[state_index]
        expected = analytic_values(*state)
        check(actual[0] == 0 and np.allclose(actual[1:28], expected, rtol=3e-13, atol=1e-10),
              f"analytic_27_values_{state_index}")
        check(np.array_equal(actual[82:], hybrid), f"analytic_all27_hybrid_rows_{state_index}")
        check(np.array_equal(actual[28:55][3:6], np.arange(4, 7) / 100) and
              np.array_equal(actual[55:82][3:6], -np.arange(4, 7) / 1000), f"legacy_aux_partials_{state_index}")
        check(np.array_equal(actual[28:55][16:], np.zeros(11)) and
              np.array_equal(actual[55:82][16:], np.zeros(11)), f"phase_fraction_partials_{state_index}")
        for axis in range(2):
            step = 2e-5
            lower, upper = list(state), list(state)
            lower[axis] -= step
            upper[axis] += step
            numerical = (analytic_values(*upper) - analytic_values(*lower)) / (2 * step)
            derivative = actual[28 + axis * 27:55 + axis * 27]
            scale = np.maximum.reduce([np.abs(numerical), np.abs(derivative), np.full(27, 1e-5)])
            error = np.max(np.abs(numerical[THERMO] - derivative[THERMO]) / scale[THERMO])
            check(error < 2e-6, f"analytic_all13_state_partials_{state_index}_{axis}", max_scaled_error=float(error))

    valid_record = [1e-7, 1e5, CRAD, *analytic_partials(math.log(1e-7), math.log(1e5))]
    invalid_records = []
    invalid_names = []
    for name, index, replacement, expected_status in [
            ("zero_density", 0, 0, 103), ("zero_temperature", 1, 0, 103),
            ("negative_radiation", 2, -CRAD, 103), ("nonfinite_potential", 3, math.nan, 14),
            ("nonpositive_pressure", 4, -1, 18),
            ("nonpositive_gas_pressure", 4, CRAD * 1e20 / (6e-7), 18),
            ("nonpositive_energy", 3, valid_record[5] - 1e7, 18),
            ("nonpositive_entropy", 5, 1e7, 18),
            ("nonpositive_cv", 8, valid_record[5] + 1e7, 18),
            ("nonpositive_chirho", 6, -valid_record[4], 18),
            ("nonpositive_chit", 7, 0, 18)]:
        altered = valid_record.copy()
        altered[index] = replacement
        invalid_records.append(" ".join(f"{value:.17e}" for value in altered))
        invalid_names.append((name, expected_status))
    invalid_data = [np.fromstring(line, sep=" ") for line in invoke("map", "\n".join(invalid_records) + "\n")]
    seeded = np.concatenate([np.arange(1, 28) / 10, np.arange(1, 28) / 100, -np.arange(1, 28) / 1000, hybrid])
    for (name, expected_status), actual in zip(invalid_names, invalid_data):
        check(actual[0] == expected_status and np.array_equal(actual[1:], seeded), f"thermodynamic_guard_{name}_atomic")

    cache = output / "analytic.cache"
    synthetic_cache(cache)
    identity = output / "analytic.identity"
    identity.write_text(identity_text(cache))
    base_logtemp = 11.59 / LN10
    base_logrho = (-11.45 + 2 * 11.59 - 12 * LN10) / LN10
    states = [(alteration, base_logrho, base_logtemp) for alteration in range(40)]
    data, unchanged = query(cache, identity, states)
    check(np.array_equal(data[:2, 0], [0, 0]), "pinned_and_reordered_composition_accepted")
    check(np.array_equal(data[0], data[1]), "network_order_does_not_change_results")
    expected_codes = {**dict.fromkeys(range(2, 13), 102), **dict.fromkeys(range(13, 17), 103),
                      **dict.fromkeys(range(17, 21), 104), **dict.fromkeys(range(21, 26), 102), 26: 103, 27: 103,
                      **dict.fromkeys(range(28, 37), 102), **dict.fromkeys(range(37, 40), 103)}
    for alteration, code in expected_codes.items():
        check(data[alteration, 0] == code and unchanged[alteration], f"query_guard_{alteration}_atomic_failure")
    check(np.array_equal(data[0, 82:], hybrid), "query_retains_all27_hybrid_rows")
    step = 2e-5
    fd_states = [(0, base_logrho, base_logtemp)]
    for axis in range(2):
        for sign in [-1, 1]:
            fd_states.append((0, base_logrho + (axis == 0) * sign * step / LN10,
                              base_logtemp + (axis == 1) * sign * step / LN10))
    fd_data, _ = query(cache, identity, fd_states)
    for axis in range(2):
        numerical = (fd_data[2 + 2 * axis, 1:28] - fd_data[1 + 2 * axis, 1:28]) / (2 * step)
        derivative = fd_data[0, 28 + axis * 27:55 + axis * 27]
        scale = np.maximum.reduce([np.abs(numerical), np.abs(derivative), np.full(27, 1e-4)])
        error = np.max(np.abs(numerical[THERMO] - derivative[THERMO]) / scale[THERMO])
        check(np.all(fd_data[:, 0] == 0) and error < 2e-4, f"cached_all13_state_partials_{axis}", max_scaled_error=float(error))

    rough_cache = output / "rough.cache"
    synthetic_cache(rough_cache, rough=True)
    rough_identity = output / "rough.identity"
    rough_identity.write_text(identity_text(rough_cache))
    face_rho = (-11.4 + 2 * 11.59 - 12 * LN10) / LN10
    face_states = []
    step = 2e-6
    for side in [-1, 1]:
        for offset in [1e-9, 1e-9 + step, 1e-9 + 2 * step]:
            face_states.append((0, face_rho + side * offset / LN10, base_logtemp))
    face_data, _ = query(rough_cache, rough_identity, face_states)
    for side_index, side in enumerate([-1, 1]):
        base, next_state, last_state = face_data[3 * side_index:3 * side_index + 3]
        numerical = side * (-3 * base[1:28] + 4 * next_state[1:28] - last_state[1:28]) / (2 * step)
        derivative = base[28:55]
        scale = np.maximum.reduce([np.abs(numerical), np.abs(derivative), np.abs(base[1:28]) * 1e-5])
        error = np.max(np.abs(numerical[THERMO] - derivative[THERMO]) / np.maximum(scale[THERMO], 1e-4))
        check(np.all(face_data[:, 0] == 0) and error < 2e-3, f"face_one_sided_all13_partials_{side}",
              max_scaled_error=float(error), inward_offset_lnRho=1e-9, step_lnRho=step)
    jump = abs(face_data[0, 28 + 7] - face_data[3, 28 + 7])
    check(jump > 1e-3, "manufactured_face_has_resolved_third_derivative_jump", chiRho_partial_jump=float(jump))
    thread_lines = invoke("threads", f"{base_logrho:.17e} {base_logtemp:.17e}\n", cache, identity)
    check(thread_lines[-1] == "THREADS 2 1000 0", "immutable_two_thread_1000_queries")
    empty, unchanged = query(cache, identity, [(0, base_logrho, base_logtemp)], mode="uninitialized")
    check(empty[0, 0] == 100 and unchanged[0], "uninitialized_rejected")
    outside, unchanged = query(cache, identity, [(0, base_logrho, 3.19), (0, -100, base_logtemp),
                                                (0, base_logrho, 6.0)])
    check(np.all(outside[:, 0] == 15) and all(unchanged), "floor_and_grid_rejections")
    invalid_cache = output / "invalid.cache"
    synthetic_cache(invalid_cache, invalid=True)
    invalid_identity = output / "invalid.identity"
    invalid_identity.write_text(identity_text(invalid_cache))
    invalid, unchanged = query(invalid_cache, invalid_identity, [(0, base_logrho, base_logtemp)])
    check(invalid[0, 0] == 16 and unchanged[0], "masked_cell_rejected")

    lines = cache.read_text().splitlines()
    mutations = {}
    changed = lines.copy()
    fields = changed[6].split()
    fields[3] = f"{float(fields[3]) * 1.000001:.17e}"
    changed[6] = " ".join(fields)
    mutations["payload_only"] = "\n".join(changed) + "\n"
    changed = lines.copy()
    changed[6] = "1 1 0 " + " ".join(["0"] * 9)
    mutations["mask_only"] = "\n".join(changed) + "\n"
    mutations["truncated"] = "\n".join(lines[:-1]) + "\n"
    mutations["trailing"] = cache.read_text() + "extra\n"
    for name, contents in mutations.items():
        corrupt = output / f"corrupt_{name}.cache"
        corrupt.write_text(contents)
        result = invoke("initialize", cache=corrupt, identity=identity)
        check(result[-1] == "INIT 101", f"authenticate_cache_{name}")

    original_identity = identity.read_text().splitlines()
    identity_mutations = {
        "expected_payload": (0, "0" * 64, 101), "grid_size": (1, "6 4", 101),
        "source": (3, "a" * 64, 12), "preprocessing": (4, "b" * 64 + " 2 mean", 101),
        "radius": (4, PROCESS_SHA + " 1 mean", 105), "route": (4, PROCESS_SHA + " 2 first", 105),
        "grid_spacing": (5, "-12 0.31 11.5 0.06", 101), "floor": (6, "3.1 " + PATTERN_SHA, 105),
        "pattern_identity": (6, "3.2 " + "a" * 64, 105), "split": (2, f"0.7 0.02 {CRAD} 0", 105),
        "radiation": (2, f"0.7 0.02 {CRAD * 1.0000044} 1", 105),
        "hydrogen": (2, f"0.6 0.02 {CRAD} 1", 105), "metals": (2, f"0.7 0.03 {CRAD} 1", 105)}
    for name, (line_index, text, code) in identity_mutations.items():
        altered = original_identity.copy()
        altered[line_index] = text
        path = output / f"wrong_{name}.identity"
        path.write_text("\n".join(altered) + "\n")
        result = invoke("initialize", cache=cache, identity=path)
        check(result[-1] == f"INIT {code}", f"identity_guard_{name}")
    bad_pattern = output / "wrong_pattern.txt"
    bad_pattern.write_text(PATTERN.read_text().replace("0.171836", "0.171835"))
    result = invoke("initialize", cache=cache, identity=identity, pattern=bad_pattern)
    check(result[-1] == "INIT 101", "pattern_payload_authenticated")

    solar_cache = WORKSPACE / "freeeos_biquintic_progress/experiments/008_standalone_table/solar_split.cache"
    check(digest(solar_cache) == CACHE_SHA and digest(PATTERN) == PATTERN_SHA, "archived_payload_and_pattern_identity")
    solar_identity = EXPERIMENT / "solar.identity"
    solar_identity.write_text(identity_text(solar_cache, CACHE_SHA))
    solar, _ = query(solar_cache, solar_identity, [(0, -5.1, 4.51), (1, -5.1, 4.51)])
    check(np.all(solar[:, 0] == 0) and np.array_equal(solar[0], solar[1]), "existing_FreeEOS_cache_query_and_reordering")
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from freeeos_potential import evaluate, thermodynamic_partials

    with solar_cache.open() as stream:
        header = [next(stream).split() for _ in range(6)]
    nq, nt = map(int, header[1])
    umin, du, tmin, dt = map(float, header[5])
    rows = np.loadtxt(solar_cache, skiprows=6)
    nodes = np.empty((3, 3, nq, nt))
    for temp_order in range(3):
        for density_order in range(3):
            nodes[density_order, temp_order] = rows[:, 3 + density_order + 3 * temp_order].reshape(nt, nq).T
    complete = rows[:, 2].reshape(nt, nq).T.astype(bool)
    logq_axis = (umin + du * np.arange(nq)) / LN10
    logt_axis = (tmin + dt * np.arange(nt)) / LN10
    sample_t = np.array([4.51, 4.93, 5.51])
    sample_rho = np.array([-5.1, -3.31, -1.11])
    sample_q = sample_rho - 2 * sample_t + 12
    native = evaluate(nodes, complete, logq_axis, logt_axis, sample_q, sample_t, add_radiation=True)
    python_state = thermodynamic_partials(native, sample_q, sample_t)
    native_data, _ = query(solar_cache, solar_identity, [(0, density, temperature)
                                                      for density, temperature in zip(sample_rho, sample_t)])
    native_fields = {0: "Pgas", 1: "E", 2: "S", 6: "grad_ad", 7: "chiRho", 8: "chiT",
                     9: "Cp", 10: "Cv", 11: "E_rho", 12: "S_T", 13: "S_rho", 14: "gamma1", 15: "gamma3"}
    parity_errors = []
    for slot, name in native_fields.items():
        expected = python_state[name]
        if slot < 3:
            expected = np.log(expected)
        scale = np.maximum(np.abs(expected), 1e-6)
        parity_errors.extend(np.abs(native_data[:, 1 + slot] - expected) / scale)
        if slot in [11, 12, 13]:
            continue
        for axis, suffix in enumerate(["_x", "_t"]):
            expected = python_state[name + suffix]
            if slot < 3:
                expected = expected / python_state[name]
            scale = np.maximum.reduce([np.abs(expected), np.abs(native_data[:, 1 + slot]) * 1e-5,
                                       np.full(len(expected), 1e-6)])
            parity_errors.extend(np.abs(native_data[:, 28 + 27 * axis + slot] - expected) / scale)
    parity_error = float(max(parity_errors))
    check(np.all(native_data[:, 0] == 0) and parity_error < 2e-4, "archived_cache_independent_Python_basis_parity",
          states=3, max_scaled_error=parity_error)
    native_fd_states = [(0, -5.1, 4.51)]
    step = 1e-5
    for axis in range(2):
        for sign in [-1, 1]:
            native_fd_states.append((0, -5.1 + (axis == 0) * sign * step / LN10,
                                     4.51 + (axis == 1) * sign * step / LN10))
    native_fd, _ = query(solar_cache, solar_identity, native_fd_states)
    for axis in range(2):
        numerical = (native_fd[2 + 2 * axis, 1:28] - native_fd[1 + 2 * axis, 1:28]) / (2 * step)
        derivative = native_fd[0, 28 + axis * 27:55 + axis * 27]
        scale = np.maximum.reduce([np.abs(numerical), np.abs(derivative), np.abs(native_fd[0, 1:28]) * 1e-5])
        error = np.max(np.abs(numerical[THERMO] - derivative[THERMO]) / np.maximum(scale[THERMO], 1e-5))
        check(np.all(native_fd[:, 0] == 0) and error < 2e-4, f"archived_cache_all13_fd_partials_{axis}",
              max_scaled_error=float(error))
    thread_lines = invoke("threads", "-5.1 4.51\n", solar_cache, solar_identity)
    check(thread_lines[-1] == "THREADS 2 1000 0", "existing_FreeEOS_cache_two_thread_1000_queries")
    results = dict(passed=all(item["passed"] for item in checks), checks=checks, count=len(checks),
                   elapsed_seconds=time.time() - started, hostname=socket.gethostname(), job=os.environ.get("SLURM_JOB_ID"),
                   threads=2, memory_limit_GiB=8, driver=str(arguments.driver.resolve()), commands=commands,
                   driver_sha256=digest(arguments.driver), thread_states_per_cache=4,
                   self_maxrss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                   child_maxrss_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                   adapter_sources={path.name: digest(path) for path in sorted(Path(__file__).parent.glob("*.f90"))},
                   source_sha256=SOURCE_SHA, preprocessing_sha256=PROCESS_SHA, cache_sha256=CACHE_SHA,
                   pattern_sha256=PATTERN_SHA, scope="Adapter unit validation, not component/final API acceptance")
    (output / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(f"{sum(item['passed'] for item in checks)}/{len(checks)} checks passed", flush=True)
    if not results["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
