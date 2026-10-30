"""Serial, logged actual-EOS client for fixed-composition experiments."""

import json
import os
import selectors
import subprocess
import time
from pathlib import Path

import numpy as np

from runtime_environment import MESA, WORKSPACE, EXPERIMENT, INLIST, MASSFILE, private_environment, sha256


CACHE = WORKSPACE / "freeeos_biquintic_progress/experiments/008_standalone_table/solar_split.cache"
IDENTITY = WORKSPACE / "freeeos_biquintic_progress/experiments/047_runtime_adapter/solar.identity"
PAYLOAD_SHA = "e2c72c6f02cb900649a6a629bfd4fc53c2e8423e5bbe6afd6a189c45cac3a204"
MANIFEST_SHA = "afc6fe676ba13bf858fc0a611268402b3663f4ad1e0ff991049713fd78788f19"
LN10 = np.log(10.)
CRAD = 7.565733250280004e-15
FIELDS = {"lnPgas": 0, "lnE": 1, "lnS": 2, "mu": 3, "lnfree_e": 4, "eta": 5,
          "grad_ad": 6, "chiRho": 7, "chiT": 8, "Cp": 9, "Cv": 10,
          "dE_dRho": 11, "dS_dT": 12, "dS_dRho": 13, "gamma1": 14, "gamma3": 15,
          "frac_FreeEOS": 21}


def verify_identity():
    metadata = json.loads(CACHE.with_suffix(".cache.json").read_text())
    if metadata["cache_sha256"] != PAYLOAD_SHA or sha256(CACHE) != PAYLOAD_SHA:
        raise ValueError("Archived potential payload differs from the frozen experiment008 identity")
    if IDENTITY.read_text().split()[0] != PAYLOAD_SHA:
        raise ValueError("Runtime identity differs from the independently archived payload digest")


def validate_success(scenario, route, result):
    if result.shape != (125,) or result[1] != 0:
        raise ValueError("Expected a successful actual-EOS protocol record")
    if not np.all(np.isfinite(result)) or result[3] != 1:
        raise ValueError("Nonfinite actual EOS output or violated hybrid retention")
    fractions = result[23:31]
    if np.any(fractions < 0) or np.any(fractions > 1) or abs(np.sum(fractions) - 1) > 1e-12:
        raise ValueError("Invalid EOS component fractions")
    fraction = result[4 + FIELDS["frac_FreeEOS"]]
    if route == "component" and fraction != 1:
        raise ValueError("Forced FreeEOS component did not retain its identity")
    expected_calls = 0 if scenario == "off" else int(fraction > 0)
    if result[2] != expected_calls:
        raise ValueError("EOS callback count disagrees with returned routing")


class EOSClient:
    def __init__(self, output, timeout=60.):
        if not np.isfinite(timeout) or timeout <= 0:
            raise ValueError("A positive finite per-query timeout is required")
        self.timeout = float(timeout)
        verify_identity()
        self.output = Path(output).resolve()
        self.output.mkdir(parents=True, exist_ok=True)
        self.environment = private_environment(self.output)
        self.query_log = (self.output / "queries.txt").open("w")
        self.result_log = (self.output / "results.txt").open("w")
        self.error_log = (self.output / "driver.stderr").open("w")
        self.command = [str(EXPERIMENT / "build/driver/candidate_driver"), str(MESA),
                        str(INLIST), str(MASSFILE), str(CACHE), str(IDENTITY)]
        self.process = subprocess.Popen(self.command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                        stderr=self.error_log, text=True, bufsize=1,
                                        cwd=self.output, env=self.environment)
        self.count = 0
        self.memo = {}
        self.memo_hits = 0

    def query(self, scenario, route, log_rho, log_temp):
        key = (scenario, route, float(log_rho), float(log_temp))
        if key in self.memo:
            self.memo_hits += 1
            return self.memo[key].copy()
        record = f"{scenario} {route} {log_rho:.17g} {log_temp:.17g}\n"
        self.query_log.write(record)
        self.process.stdin.write(record)
        self.process.stdin.flush()
        response = self.read_response()
        self.result_log.write(response)
        values = np.fromstring(response, sep=" ")
        self.count += 1
        if len(values) != 125 or values[0] != self.count:
            raise ValueError("Actual-EOS protocol shape or sequence mismatch")
        if values[1] == 0:
            validate_success(scenario, route, values)
        self.memo[key] = values
        return values.copy()

    def read_response(self):
        deadline = time.monotonic() + self.timeout
        response = bytearray()
        with selectors.DefaultSelector() as selector:
            selector.register(self.process.stdout, selectors.EVENT_READ)
            while not response.endswith(b"\n"):
                remaining = deadline - time.monotonic()
                if remaining <= 0 or not selector.select(remaining):
                    raise TimeoutError("Actual-EOS query exceeded its response deadline")
                chunk = os.read(self.process.stdout.fileno(), 4096)
                if not chunk:
                    raise RuntimeError("EOS driver closed its response stream")
                response.extend(chunk)
                if len(response) > 65536 or response.count(b"\n") > 1:
                    raise ValueError("Unexpected actual-EOS response framing")
        return response.decode("ascii")

    def state(self, model, route, density_log, temperature_log):
        scenario = {"legacy": "off", "candidate": "candidate"}[model]
        result = self.query(scenario, route, density_log / LN10, temperature_log / LN10)
        if result[1] != 0:
            raise ValueError(f"Actual EOS status {int(result[1])} at {density_log}, {temperature_log}")
        values, density, thermal = result[4:31], result[31:58], result[58:85]
        validate_success(scenario, route, result)
        try:
            with np.errstate(over="raise", invalid="raise", divide="raise"):
                temperature = np.exp(temperature_log)
                gas = np.exp(values[FIELDS["lnPgas"]])
                radiation = CRAD * temperature**4 / 3
                state = {"T": temperature, "rho": np.exp(density_log), "Pgas": gas,
                         "P": gas + radiation, "P_x": gas * density[0],
                         "P_t": gas * thermal[0] + 4 * radiation,
                         "frac_FreeEOS": values[FIELDS["frac_FreeEOS"]], "callback_count": result[2]}
        except FloatingPointError as error:
            raise ValueError("Nonfinite reconstructed actual-EOS state") from error
        for field in ("Cv", "Cp", "chiRho", "chiT", "grad_ad"):
            column = FIELDS[field]
            state[field] = values[column]
            state[field + "_x"] = density[column]
            state[field + "_t"] = thermal[column]
        if not all(np.isfinite(value) for value in state.values()):
            raise ValueError("Nonfinite reconstructed actual-EOS state or partial")
        if min(state["P"], state["Cv"], state["Cp"], state["chiRho"]) <= 0:
            raise ValueError("Nonpositive actual EOS pressure or response")
        return state

    def close(self):
        if self.process.stdin and not self.process.stdin.closed:
            self.process.stdin.close()
        try:
            status = self.process.wait(timeout=30)
        except subprocess.TimeoutExpired:
            self.process.terminate()
            try:
                self.process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                self.process.kill()
                self.process.wait(timeout=5)
            raise
        finally:
            self.process.stdout.close()
            self.query_log.close()
            self.result_log.close()
            self.error_log.close()
        if status:
            raise RuntimeError(f"EOS driver exited with status {status}")
        (self.output / "client.json").write_text(json.dumps({"command": self.command,
            "actual_queries": self.count, "identical_query_memo_hits": self.memo_hits,
            "per_query_timeout_seconds": self.timeout,
            "cache_sha256": PAYLOAD_SHA, "identity_sha256": sha256(IDENTITY),
            "executable_sha256": sha256(self.command[0]),
            "semantics": "serial logged experiment driver; memoization only for exactly identical queries"}, indent=2) + "\n")

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception, traceback):
        self.close()
