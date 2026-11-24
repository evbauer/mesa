"""Matched component/final identity checks at saved stellar-profile states."""

import argparse
import json
import os
from pathlib import Path
import shlex
import sys

import numpy as np


CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE.parent / "runtime_bridge"))
from api_client import EOSClient, FIELDS, LN10, CRAD
from check_blends_faces import balances
from runtime_environment import WORKSPACE, installed_inventory, sha256, verify_private_cache


LABELS = ("Dpe", "Dse", "Dsp", "energy_Cv")


def profile(path):
    lines = path.read_text().splitlines()
    metadata = dict(zip(lines[1].split(), shlex.split(lines[2])))
    names = lines[5].split()
    values = np.loadtxt(path, skiprows=6, ndmin=2)
    if values.shape[1] != len(names) or not np.all(np.isfinite(values)):
        raise ValueError("Malformed or nonfinite stellar profile")
    return metadata, {name: values[:, column] for column, name in enumerate(names)}


def diagnostics(result, point):
    values, density, thermal = result[4:31], result[31:58], result[58:85]
    rho, temperature = np.exp(np.asarray(point) * LN10)
    gas, energy, entropy = np.exp(values[:3])
    pressure = gas + CRAD * temperature**4 / 3
    pressure_t = gas * thermal[0] + 4 * CRAD * temperature**4 / 3
    terms = ((rho * energy * density[1], pressure_t, -pressure),
             (temperature * entropy * thermal[2], -energy * thermal[1]),
             (-rho * temperature * entropy * density[2], -pressure_t),
             (energy * thermal[1], -temperature * values[FIELDS["Cv"]]))
    scaled = [sum(group) / max(sum(abs(value) for value in group), 1.) for group in terms]
    ratios = balances(result, point)
    if not np.all(np.isfinite(ratios + scaled)):
        raise ValueError("Nonfinite identity diagnostic")
    return {"ratios": dict(zip(LABELS, ratios)), "scaled_balances": dict(zip(LABELS, scaled))}


def distribution(records):
    result = {}
    for key in ("off_component", "off_final", "candidate_component", "candidate_final"):
        result[key] = {}
        for population, selected in (
                ("all", records),
                ("pure", [record for record in records if record["fraction_FreeEOS"] == 1.]),
                ("blend", [record for record in records if 0 < record["fraction_FreeEOS"] < 1.])):
            summaries = {"count": len(selected)}
            if selected:
                for route in ("ratios", "scaled_balances"):
                    summaries[route] = {name: dict(zip(("median_abs", "p99_abs", "max_abs"),
                        np.quantile([abs(record[key][route][name]) for record in selected], [.5, .99, 1.]).tolist()))
                        for name in LABELS}
            result[key][population] = summaries
    return result


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if (output / "summary.json").exists():
        raise ValueError("Use a new output directory; retain the earlier audit")
    before = installed_inventory()
    inputs = {str(path.resolve().relative_to(WORKSPACE)): sha256(path) for path in args.profiles}
    all_records, profiles = [], []
    with EOSClient(output / "api") as client:
        for path in args.profiles:
            metadata, columns = profile(path)
            records = []
            for index, point in enumerate(zip(columns["logRho"], columns["logT"])):
                returns = {(model, route): client.query(model, route, *point)
                           for model in ("off", "candidate") for route in ("component", "final")}
                if any(values[1] != 0 for values in returns.values()):
                    raise ValueError(f"EOS failure at {path}, zone {columns['zone'][index]}")
                fractions = [returns[model, "final"][4 + FIELDS["frac_FreeEOS"]]
                             for model in ("off", "candidate")]
                if fractions[0] != fractions[1] or abs(fractions[1] - columns["eos_frac_FreeEOS"][index]) > 2e-12:
                    raise ValueError("Profile/replay or matched-backend routing disagreement")
                record = {"zone": int(columns["zone"][index]), "logRho": float(point[0]),
                          "logT": float(point[1]), "fraction_FreeEOS": float(fractions[1])}
                for (model, route), values in returns.items():
                    record[model + "_" + route] = diagnostics(values, point)
                records.append(record)
            profiles.append({"input": str(path.resolve().relative_to(WORKSPACE)), "metadata": metadata,
                             "distributions": distribution(records), "records": records})
            all_records.extend(records)
        query_count = client.count
    verify_private_cache(output / "api")
    after = installed_inventory()
    if before != after or any(sha256(WORKSPACE / path) != digest for path, digest in inputs.items()):
        raise ValueError("Protected inputs changed")
    result = {"scope": "EOS-only matched-state replay of saved candidate stellar profiles, not another evolution run",
              "sampling": "all saved zones of explicitly named profiles; unweighted, not population or mass-weighted statistics",
              "models": "both backends and both APIs at identical candidate-profile rho/T; no temporal or spatial EOS-model matching",
              "metrics": "returned-function identities in order Dpe,Dse,Dsp and e_t/(T Cv)-1; total radiation convention",
              "scaled_balance": "signed identity numerator divided by max(1,sum(abs(additive terms)))",
              "command": sys.argv, "host": os.uname().nodename, "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
              "inputs": inputs, "driver_sha256": sha256(CODE.parents[4] / "freeeos_biquintic_progress/experiments/050_runtime_api/build/driver/candidate_driver"),
              "code_sha256": {str(path.relative_to(WORKSPACE)): sha256(path) for path in
                  (Path(__file__), CODE.parent / "runtime_bridge/api_client.py", CODE.parent / "runtime_bridge/check_blends_faces.py")},
              "installed_unchanged": before == after, "actual_queries": query_count,
              "profile_count": len(profiles), "zone_records": len(all_records),
              "distributions": distribution(all_records), "profiles": profiles}
    (output / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({name: result[name] for name in ("profile_count", "zone_records", "actual_queries", "installed_unchanged")}, indent=2))


if __name__ == "__main__":
    run()
