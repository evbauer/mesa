"""Run unchanged051 consumer recipes against a separately linked053 variant."""

from datetime import datetime
import importlib
import json
from pathlib import Path
import sys

DIRECTORY = Path(__file__).resolve().parent
CODE = DIRECTORY.parent
WORKSPACE = DIRECTORY.parents[4]
OUTPUT = WORKSPACE / "freeeos_biquintic_progress/experiments/053_arithmetic_control"
sys.path.insert(0, str(CODE))
from runtime_consumers import run_consumers as runner
import api_client
from runtime_environment import sha256


def main():
    output = OUTPUT / "consumer_variant"
    if (output / "full/trials.jsonl").exists():
        raise RuntimeError("Preserve existing053 trials; use a separately named experiment for a rerun")
    original_output = runner.OUTPUT
    api_client.EXPERIMENT = OUTPUT
    runner.OUTPUT = output
    runner.main()
    sys.modules["run_consumers"] = runner
    analysis = importlib.import_module("runtime_consumers.analyze")
    analysis.main()
    runner.write_json(output / "variant_identity.json", dict(
        completed=datetime.now().isoformat(), original051_output=str(original_output),
        control="Separate arithmetic-only real native-cell evaluator; unchanged nodes, metadata, coordinate mapping, adapter, EOS and solver",
        baseline_driver_sha256="62cbeea83a1fa3bc833a811e0f0586e4ab9fbd8f10e1bb15b6dad2177d287b16",
        source_sha256={str(path.relative_to(WORKSPACE)): sha256(path) for path in
                      (WORKSPACE / "mesa/eos/private/arithmetic_weights.f90", DIRECTORY / "freeeos_biquintic.f90", Path(__file__))},
        scientific_baseline_rerun=False, unchanged051_failures=True))
    prior = {tuple(record[key] for key in ("consumer", "route", "trial_id")): record for record in
             (json.loads(line) for line in (original_output / "full/trials.jsonl").read_text().splitlines())}
    differences = []
    for line in (output / "full/trials.jsonl").read_text().splitlines():
        record = json.loads(line)
        original = prior[tuple(record[key] for key in ("consumer", "route", "trial_id"))]
        before, after = original["newton"], record["newton"]
        if original["status"] != record["status"] or before["iterations"] != after["iterations"]:
            differences.append(dict(consumer=record["consumer"], model=record["model"], route=record["route"],
                                    trial_id=record["trial_id"], case_id=record["case_id"], family=record["family"],
                                    original_status=original["status"], variant_status=record["status"],
                                    original_iterations=before["iterations"], variant_iterations=after["iterations"]))
    runner.write_json(output / "differences_from051.json", differences)
    print(json.dumps(differences, indent=2))


if __name__ == "__main__":
    main()
