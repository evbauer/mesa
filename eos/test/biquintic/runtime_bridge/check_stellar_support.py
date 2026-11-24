"""Replay a logged legacy stellar trajectory through the candidate EOS only."""

import json

import numpy as np

from api_client import EOSClient, FIELDS
from runtime_environment import WORKSPACE, installed_inventory, sha256, verify_private_cache


def run():
    source = WORKSPACE / "freeeos_biquintic_progress/experiments/049_stellar_preparation/control_legacy/queries.tsv"
    output = WORKSPACE / "freeeos_biquintic_progress/experiments/052_stellar_api/support_preflight"
    unique = {}
    row_count = 0
    for record in source.read_text().splitlines():
        if record.startswith("#") or not record.strip():
            continue
        columns = record.split()
        point = (float(columns[4]), float(columns[5]))
        fractions = np.asarray(columns[6:14], dtype=float)
        if point in unique and not np.allclose(unique[point]["fractions"], fractions, rtol=0, atol=2e-12):
            raise ValueError("Legacy trace reports inconsistent routing at identical state")
        unique.setdefault(point, {"fractions": fractions, "phases": set()})["phases"].add(columns[3])
        row_count += 1
    before = installed_inventory()
    states = []
    with EOSClient(output) as client:
        for point, previous in unique.items():
            values = client.query("candidate", "final", *point)
            status = int(values[1])
            state = {"logRho": point[0], "logT": point[1], "status": status,
                     "legacy_fractions": previous["fractions"].tolist(), "phases": sorted(previous["phases"])}
            if status == 0:
                fractions = values[4 + 19:4 + 27]
                fraction = values[4 + FIELDS["frac_FreeEOS"]]
                state.update(fractions=fractions.tolist(), callback_count=int(values[2]),
                    hybrid_retained=int(values[3]), finite=bool(np.all(np.isfinite(values[4:]))),
                    fraction_difference=float(np.max(abs(fractions - previous["fractions"]))))
                if values[2] != int(fraction > 0):
                    raise ValueError("Unexplained callback count on stellar replay")
            states.append(state)
    after = installed_inventory()
    verify_private_cache(output)
    failures = [state for state in states if state["status"] != 0]
    successful = [state for state in states if state["status"] == 0]
    summary = {"scope": "EOS replay of legacy stellar states, NOT candidate stellar evolution; fixed audited composition client",
        "source": str(source.relative_to(WORKSPACE)), "source_sha256": sha256(source),
        "logged_rows": row_count, "unique_states": len(unique), "failures": failures,
        "min_logT": min(point[1] for point in unique), "max_logT": max(point[1] for point in unique),
        "min_FreeEOS_fraction": min(state["fractions"][2] for state in successful),
        "max_FreeEOS_fraction": max(state["fractions"][2] for state in successful),
        "maximum_fraction_difference": max(state["fraction_difference"] for state in successful),
        "hybrid_and_finite": all(state["hybrid_retained"] == 1 and state["finite"] for state in successful),
        "installed_unchanged": before == after,
        "limitation": "All-state actual stellar composition validation belongs to049; altered candidate trajectories are not covered"}
    (output / "states.json").write_text(json.dumps(states, indent=2) + "\n")
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if failures or not summary["hybrid_and_finite"] or before != after or summary["maximum_fraction_difference"] > 2e-12:
        raise AssertionError("Stellar-support EOS replay gate failed")


if __name__ == "__main__":
    run()
