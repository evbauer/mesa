"""Export an isolated, metadata-bearing standalone cache; never a MESA cache."""

import argparse
import json
from pathlib import Path

import numpy as np

from freeeos_potential import CRAD, LN10, checksum, construct_twists, read_table, source_constraints


def export_cache(source, output, radius=2, route="mean", split_radiation=True):
    header, logq, logt, rows = read_table(source)
    low, valid, _ = source_constraints(rows, logq, logt, split_radiation=split_radiation)
    nodes, complete, _ = construct_twists(low, valid, (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10, radius, route)
    metadata = {"format": "FREEEOS_BIQUINTIC_V1", "source_sha256": checksum(source),
                "X": float(header[1]), "Z": float(header[2]), "radiation_constant": float(CRAD),
                "split_radiation": split_radiation, "radius": radius, "route": route,
                "preprocessing_sha256": checksum(Path(__file__).with_name("freeeos_potential.py")),
                "minimum_logT": 3.2, "valid_nodes": int(complete.sum()),
                "limitations": "strict standalone prototype; complete native stencils only; no legacy fallback; not production source accuracy"}
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as stream:
        stream.write("FREEEOS_BIQUINTIC_V1\n")
        stream.write(f"{len(logq)} {len(logt)}\n")
        stream.write(f"{header[1]:.17e} {header[2]:.17e} {CRAD:.17e} {int(split_radiation)}\n")
        stream.write(metadata["source_sha256"] + "\n")
        stream.write(f"{metadata['preprocessing_sha256']} {radius} {route}\n")
        stream.write(" ".join(f"{value:.17e}" for value in (logq[0] * LN10, (logq[1] - logq[0]) * LN10,
                                                            logt[0] * LN10, (logt[1] - logt[0]) * LN10)) + "\n")
        for temperature_index in range(len(logt)):
            for density_index in range(len(logq)):
                accepted = complete[density_index, temperature_index]
                corner = nodes[:, :, density_index, temperature_index].ravel(order="F") if accepted else np.zeros(9)
                stream.write(f"{density_index + 1} {temperature_index + 1} {int(accepted)} "
                             + " ".join(f"{value:.17e}" for value in corner) + "\n")
    metadata["cache_sha256"] = checksum(output)
    output.with_suffix(output.suffix + ".json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def main():
    workspace = Path(__file__).resolve().parents[4]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=workspace / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--radius", type=int, choices=(1, 2, 3), default=2)
    parser.add_argument("--route", choices=("mean", "first", "second"), default="mean")
    parser.add_argument("--total-potential", action="store_true")
    arguments = parser.parse_args()
    print(json.dumps(export_cache(arguments.source, arguments.output, arguments.radius,
                                  arguments.route, not arguments.total_potential), indent=2))


if __name__ == "__main__":
    main()
