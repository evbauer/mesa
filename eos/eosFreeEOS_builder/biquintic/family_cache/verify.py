"""Independently authenticate a completed standalone family and its sources."""

import argparse
import json
from pathlib import Path

from produce import SOURCE_DIR, DEFAULT_OUTPUT, checksum, source_identity


def verify(source_dir, directory):
    sources = {source.stem: source for source in source_dir.glob("mesa-FreeEOS_*.data")}
    if len(sources) != 89:
        raise ValueError(f"Expected 89 distinct source tables, found {len(sources)}")
    lines = (directory / "family.manifest").read_text().splitlines()
    if lines[:2] != ["FREEEOS_BIQUINTIC_FAMILY_V1", "89"] or len(lines) != 91:
        raise ValueError("Manifest magic, count or row count is wrong")
    audit = json.loads((directory / "family.audit.json").read_text())
    if audit["entry_count"] != 89 or len(audit["tables"]) != 89:
        raise ValueError("Family audit entry count is wrong")
    seen = set()
    available = 0
    for line in lines[2:]:
        columns = line.split()
        if len(columns) != 6:
            raise ValueError("Manifest row must have exactly six fields")
        hydrogen, metals, flag, filename, source_sha, payload_sha = columns
        if flag not in ("0", "1"):
            raise ValueError("Invalid availability flag")
        if flag == "0":
            if filename != "-" or payload_sha != "-":
                raise ValueError("Unavailable row violates contract")
            matches = [source for source in sources.values()
                       if source_identity(source)[:2] == (float(hydrogen), float(metals))]
            if len(matches) != 1:
                raise ValueError("Unavailable row cannot be identified uniquely")
            source = matches[0]
        else:
            if not filename.endswith(".bq") or Path(filename).name != filename:
                raise ValueError("Payload filename is not a local .bq file")
            source = sources.get(filename[:-3])
            if source is None:
                raise ValueError(f"Payload has no source: {filename}")
            payload = directory / filename
            if not payload.is_file() or checksum(payload) != payload_sha:
                raise ValueError(f"Payload SHA mismatch: {filename}")
            with payload.open() as stream:
                header = [next(stream).split() for _ in range(6)]
            if (header[0] != ["FREEEOS_BIQUINTIC_V1"] or len(header[1]) != 2 or
                    header[2][3] != "1" or header[3] != [source_sha] or
                    header[4][1:] != ["2", "mean"]):
                raise ValueError(f"Payload contract mismatch: {filename}")
            if abs(float(header[2][0]) - float(hydrogen)) > 1e-12 or abs(float(header[2][1]) - float(metals)) > 1e-12:
                raise ValueError(f"Payload composition mismatch: {filename}")
            available += 1
        if source.stem in seen:
            raise ValueError(f"Source appears more than once: {source.name}")
        seen.add(source.stem)
        if source_identity(source)[:2] != (float(hydrogen), float(metals)):
            raise ValueError(f"Manifest composition differs from filename: {source.name}")
        if checksum(source) != source_sha:
            raise ValueError(f"Source SHA mismatch: {source.name}")
        per_table = json.loads((directory / (source.stem + ".audit.json")).read_text())
        if per_table["source_sha256"] != source_sha or per_table["available"] != int(flag):
            raise ValueError(f"Per-source audit mismatch: {source.name}")
    if seen != set(sources) or available != audit["available"]:
        raise ValueError("Family coverage or availability disagrees with audit")
    print(f"PASS {len(seen)} source SHA checks, {available} payload SHA/header checks, "
          f"{89 - available} unavailable entries")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=SOURCE_DIR)
    parser.add_argument("--directory", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    verify(arguments.source_dir, arguments.directory)


if __name__ == "__main__":
    main()
