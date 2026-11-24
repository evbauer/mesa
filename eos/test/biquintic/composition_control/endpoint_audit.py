"""Independent printed-grid evidence for a narrowly authorized metadata resolution."""

from pathlib import Path

import numpy as np


HEADER_FIELDS = ("version", "X", "Z", "num_logT", "logT_min", "logT_max", "delta_logT",
                 "num_logQ", "logQ_min", "logQ_max", "delta_logQ")
ENDPOINT_FIELDS = ("logQ_max", "logT_max")
COORDINATE_TOLERANCE = 1e-10


def comparison_header(header, ignore_endpoint_metadata=False):
    ignored = {"X", "Z"}
    if ignore_endpoint_metadata:
        ignored.update(ENDPOINT_FIELDS)
    return np.asarray([value for name, value in zip(HEADER_FIELDS, header) if name not in ignored])


def header_grid(header):
    if header.shape != (11,) or not np.all(np.isfinite(header)) or header[0] != 51:
        raise ValueError("Audited resolution requires a finite format-51 header")
    for index in (3, 7):
        if header[index] != int(header[index]) or header[index] < 2:
            raise ValueError("Grid counts must be exact positive integers, not endpoint overrides")
    if header[6] <= 0 or header[10] <= 0:
        raise ValueError("Grid spacings must be positive")
    temperature_count, density_count = int(header[3]), int(header[7])
    derived_q = header[8] + np.arange(density_count) * header[10]
    derived_t = header[4] + np.arange(temperature_count) * header[6]
    return derived_q, derived_t


def read_endpoint_table(path, allow_audited_endpoint_metadata=False):
    with Path(path).open() as stream:
        next(stream)
        header = np.fromstring(next(stream), sep=" ")
        logq, logt = header_grid(header)
        if allow_audited_endpoint_metadata and not np.isclose(header[2], .02, rtol=0, atol=1e-12):
            raise ValueError("Independent endpoint override is authorized only at Z=.02")
        records = []
        for line in stream:
            columns = line.split()
            if len(columns) == 22 and columns[0][0] in "+-0123456789":
                records.append(np.asarray(columns, dtype=float))
    expected_records = len(logq) * len(logt)
    if len(records) != expected_records:
        raise ValueError("Source record count differs from declared num_logQ*num_logT; endpoint overrides cannot change counts")
    if not allow_audited_endpoint_metadata:
        for name, declared, derived in (("logQ_max", header[9], logq[-1]), ("logT_max", header[5], logt[-1])):
            if not np.isclose(declared, derived, rtol=0, atol=COORDINATE_TOLERANCE):
                raise ValueError(f"Strict independent loader rejects inconsistent {name}; require explicit audited endpoint option")
    rows = np.asarray(records).reshape(len(logq), len(logt), 22)
    with Path(path).open() as stream:
        audit = audit_printed_grid(stream, header, logq, logt)
    audit["loader"] = "independent sidecar strict/explicit-audited loader; no core read_table call"
    audit["allow_audited_endpoint_metadata"] = allow_audited_endpoint_metadata
    return header, logq, logt, rows, audit


def audit_printed_grid(stream, header, logq, logt):
    derived_q, derived_t = header_grid(header)
    density_count, temperature_count = len(derived_q), len(derived_t)
    if not np.array_equal(logq, derived_q) or not np.array_equal(logt, derived_t):
        raise ValueError("Actual grid/count/min/step differs; endpoint metadata cannot authorize resampling")
    next(stream)
    printed_header = np.fromstring(next(stream), sep=" ")
    if not np.array_equal(printed_header, header):
        raise ValueError("Independently parsed header differs from the table loader")
    labels, block_counts = [], []
    expect_label = False
    maximum_temperature_error = 0.0
    maximum_density_error = 0.0
    known_density_errors = []
    records = 0
    for line in stream:
        stripped = line.strip()
        if stripped == "logQ = logRho - 2*logT + 12":
            if expect_label:
                raise ValueError("Missing Q-block label")
            expect_label = True
            continue
        if expect_label:
            if not stripped:
                continue
            if len(stripped.split()) != 1:
                raise ValueError("Q-block label must be a scalar")
            label = float(stripped)
            density_index = len(labels)
            if density_index >= density_count:
                raise ValueError("More printed Q blocks than declared count")
            if not np.isfinite(label) or abs(label - derived_q[density_index]) > COORDINATE_TOLERANCE:
                raise ValueError(f"Printed Q-block label differs from count/min/step at iq={density_index}")
            if labels and label <= labels[-1]:
                raise ValueError("Q-block labels are not strictly monotonic")
            labels.append(label)
            block_counts.append(0)
            expect_label = False
            continue
        columns = stripped.split()
        if len(columns) != 22 or columns[0][0] not in "+-0123456789":
            continue
        if not labels:
            raise ValueError("Printed records appear before a Q-block label")
        density_index, temperature_index = len(labels) - 1, block_counts[-1]
        if temperature_index >= temperature_count:
            raise ValueError("Printed Q-block record count exceeds declared num_logT")
        printed_temperature, printed_density = float(columns[0]), float(columns[18])
        if not np.isfinite(printed_temperature) or not np.isfinite(printed_density):
            raise ValueError("Nonfinite printed coordinate")
        temperature_error = abs(printed_temperature - derived_t[temperature_index])
        if temperature_error > COORDINATE_TOLERANCE:
            raise ValueError(f"Printed T differs from the authoritative grid at iq={density_index}, it={temperature_index}")
        density_error = printed_density - (derived_q[density_index] + 2 * derived_t[temperature_index] - 12)
        if abs(density_error) > COORDINATE_TOLERANCE:
            known_exclusion = (np.allclose(header[1:3], (.7, .02), rtol=0, atol=1e-12)
                               and density_index == 520 and abs(density_error + 1e-4) <= COORDINATE_TOLERANCE)
            if not known_exclusion:
                raise ValueError(f"Printed rho differs outside the known X=.70/Z=.02 iq=520 exclusion at iq={density_index}, it={temperature_index}")
            known_density_errors.append(density_error)
        else:
            maximum_density_error = max(maximum_density_error, abs(density_error))
        maximum_temperature_error = max(maximum_temperature_error, temperature_error)
        block_counts[-1] += 1
        records += 1
    if expect_label or len(labels) != density_count or any(count != temperature_count for count in block_counts):
        raise ValueError("Printed Q-block label/record count differs from declared num_logQ*num_logT")
    if records != density_count * temperature_count:
        raise ValueError("Printed record count differs from declared num_logQ*num_logT")
    if known_density_errors and len(known_density_errors) != temperature_count:
        raise ValueError("Known iq=520 density discrepancy must cover exactly its already-excluded block")
    return {"status": "PASSED", "independent_Q_block_parser": True,
            "printed_record_count": records, "declared_record_count": density_count * temperature_count,
            "Q_block_count": len(labels), "records_per_block": temperature_count,
            "Q_labels_strictly_monotonic": True, "Q_labels_match_count_min_step": True,
            "Q_label_min": labels[0], "Q_label_max": labels[-1],
            "derived_logQ_endpoint": float(derived_q[-1]), "derived_logT_endpoint": float(derived_t[-1]),
            "printed_temperature_max_grid_error": maximum_temperature_error,
            "printed_density_max_grid_error_outside_known_exclusion": maximum_density_error,
            "known_coordinate_exclusions": [{"X": .7, "Z": .02, "iq": 520,
                                              "record_count": len(known_density_errors),
                                              "signed_density_error_min": min(known_density_errors),
                                              "signed_density_error_max": max(known_density_errors)}]
                                           if known_density_errors else [],
            "authoritative_grid": "Declared count/min/step independently confirmed by every Q label, block count and printed rho/T",
            "display_endpoint_policy": "Only logQ_max/logT_max may be ignored; no source edits, count changes, coordinate repair or resampling"}
