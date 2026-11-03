"""Shared definitions for the 070 native-face comparison."""

import hashlib
import json
import os
import socket
from pathlib import Path

import numpy as np


CODE = Path(__file__).resolve().parent
BIQUINTIC = CODE.parent
WORKSPACE = BIQUINTIC.parents[3]
OUTPUT = WORKSPACE / "freeeos_biquintic_progress/experiments/070_freeeos_helm_map"
OLD = WORKSPACE / "freeeos_biquintic_progress/experiments/060_pure_freeeos_faces"
HELM048 = BIQUINTIC / "helm_control"
SOURCE = WORKSPACE / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
MASSFILE = WORKSPACE / "mesa/eos/eosFreeEOS_builder/mass_frac_solar.txt"
JOB = "5549798"
NODE = "rzwhippet25"
DISPLAY_FLOOR = 1e-17


def placement():
    if socket.gethostname().split(".")[0] != NODE or os.environ.get("SLURM_JOB_ID") != JOB:
        raise RuntimeError("TASK070 requires the inherited rzwhippet25/job5549798 allocation")


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def scaled_jump(response_left, response_right, partial_left, partial_right):
    """Unit-invariant signed response-scaled jump; all-zero inputs give zero."""
    arrays = np.broadcast_arrays(*[np.asarray(value) for value in
                                   (response_left, response_right, partial_left, partial_right)])
    if not all(np.isfinite(value).all() for value in arrays):
        raise ValueError("nonfinite one-sided trace")
    denominator = np.maximum.reduce([np.abs(value) for value in arrays])
    difference = arrays[2] - arrays[3]
    return np.divide(difference, denominator, out=np.zeros_like(difference, dtype=float),
                     where=denominator != 0), denominator


def summary(values):
    values = np.asarray(values, dtype=float)
    if not len(values):
        return {"count": 0}
    finite = values[np.isfinite(values)]
    if not len(finite):
        return {"count": len(values), "finite": 0}
    magnitude = np.abs(finite)
    return {"count": int(len(values)), "finite": int(len(finite)),
            "zero": int(np.count_nonzero(magnitude == 0)),
            "median": float(np.median(magnitude)), "p90": float(np.quantile(magnitude, .9)),
            "p99": float(np.quantile(magnitude, .99)), "maximum": float(np.max(magnitude)),
            "below_display": int(np.count_nonzero(magnitude < 1e-5)),
            "above_display": int(np.count_nonzero(magnitude > 10))}
