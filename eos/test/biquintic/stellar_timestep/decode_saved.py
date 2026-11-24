"""Lead-cleared read-only correction of the saved half-baseline profile index."""

import time

import timestep as task
from timestep_audit import audit_run


task.frozen.placement()
if time.time() >= task.NUMERICAL_STOP:
    raise RuntimeError("past numerical-analysis cutoff")
if task.frozen.digest(task.EXPERIMENT / "CONTINUATION_STOP_half_baseline.json") != "36ec3412d58530e412ac3f8f48fd5e5aa8453919e347bd68bfdaef3257791f08":
    raise ValueError("original decoder failure record changed")
task.verify_frozen()
audit_run("half_baseline")
