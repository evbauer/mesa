"""Archived, no-refit active FreeEOS slice with one measured MESA selection stream."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

sys.dont_write_bytecode = True
CODE = Path(__file__).resolve().parents[1]
WORKSPACE = CODE.parents[3]
OUTPUT = WORKSPACE / "freeeos_biquintic_progress/experiments/023_active_slice"
ARCHIVE = WORKSPACE / "freeeos_biquintic_progress/experiments/015_smooth_validation/active_chiT"
sys.path.insert(0, str(CODE))
sys.path.insert(0, str(WORKSPACE / "freeeos_report"))

import numpy as np

import audit_freeeos as baseline
from freeeos_potential import (
    CRAD, LN10, checksum, construct_twists, evaluate, identity_residuals, read_table,
    source_constraints, thermodynamic_partials, valid_cells, query_cells,
)
from smooth_patch import evaluate_fit


SLICE_LOGT = 6.756476983450105
CANDIDATES = ("thermo_0.01", "thermo_0.1")
REGION_NAMES = ("Pure FreeEOS", "FreeEOS/Skye blend", "FreeEOS inactive")
FRACTION_NAMES = ("HELM", "OPAL_SCVH", "FreeEOS", "PC", "Skye", "CMS", "ideal", "WDM_AV")


def slice_queries():
    sample_q = np.linspace(.50, 1.10, 1001)
    sample_t = np.full(sample_q.shape, SLICE_LOGT)
    return sample_q, sample_t


def fraction_regions(freeeos):
    if not np.all(np.isfinite(freeeos)) or np.any((freeeos < -1e-10) | (freeeos > 1 + 1e-10)):
        raise ValueError("Invalid measured FreeEOS fraction")
    return np.where(freeeos <= 1e-12, 2, np.where(freeeos >= 1 - 1e-12, 0, 1))


def region_spans(sample_q, regions):
    edges = np.r_[sample_q[0], (sample_q[:-1] + sample_q[1:]) / 2, sample_q[-1]]
    changes = np.r_[0, np.flatnonzero(np.diff(regions)) + 1, len(regions)]
    return [{"label": REGION_NAMES[int(regions[start])], "code": int(regions[start]),
             "sample_count": int(stop - start), "sample_logQ_min": float(sample_q[start]),
             "sample_logQ_max": float(sample_q[stop - 1]),
             "shading_logQ_bounds_midpoint_estimate": [float(edges[start]), float(edges[stop])]}
            for start, stop in zip(changes[:-1], changes[1:])]


def legacy_state(cache, logq, logt, sample_q, sample_t):
    values = baseline.evaluate_spline(cache, logq, logt, sample_q, sample_t)
    density_partial = baseline.evaluate_spline(cache, logq, logt, sample_q, sample_t, (1, 0))
    temperature_partial = baseline.evaluate_spline(cache, logq, logt, sample_q, sample_t, (0, 1))
    field = baseline.FIELDS["Cv"]
    return {"Cv": values[:, field], "Cv_t": (temperature_partial[:, field] - 2 * density_partial[:, field]) / LN10,
            "Pgas": np.exp(values[:, baseline.FIELDS["lnPgas"]])}


def summary(values):
    values = np.asarray(values)
    finite = np.isfinite(values)
    result = {"count": len(values), "nonfinite_count": int((~finite).sum())}
    if finite.any():
        result.update(zip(("abs_median", "abs_p99", "abs_max"),
                          map(float, np.quantile(abs(values[finite]), (.5, .99, 1)))))
    return result


def cache_inventory(directory):
    return {path.name: {"size": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns}
            for path in sorted(directory.glob("*.bin"))}


def measured_fractions(output, sample_q, sample_t):
    mesa = WORKSPACE / "mesa"
    sdk = WORKSPACE / "mesasdk-x86_64-linux-26.6.1"
    probe_dir = CODE / "mesa_probe"
    executable = probe_dir / "build/mesa_probe"
    inlist = probe_dir / "inlist_probe"
    massfile = mesa / "eos/eosFreeEOS_builder/mass_frac_solar.txt"
    if checksum(massfile) != "bd62922804e073a2a410480994b8a28e85b0bc53e469729d6ed5136ea9da10a6":
        raise ValueError("Declared scaled-solar mass pattern changed")
    input_text = "".join(f"{density:.17e} {temperature:.17e}\n"
                         for density, temperature in zip(sample_q + 2 * sample_t - 12, sample_t))
    input_sha = hashlib.sha256(input_text.encode()).hexdigest()
    metadata_path, table_path = output / "probe_metadata.json", output / "probe.tsv"
    fingerprints = {str(path): checksum(path) for path in (executable, inlist, massfile)}
    if metadata_path.exists() and table_path.exists():
        metadata = json.loads(metadata_path.read_text())
        if metadata["input_sha256"] != input_sha or metadata["files_sha256"] != fingerprints:
            raise ValueError("Existing single-call probe does not match this slice; no second probe is permitted")
        if checksum(table_path) != metadata["output_sha256"] or metadata["returncode"] != 0:
            raise ValueError("Saved probe failed its integrity check")
        if cache_inventory(mesa / "data/eosDT_data/cache") != metadata["installed_cache_inventory_before"]:
            raise ValueError("Installed cache inventory changed since the saved single-call probe")
        if any(not path.is_symlink() for path in (Path(metadata["cache_root"]) / "eosDT_cache").glob("*.bin")):
            raise ValueError("Saved probe cache links were replaced")
        metadata = dict(metadata, reused_saved_single_call=True)
    else:
        if metadata_path.exists() or table_path.exists():
            raise ValueError("Partial prior probe artifacts: do not silently invoke a second process")
        source_cache = mesa / "data/eosDT_data/cache"
        before = cache_inventory(source_cache)
        if not before:
            raise ValueError("Existing installed caches are required")
        cache_root, temp_root = output / "runtime_cache", output / "runtime_temp_cache"
        cache_dir = cache_root / "eosDT_cache"
        cache_dir.mkdir(parents=True, exist_ok=True)
        temp_root.mkdir(parents=True, exist_ok=True)
        for name in before:
            destination = cache_dir / name
            target = source_cache / name
            if destination.is_symlink():
                if destination.resolve() != target.resolve():
                    raise ValueError("Unexpected private cache symlink")
            elif destination.exists():
                raise ValueError("Refuse private replacement of an installed cache link")
            else:
                destination.symlink_to(target)
        environment = dict(os.environ)
        environment.update(MESA_DIR=str(mesa), MESASDK_ROOT=str(sdk), OMP_NUM_THREADS="1",
                           MESA_CACHES_DIR=str(cache_root), MESA_TEMP_CACHES_DIR=str(temp_root),
                           LD_LIBRARY_PATH=f"{sdk}/lib:{sdk}/lib64:" + environment.get("LD_LIBRARY_PATH", ""))
        environment.pop("MESA_TEMP_CACHES_DISABLE", None)
        command = [str(executable), str(mesa), str(inlist), str(massfile)]
        (output / "probe_input.txt").write_text(input_text)
        metadata = {"command": command, "input_sha256": input_sha, "files_sha256": fingerprints,
                    "process_invocations": 1, "queries": len(sample_q), "reused_saved_single_call": False,
                    "installed_cache_inventory_before": before, "cache_root": str(cache_root),
                    "temporary_cache_root": str(temp_root), "WDM_AV": False, "X_Y_Z": [.7, .28, .02]}
        metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
        with table_path.open("w") as stdout, (output / "probe.stderr").open("w") as stderr:
            completed = subprocess.run(command, input=input_text, text=True, cwd=output, env=environment,
                                       stdout=stdout, stderr=stderr, timeout=120)
        after = cache_inventory(source_cache)
        replacements = [path.name for path in cache_dir.glob("*.bin") if not path.is_symlink()]
        metadata.update(returncode=completed.returncode, output_sha256=checksum(table_path),
                        installed_cache_inventory_unchanged=before == after, replaced_cache_links=replacements)
        metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
        if completed.returncode or before != after or replacements:
            raise ValueError("Single-call MESA probe/cache safety failure; see probe metadata and stderr")
    data = np.genfromtxt(table_path, names=True)
    if len(data) != len(sample_q) or np.any(data["ierr"] != 0):
        raise ValueError("Probe query count/status mismatch")
    if not np.allclose(data["logRho"], sample_q + 2 * sample_t - 12, rtol=0, atol=1e-12) or not np.allclose(data["logT"], sample_t, rtol=0, atol=1e-12):
        raise ValueError("Probe coordinates differ from the plotted queries")
    fractions = {name: data["frac_" + name] for name in FRACTION_NAMES}
    if not np.allclose(sum(fractions.values()), 1, rtol=0, atol=1e-10):
        raise ValueError("Measured fractions do not sum to one")
    if not np.allclose(fractions["FreeEOS"] + fractions["Skye"], 1, rtol=0, atol=1e-10):
        raise ValueError("Slice is not solely a FreeEOS/Skye transition; do not mislabel the measured blend")
    return data, fractions, metadata


def plot_curves(output, sample_q, states, fractions, spans, source_knots):
    os.environ["MPLCONFIGDIR"] = str(output / ".mplconfig")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(3, 1, figsize=(11.5, 9), sharex=True, constrained_layout=True)
    colors = {"BQ_exact": "tab:blue", "legacy": "tab:orange", "thermo_0.01": "tab:green", "thermo_0.1": "tab:red"}
    labels = {"BQ_exact": "Exact-constraint BQ", "legacy": "Legacy independent Cv", "thermo_0.01": "Archived C4 thermo 0.01", "thermo_0.1": "Archived C4 thermo 0.1"}
    for name, state in states.items():
        axes[0].plot(sample_q, state["Cv"] / 1e8, color=colors[name], label=labels[name], linewidth=1.4)
        axes[1].plot(sample_q, state["Cv_t"] / 1e8, color=colors[name], linewidth=1.4)
    shading = ("#b7dfba", "#ffe0a3", "#d6d6d6")
    for axis in axes:
        for span in spans:
            axis.axvspan(*span["shading_logQ_bounds_midpoint_estimate"], color=shading[span["code"]], alpha=.28, zorder=0)
        for knot in source_knots:
            axis.axvline(knot, color=".6", linewidth=.5, alpha=.45, zorder=0)
        for face in (.77, .80):
            axis.axvline(face, color=".25", linestyle=":", linewidth=1)
        axis.grid(axis="y", alpha=.2)
    axes[0].set_ylabel(r"$C_v$ [$10^8$ erg g$^{-1}$ K$^{-1}$]")
    axes[1].set_ylabel(r"$\partial C_v/\partial\ln T|_\rho$" + "\n" + r"[$10^8$ erg g$^{-1}$ K$^{-1}$]")
    axes[0].legend(loc="best", fontsize=9, ncol=2)
    axes[0].set_title(f"Actual active-FreeEOS slice: logT = {SLICE_LOGT:.15f}; X=.70, Z=.02\nOffline component curves, not runtime-blended Cv; archived C4 coefficients only, no refits", fontsize=11)
    axes[2].plot(sample_q, fractions["FreeEOS"], color="tab:blue", label="Actual MESA fracFreeEOS", linewidth=2)
    axes[2].plot(sample_q, fractions["Skye"], color="tab:purple", label="Actual MESA fracSkye", linewidth=2)
    axes[2].set_ylabel("Actual MESA fractions")
    axes[2].set_ylim(-.05, 1.14)
    axes[2].legend(loc="center left", fontsize=9)
    for span in spans:
        center = sum(span["shading_logQ_bounds_midpoint_estimate"]) / 2
        axes[2].text(center, 1.07, span["label"], ha="center", va="center", fontsize=9)
    axes[2].set_xlabel(r"$\log_{10} Q = \log_{10}\rho - 2\log_{10}T + 12$; thin guides: native Q knots")
    axes[2].set_xlim(.50, 1.10)
    figure.savefig(output / "active_freeeos_slice.png", dpi=180)
    figure.savefig(output / "active_freeeos_slice.pdf")
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    if output != OUTPUT.resolve():
        parser.error("Exclusive outputs must be experiments/023_active_slice")
    output.mkdir(parents=True, exist_ok=True)
    source = WORKSPACE / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
    cache_path = WORKSPACE / "mesa/data/eosDT_data/cache/mesa-FreeEOS_02z70x.bin"
    archive_path = ARCHIVE / "fits_and_queries.npz"
    archive_metadata_path = ARCHIVE / "results.json"
    metadata = json.loads(archive_metadata_path.read_text())
    if checksum(source) != metadata["provenance"]["files"][str(source)]:
        raise ValueError("Source differs from archived 015 active-patch provenance")
    baseline_provenance = json.loads((WORKSPACE / "freeeos_biquintic_progress/experiments/001_native_prototype/results.json").read_text())["provenance"]
    if checksum(cache_path) != baseline_provenance["cache_sha256"]:
        raise ValueError("Existing baseline cache SHA differs; do not rebuild it")
    files = [source, cache_path, archive_path, archive_metadata_path, Path(__file__),
             Path(__file__).with_name("test_slice.py"), CODE / "freeeos_potential.py", CODE / "smooth_patch.py",
             WORKSPACE / "freeeos_report/audit_freeeos.py", Path(__file__).with_name("README.md"),
             CODE / "mesa_probe/probe.f90", CODE / "mesa_probe/inlist_probe"]
    fingerprints = {str(path): checksum(path) for path in files}
    header, logq, logt, rows = read_table(source)
    if not np.allclose(header[1:3], (.7, .02), rtol=0, atol=1e-12):
        raise ValueError("Wrong composition")
    low, valid, _ = source_constraints(rows, logq, logt, split_radiation=False)
    nodes, complete, _ = construct_twists(low, valid, (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10, radius=2, route="mean")
    sample_q, sample_t = slice_queries()
    if np.min(abs(logt - SLICE_LOGT)) <= 1e-10:
        raise ValueError("Required slice temperature must be strictly off all source T knots")
    density_index, temperature_index, inside = query_cells(logq, logt, sample_q, sample_t)
    if not np.all(inside & valid_cells(complete)[density_index, temperature_index]):
        raise ValueError("Slice lacks complete clean exact-BQ stencils")
    _, cache = baseline.read_cache(cache_path, len(logq), len(logt))
    states = {"BQ_exact": thermodynamic_partials(evaluate(nodes, complete, logq, logt, sample_q, sample_t), sample_q, sample_t),
              "legacy": legacy_state(cache, logq, logt, sample_q, sample_t)}
    models = {}
    archive_checks = {}
    with np.load(archive_path, allow_pickle=False) as archive:
        patch_q, patch_t = archive["patch_logQ"], archive["patch_logT"]
        if not (patch_q[0] < sample_q.min() < sample_q.max() < patch_q[-1] and patch_t[0] < SLICE_LOGT < patch_t[-1]):
            raise ValueError("Every plotted query, including start/end, must be strictly inside the archived patch")
        metrics = archive["metrics"].tolist()
        saved_q, saved_t = archive["logQ"][:256], archive["logT"][:256]
        reproduced_bq = thermodynamic_partials(evaluate(nodes, complete, logq, logt, saved_q, saved_t), saved_q, saved_t)
        archive_checks["BQ_exact"] = {field: summary((reproduced_bq[field] - archive["biquintic_values"][metrics.index(field), :256]) /
                                                    np.maximum(1, abs(archive["biquintic_values"][metrics.index(field), :256]))) for field in ("Cv", "Cv_t")}
        archive_checks["BQ_exact"]["Cv_t_field_scaled"] = summary((reproduced_bq["Cv_t"] - archive["biquintic_values"][metrics.index("Cv_t"), :256]) /
            np.maximum.reduce([np.ones(256), abs(archive["biquintic_values"][metrics.index("Cv_t"), :256]), abs(archive["biquintic_values"][metrics.index("Cv"), :256])]))
        for name in CANDIDATES:
            coefficients = archive[name + "_full_coefficients"].copy()
            models[name] = coefficients
            native = evaluate_fit(coefficients, patch_q * LN10, patch_t * LN10, sample_q * LN10, sample_t * LN10)
            states[name] = thermodynamic_partials(native, sample_q, sample_t)
            reproduced = thermodynamic_partials(evaluate_fit(coefficients, patch_q * LN10, patch_t * LN10, saved_q * LN10, saved_t * LN10), saved_q, saved_t)
            archive_checks[name] = {field: summary((reproduced[field] - archive[name + "_values"][metrics.index(field), :256]) /
                                                 np.maximum(1, abs(archive[name + "_values"][metrics.index(field), :256]))) for field in ("Cv", "Cv_t")}
            archive_checks[name]["Cv_t_field_scaled"] = summary((reproduced["Cv_t"] - archive[name + "_values"][metrics.index("Cv_t"), :256]) /
                np.maximum.reduce([np.ones(256), abs(archive[name + "_values"][metrics.index("Cv_t"), :256]), abs(archive[name + "_values"][metrics.index("Cv"), :256])]))
    (output / "archive_reproduction.json").write_text(json.dumps(archive_checks, indent=2, allow_nan=False) + "\n")
    if any(values[field]["abs_max"] > 1e-8 or values[field]["nonfinite_count"] for values in archive_checks.values() for field in ("Cv", "Cv_t_field_scaled")):
        raise ValueError("Current evaluators fail archived coefficient/query reproduction: " + json.dumps(archive_checks))
    positivity = {}
    for name, state in states.items():
        positivity[name] = {field: int((~np.isfinite(state[field]) | (state[field] <= 0)).sum())
                            for field in ("Pgas", "Cv")}
        if any(positivity[name].values()) or not np.all(np.isfinite(state["Cv_t"])):
            raise ValueError(f"Invalid/negative gas pressure or heat capacity in {name}; do not refit or hide it")
    safe = np.min(abs(sample_q[:, None] - logq[None, :]), axis=1) > 2e-4 / LN10 + 1e-12
    fd = {}
    for name, state in states.items():
        fd[name] = {}
        for step in (1e-4, 3e-5, 1e-5):
            plus_q, minus_q = sample_q[safe] - 2 * step / LN10, sample_q[safe] + 2 * step / LN10
            plus_t, minus_t = sample_t[safe] + step / LN10, sample_t[safe] - step / LN10
            if name == "legacy":
                plus = legacy_state(cache, logq, logt, plus_q, plus_t)["Cv"]
                minus = legacy_state(cache, logq, logt, minus_q, minus_t)["Cv"]
            elif name == "BQ_exact":
                plus = thermodynamic_partials(evaluate(nodes, complete, logq, logt, plus_q, plus_t), plus_q, plus_t)["Cv"]
                minus = thermodynamic_partials(evaluate(nodes, complete, logq, logt, minus_q, minus_t), minus_q, minus_t)["Cv"]
            else:
                plus = thermodynamic_partials(evaluate_fit(models[name], patch_q * LN10, patch_t * LN10, plus_q * LN10, plus_t * LN10), plus_q, plus_t)["Cv"]
                minus = thermodynamic_partials(evaluate_fit(models[name], patch_q * LN10, patch_t * LN10, minus_q * LN10, minus_t * LN10), minus_q, minus_t)["Cv"]
            scale = np.maximum.reduce([np.ones(safe.sum()), abs(state["Cv"][safe]), abs(state["Cv_t"][safe])])
            fd[name][str(step)] = summary(((plus - minus) / (2 * step) - state["Cv_t"][safe]) / scale)
    fd_guard_passed = all(metric["nonfinite_count"] == 0 and metric["abs_max"] <= 1e-3
                          for steps in fd.values() for metric in steps.values())
    if not fd_guard_passed:
        raise ValueError("Fixed-rho Cv finite-difference numerical guard failed; no refit or curve suppression is permitted")
    probe_data, fractions, probe_metadata = measured_fractions(output, sample_q, sample_t)
    regions = fraction_regions(fractions["FreeEOS"])
    spans = region_spans(sample_q, regions)
    source_knots = logq[(logq >= .5 - 1e-12) & (logq <= 1.1 + 1e-12)]
    artifact = {"logQ": sample_q, "logT": sample_t, "logRho": sample_q + 2 * sample_t - 12,
                "region_code": regions, "fd_safe_mask": safe, "native_logQ_guides": source_knots}
    for name, state in states.items():
        for field in ("Cv", "Cv_t", "Pgas"):
            artifact[name + "__" + field] = state[field]
    for name, values in fractions.items():
        artifact["MESA_frac_" + name] = values
    for field in probe_data.dtype.names:
        artifact["MESA_probe__" + field] = probe_data[field]
    np.savez_compressed(output / "queries.npz", **artifact)
    plot_curves(output, sample_q, states, fractions, spans, source_knots)
    result = {"status": "Active FreeEOS slice; offline archived-potential curves, NOT fresh accuracy or a runtime C4 implementation",
              "command": [sys.executable, *sys.argv], "command_shell": shlex.join([sys.executable, *sys.argv]),
              "python": sys.version, "numpy": np.__version__, "radiation_constant": float(CRAD),
              "source_sha256": fingerprints[str(source)], "legacy_cache_sha256": fingerprints[str(cache_path)],
              "archived_coefficients_file_sha256": fingerprints[str(archive_path)],
              "files_sha256": fingerprints, "X": .7, "Z": .02, "logT": SLICE_LOGT, "query_count": len(sample_q),
              "logQ_range": [.5, 1.1], "logT_distance_to_nearest_native_knot": float(np.min(abs(logt - SLICE_LOGT))),
              "archived_patch_bounds": {"logQ": patch_q[[0, -1]].tolist(), "logT": patch_t[[0, -1]].tolist()},
              "start_and_end_strict_patch_interior": True, "C4_coefficients": {name: name + "_full_coefficients" for name in CANDIDATES},
              "representation": "Exact constraint radius2/mean total-potential BQ; archived 015 TOTAL-potential C4 full coefficients, no analytic radiation added again; no refits",
              "legacy": "Existing independently interpolated Cv cache; physical Cv_t=(decimal tau partial-2 decimal q partial)/ln10, not MESA runtime-blended Cv_t",
              "physical_derivative": "Cv_t = dCv/dlnT at fixed rho, NOT the slope of this fixed-T density sweep",
              "archive_reproduction_256_saved_queries": archive_checks,
              "archive_reproduction_guard": "Cv relative and Cv_t/max(1,abs(Cv),abs(Cv_t)) errors <=1e-8; raw Cv_t relative errors retained to expose third-derivative cancellation arithmetic, not physical accuracy",
              "positivity_guards": positivity,
              "finite_difference": {"fixed_rho_log_shift": "delta logQ=-2h/ln10, delta logT=h/ln10", "safe_count": int(safe.sum()),
                                    "excluded_knot_centers": int((~safe).sum()), "reason": "Avoid crossing BQ normal-third discontinuities at source Q knots", "metrics_scaled_by_max_Cv_partial_floor1": fd},
              "finite_difference_guard": {"passed": fd_guard_passed, "scaled_max_limit": 1e-3,
                                          "meaning": "Returned-function arithmetic guard, not a physical-accuracy or runtime tolerance"},
              "potential_balanced_identities": {name: {field: summary(values) for field, values in identity_residuals(state, balanced=True).items()}
                                                for name, state in states.items() if name != "legacy"},
              "MESA_probe": probe_metadata, "measured_region_counts": {label: int((regions == index).sum()) for index, label in enumerate(REGION_NAMES)},
              "measured_region_spans": spans,
              "marked_face_fractions": {str(face): {"logQ": float(sample_q[index]), "FreeEOS": float(fractions["FreeEOS"][index]),
                                                   "Skye": float(fractions["Skye"][index]), "label": REGION_NAMES[int(regions[index])]}
                                        for face in (.77, .80) for index in [int(np.argmin(abs(sample_q - face)))]},
              "limitations": ["No refits, source/cache edits or main-report edits; lead integrates plots", "C4 is restricted to the archived patch interior, no global C4 or external blend claim",
                              "Offline curves are shown even where FreeEOS is inactive, explicitly shaded; measured MESA fractions are not applied to these component curves",
                              "No fresh physical accuracy; archived gauge/metal-pattern assumptions remain unverified", "Shading transition bounds are sample midpoint estimates, not exact runtime blend boundaries"]}
    result["output_artifacts_sha256"] = {name: checksum(output / name) for name in
                                        ("queries.npz", "active_freeeos_slice.png", "active_freeeos_slice.pdf", "probe.tsv", "archive_reproduction.json")}
    for path, digest in fingerprints.items():
        if checksum(path) != digest:
            raise ValueError(f"Input/code changed while generating the slice: {path}")
    (output / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    lines = ["# Actual active-FreeEOS slice", "", result["status"], "",
             f"1001 queries: logQ=[.50,1.10], logT={SLICE_LOGT:.15f}, X=.70, Z=.02. No refits.", "",
             "`active_freeeos_slice.png` / `.pdf`: Cv, physical Cv_t at fixed rho, and actual MESA fracFreeEOS/Skye in three panels.",
             "Native Q knots are vertical guides; .77/.80 faces are emphasized. Start/end are strictly inside the archived patch.", "",
             "## Measured selection", ""]
    lines.extend(f"- {span['label']}: {span['sample_count']} points, sampled logQ=[{span['sample_logQ_min']:.6f},{span['sample_logQ_max']:.6f}]." for span in spans)
    lines.extend(["", "Faces: " + json.dumps(result["marked_face_fractions"], sort_keys=True), "",
                  "## Guards and interpretation", "", "Gas pressure and Cv are positive and Cv_t finite for all four offline curves.",
                  f"Fixed-rho Cv finite differences use {int(safe.sum())} common off-knot centers; {int((~safe).sum())} source-knot centers are excluded.",
                  "All three finite-difference steps pass the 1e-3 scaled numerical guard; the limit is not a physical-accuracy tolerance.",
                  "Archived C4/BQ checks use 256 saved 015 queries; Cv-relative and Cv-scaled derivative arithmetic errors are bounded, with raw Cv_t-relative cancellation errors also retained.",
                  "All source/cache/archive/code SHA-256s are recorded and reverified; no coefficient is changed to pass these checks.",
                  "The probe uses one existing executable invocation and exclusively 023-owned cache links/temp directories. Installed cache inventories remain unchanged.",
                  "Offline component curves are not runtime-blended MESA Cv. Active/blended/inactive labels come only from measured fractions.",
                  "No fresh physical accuracy, no refits and no global C4 claim. Prior dense Skye-only slices must not substitute for active FreeEOS evidence.", "",
                  "## Reproduce", "", "```bash", result["command_shell"], "```", "",
                  "A repeated plotting run reuses the integrity-checked saved single-call probe; it does not make a second MESA call.", ""])
    (output / "README.md").write_text("\n".join(lines))
    print(json.dumps({"region_counts": result["measured_region_counts"], "faces": result["marked_face_fractions"],
                      "output": str(output)}, indent=2), flush=True)


if __name__ == "__main__":
    main()
