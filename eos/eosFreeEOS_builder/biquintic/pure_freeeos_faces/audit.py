"""Survey default-MESA routing and original-014 jumps on every native face."""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

BIQUINTIC = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BIQUINTIC))
from freeeos_potential import (CRAD, LN10, checksum, construct_twists, evaluate, read_table,
                               source_constraints, thermodynamic_partials, valid_cells)

WORKSPACE = BIQUINTIC.parents[3]
SOURCE = WORKSPACE / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
EXPERIMENT = WORKSPACE / "freeeos_biquintic_progress/experiments/060_pure_freeeos_faces"
ARCHIVE014 = WORKSPACE / "freeeos_biquintic_progress/experiments/014_active_domain"
PROBE = BIQUINTIC / "mesa_probe/run_probe.py"
TOLERANCE = 1e-12
BATCH = 4096


def native_faces(logq, logt, axis):
    """Return every internal shared face in row-major grid order."""
    if axis == 0:
        density_index, temperature_index = np.indices((len(logq) - 2, len(logt) - 1))
        face_q = logq[density_index + 1]
        face_t = (logt[temperature_index] + logt[temperature_index + 1]) / 2
    else:
        density_index, temperature_index = np.indices((len(logq) - 1, len(logt) - 2))
        face_q = (logq[density_index] + logq[density_index + 1]) / 2
        face_t = logt[temperature_index + 1]
    return tuple(array.ravel() for array in (density_index, temperature_index, face_q, face_t))


def face_support(complete, density_index, temperature_index, axis):
    cells = valid_cells(complete)
    if axis == 0:
        return cells[density_index, temperature_index] & cells[density_index + 1, temperature_index]
    return cells[density_index, temperature_index] & cells[density_index, temperature_index + 1]


def nearby_nodes(mask, density_index, temperature_index, axis):
    result = np.zeros(len(density_index), dtype=bool)
    density_width, temperature_width = (3, 2) if axis == 0 else (2, 3)
    for density_offset in range(density_width):
        for temperature_offset in range(temperature_width):
            result |= mask[density_index + density_offset, temperature_index + temperature_offset]
    return result


def support_causes(rows, logt, source_valid, complete, density_index, temperature_index, axis):
    supported = face_support(complete, density_index, temperature_index, axis)
    cold = nearby_nodes(np.broadcast_to(logt[None, :] < 3.1 - 1e-12, source_valid.shape),
                        density_index, temperature_index, axis)
    floor = nearby_nodes(np.broadcast_to((logt[None, :] >= 3.1 - 1e-12) &
                                                   (logt[None, :] < 3.2 - 1e-12), source_valid.shape),
                         density_index, temperature_index, axis)
    source_failure = nearby_nodes((rows[:, :, 17] != 0) & (logt[None, :] >= 3.2 - 1e-12),
                                  density_index, temperature_index, axis)
    source_failure |= nearby_nodes((~source_valid) & (logt[None, :] >= 3.2 - 1e-12),
                                   density_index, temperature_index, axis)
    cause = np.zeros(len(density_index), dtype=np.uint8)
    cause[~supported] = 4
    cause[~supported & source_failure] = 3
    cause[~supported & floor] = 2
    cause[~supported & cold] = 1
    return supported, cause


def classify_probe(probe):
    status = probe[:, 2]
    fractions = probe[:, 3:11]
    finite = np.all(np.isfinite(probe), axis=1)
    successful = (status == 0) & finite
    sums = np.sum(fractions, axis=1)
    fraction = fractions[:, 2]
    good_fractions = (np.abs(sums - 1) <= TOLERANCE) & np.all(fractions >= -TOLERANCE, axis=1)
    good_fractions &= np.all(fractions <= 1 + TOLERANCE, axis=1)
    pure = successful & good_fractions & (np.abs(fraction - 1) <= TOLERANCE)
    pure &= np.sum(np.abs(np.delete(fractions, 2, axis=1)), axis=1) <= TOLERANCE
    return successful, pure, good_fractions


def locations(mask, face_q, face_t, limit=8):
    return [{"logQ": float(face_q[index]), "logT": float(face_t[index]),
             "logRho": float(face_q[index] + 2 * face_t[index] - 12)}
            for index in np.flatnonzero(mask)[:limit]]


def coordinate_range(mask, face_q, face_t):
    if not np.any(mask):
        return {}
    density = face_q[mask] + 2 * face_t[mask] - 12
    return {"logQ": [float(np.min(face_q[mask])), float(np.max(face_q[mask]))],
            "logT": [float(np.min(face_t[mask])), float(np.max(face_t[mask]))],
            "logRho": [float(np.min(density)), float(np.max(density))]}


def summarize_jump(jump, mask, face_q, face_t):
    selected = np.abs(jump[mask])
    if not len(selected):
        return {"count": 0}
    worst = np.flatnonzero(mask)[np.argmax(selected)]
    return {"count": int(len(selected)), "median": float(np.median(selected)),
            "p90": float(np.quantile(selected, .9)), "p99": float(np.quantile(selected, .99)),
            "max": float(np.max(selected)), "exact_zero": int(np.count_nonzero(jump[mask] == 0)),
            "below_display_min_1e-7": int(np.count_nonzero(selected < 1e-7)),
            "above_display_max_10": int(np.count_nonzero(selected > 10)),
            "worst": {**locations(np.arange(len(jump)) == worst, face_q, face_t, 1)[0],
                      "signed_scaled_jump": float(jump[worst])}}


def prepare(output):
    _, logq, logt, _ = read_table(SOURCE)
    arrays = {}
    with (output / "face_queries.txt").open("w") as stream:
        for axis, direction in ((0, "u"), (1, "tau")):
            density_index, temperature_index, face_q, face_t = native_faces(logq, logt, axis)
            arrays.update({f"{direction}_{name}": array for name, array in
                           (("density_index", density_index), ("temperature_index", temperature_index),
                            ("logQ", face_q), ("logT", face_t))})
            for density, temperature in zip(face_q + 2 * face_t - 12, face_t):
                stream.write(f"{density:.17g} {temperature:.17g}\n")
    np.savez_compressed(output / "face_manifest.npz", **arrays)
    print({"u": len(arrays["u_logQ"]), "tau": len(arrays["tau_logQ"]),
           "total": len(arrays["u_logQ"]) + len(arrays["tau_logQ"])})


def probe_faces(output):
    command = [sys.executable, "-B", str(PROBE), "--mesa-dir", str(WORKSPACE / "mesa"),
               "--output", str(output)]
    (output / "probe_command.json").write_text(json.dumps(command, indent=2) + "\n")
    with (output / "face_queries.txt").open() as input_stream, \
         (output / "face_probe.tsv").open("w") as output_stream, \
         (output / "face_probe.stderr").open("w") as error_stream:
        completed = subprocess.run(command, stdin=input_stream, stdout=output_stream,
                                   stderr=error_stream, check=False)
    (output / "probe_returncode.txt").write_text(f"{completed.returncode}\n")
    if completed.returncode:
        raise RuntimeError(f"Probe returned {completed.returncode}")


def overlap_control(nodes, complete, logq, logt):
    frozen = np.load(ARCHIVE014 / "active_faces.npz")
    result = {}
    for axis, direction in ((0, "u"), (1, "tau")):
        face_q, face_t = frozen[f"{direction}_q"], frozen[f"{direction}_t"]
        density_index = np.searchsorted(logq, face_q, side="left") - 1
        temperature_index = np.searchsorted(logt, face_t, side="left") - 1
        density_index = np.clip(density_index, 0, len(logq) - 2)
        temperature_index = np.clip(temperature_index, 0, len(logt) - 2)
        supported = face_support(complete, density_index, temperature_index, axis)
        chosen = np.flatnonzero(supported)[:128]
        if not len(chosen):
            raise ValueError("No frozen014 overlap is supported")
        left_index = density_index[chosen], temperature_index[chosen]
        right_index = (density_index[chosen] + (axis == 0), temperature_index[chosen] + (axis == 1))
        left = thermodynamic_partials(evaluate(nodes, complete, logq, logt, face_q[chosen], face_t[chosen],
                                               cell_indices=left_index), face_q[chosen], face_t[chosen])
        right = thermodynamic_partials(evaluate(nodes, complete, logq, logt, face_q[chosen], face_t[chosen],
                                                cell_indices=right_index), face_q[chosen], face_t[chosen])
        scale = np.maximum.reduce([abs(left["Cv"]), abs(right["Cv"]), abs(left["Cv_t"]),
                                   abs(right["Cv_t"]), np.ones(len(chosen))])
        current = (left["Cv_t"] - right["Cv_t"]) / scale
        old = frozen[f"{direction}_Cv_t_jump"][chosen]
        frozen_pure = frozen[f"{direction}_probe"][chosen, 5] == 1
        if not np.any(frozen_pure):
            raise ValueError("Frozen014 overlap lacks exactly pure routing")
        result[direction] = {"checked": int(len(chosen)), "max_absolute_difference": float(np.max(abs(current - old))),
                             "median_absolute_difference": float(np.median(abs(current - old))),
                             "frozen_exact_pure_checked": int(frozen_pure.sum()),
                             "max_absolute_difference_on_frozen_exact_pure": float(np.max(abs(current[frozen_pure] - old[frozen_pure]))),
                             "max_frozen_jump": float(np.max(abs(old)))}
    return result


def analyze(output):
    _, logq, logt, rows = read_table(SOURCE)
    low, source_valid, _ = source_constraints(rows, logq, logt, split_radiation=False, minimum_logT=3.2)
    nodes, complete, _ = construct_twists(low, source_valid, (logq[1] - logq[0]) * LN10,
                                          (logt[1] - logt[0]) * LN10, radius=2, route="mean")
    manifest = np.load(output / "face_manifest.npz")
    probe = np.loadtxt(output / "face_probe.tsv", skiprows=1, ndmin=2)
    nominal = sum(len(manifest[f"{direction}_logQ"]) for direction in ("u", "tau"))
    if probe.shape != (nominal, 14):
        raise ValueError(f"Probe shape {probe.shape}, expected ({nominal}, 14)")
    all_q = np.concatenate([manifest[f"{direction}_logQ"] for direction in ("u", "tau")])
    all_t = np.concatenate([manifest[f"{direction}_logT"] for direction in ("u", "tau")])
    if not np.allclose(probe[:, :2], np.column_stack((all_q + 2 * all_t - 12, all_t)), atol=1e-12, rtol=0):
        raise ValueError("Probe output does not follow complete manifest")
    summary = {"nominal_faces": nominal, "tolerance": TOLERANCE,
               "method": "all internal faces, tangential midpoint, actual default-MESA fractions before support mask",
               "interpolant": "original014 total-potential radius2 mean, source floor 3.2, split_radiation=False",
               "support_codes": {"0": "supported", "1": "cold logT<3.1", "2": "source floor 3.1<=logT<3.2",
                                 "3": "other source invalid", "4": "centered stencil/edge margin"},
               "overlap014": overlap_control(nodes, complete, logq, logt), "faces": {},
               "provenance": {"source_sha256": checksum(SOURCE), "audit_script_sha256": checksum(__file__),
                              "preprocessing_sha256": checksum(BIQUINTIC / "freeeos_potential.py"),
                              "probe_script_sha256": checksum(PROBE), "probe_executable_sha256": checksum(BIQUINTIC / "mesa_probe/build/mesa_probe"),
                              "massfile_sha256": checksum(WORKSPACE / "mesa/eos/eosFreeEOS_builder/mass_frac_solar.txt"),
                              "queries_sha256": checksum(output / "face_queries.txt"), "manifest_sha256": checksum(output / "face_manifest.npz"),
                              "probe_output_sha256": checksum(output / "face_probe.tsv"),
                              "old014_archive_sha256": checksum(ARCHIVE014 / "active_faces.npz")}}
    archive = {}
    offset = 0
    for axis, direction in ((0, "u"), (1, "tau")):
        face_q, face_t = manifest[f"{direction}_logQ"], manifest[f"{direction}_logT"]
        density_index, temperature_index = manifest[f"{direction}_density_index"], manifest[f"{direction}_temperature_index"]
        count = len(face_q)
        measured = probe[offset:offset + count]
        offset += count
        successful, pure, good_fractions = classify_probe(measured)
        supported, cause = support_causes(rows, logt, source_valid, complete, density_index, temperature_index, axis)
        selected = pure & supported
        if axis == 0:
            face_log_pgas = (rows[density_index + 1, temperature_index, 1] +
                             rows[density_index + 1, temperature_index + 1, 1]) / 2
        else:
            face_log_pgas = (rows[density_index, temperature_index + 1, 1] +
                             rows[density_index + 1, temperature_index + 1, 1]) / 2
        log_radiation_ratio = np.log10(CRAD / 3) + 4 * face_t - face_log_pgas
        jumps = np.full(count, np.nan)
        for start in range(0, count, BATCH):
            chosen = np.flatnonzero(selected[start:start + BATCH]) + start
            if not len(chosen):
                continue
            left_index = density_index[chosen], temperature_index[chosen]
            right_index = (density_index[chosen] + (axis == 0), temperature_index[chosen] + (axis == 1))
            left = thermodynamic_partials(evaluate(nodes, complete, logq, logt, face_q[chosen], face_t[chosen],
                                                   cell_indices=left_index), face_q[chosen], face_t[chosen])
            right = thermodynamic_partials(evaluate(nodes, complete, logq, logt, face_q[chosen], face_t[chosen],
                                                    cell_indices=right_index), face_q[chosen], face_t[chosen])
            scale = np.maximum.reduce([abs(left["Cv"]), abs(right["Cv"]), abs(left["Cv_t"]),
                                       abs(right["Cv_t"]), np.ones(len(chosen))])
            jumps[chosen] = (left["Cv_t"] - right["Cv_t"]) / scale
        finite_jump = selected & np.isfinite(jumps)
        status, status_count = np.unique(measured[:, 2], return_counts=True)
        summary["faces"][direction] = {"nominal": count, "probe_status_counts": {str(int(key)): int(value) for key, value in zip(status, status_count)},
                                       "successful_finite": int(successful.sum()), "invalid_fraction_sum_or_range": int((successful & ~good_fractions).sum()),
                                       "pure": int(pure.sum()), "pure_exact_one": int((pure & (measured[:, 5] == 1)).sum()),
                                       "pure_tolerance_only": int((pure & (measured[:, 5] != 1)).sum()),
                                       "pure_supported": int(selected.sum()), "nonfinite_jump": int((selected & ~np.isfinite(jumps)).sum()),
                                       "source_nodal_logPrad_over_Pgas_midpoint_proxy": {
                                           "supported_pure_ratio_gt_1": int((selected & (log_radiation_ratio > 0)).sum()),
                                           "supported_pure_ratio_gt_100": int((selected & (log_radiation_ratio > 2)).sum()),
                                           "supported_pure_max_log10_ratio": float(np.max(log_radiation_ratio[selected]))},
                                       "jump": summarize_jump(jumps, finite_jump, face_q, face_t),
                                       "pure_support_gaps": {str(code): {"count": int((pure & (cause == code)).sum()),
                                                                          "coordinate_range": coordinate_range(pure & (cause == code), face_q, face_t),
                                                                          "examples": locations(pure & (cause == code), face_q, face_t)}
                                                             for code in range(1, 5)},
                                       "all_support_gaps": {str(code): int((cause == code).sum()) for code in range(1, 5)},
                                       "pure_q_range": [float(np.min(face_q[pure])), float(np.max(face_q[pure]))] if pure.any() else [],
                                       "pure_t_range": [float(np.min(face_t[pure])), float(np.max(face_t[pure]))] if pure.any() else []}
        archive.update({f"{direction}_probe": measured, f"{direction}_pure": pure, f"{direction}_support_code": cause,
                        f"{direction}_Cv_t_jump": jumps, f"{direction}_source_logPrad_over_Pgas_proxy": log_radiation_ratio})
    np.savez_compressed(output / "face_results.npz", **archive)
    (output / "results.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({direction: {key: value for key, value in data.items() if key in
                       ("nominal", "pure", "pure_supported", "pure_support_gaps", "jump")}
                      for direction, data in summary["faces"].items()}, indent=2))


def plot(output):
    os.environ.setdefault("MPLCONFIGDIR", str(output / ".mplconfig"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    manifest = np.load(output / "face_manifest.npz")
    results = np.load(output / "face_results.npz")
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharex=True, sharey=True, constrained_layout=True)
    norm = TwoSlopeNorm(vmin=-5, vcenter=-2, vmax=1)
    for axis, direction in zip(axes, ("u", "tau")):
        face_q, face_t = manifest[f"{direction}_logQ"], manifest[f"{direction}_logT"]
        density = face_q + 2 * face_t - 12
        pure = results[f"{direction}_pure"]
        support = results[f"{direction}_support_code"] == 0
        jump = results[f"{direction}_Cv_t_jump"]
        gap = pure & ~support
        chosen = pure & support & np.isfinite(jump)
        axis.scatter(density[gap], face_t[gap], s=2.2, c=".55", rasterized=True, label="Pure routing; no support")
        image = axis.scatter(density[chosen], face_t[chosen], s=2.2,
                             c=np.log10(np.maximum(abs(jump[chosen]), 1e-17)),
                             cmap="coolwarm", norm=norm, rasterized=True)
        axis.set_title(f"{direction} faces: {chosen.sum():,} supported; {gap.sum():,} gaps")
        axis.set_xlabel(r"$\log_{10}\rho$ [g cm$^{-3}$]")
        axis.grid(alpha=.15)
    axes[0].set_ylabel(r"$\log_{10}T$ [K]")
    axes[1].legend(loc="lower right", fontsize=8, markerscale=3)
    figure.colorbar(image, ax=axes, label=r"$\log_{10}$ scaled $|\Delta(\partial C_V/\partial\ln T)|$", extend="both")
    figure.suptitle("FreeEOS-only one-sided jumps where default MESA is fully FreeEOS\n"
                    "X=.70, Y=.28, Z=.02; all native faces sampled at tangential midpoints", fontsize=11)
    presentation = {"cmap": "coolwarm", "vmin": -5, "vcenter": -2, "vmax": 1,
                    "results_sha256": checksum(output / "face_results.npz"),
                    "below_display": {direction: int(np.count_nonzero(
                        results[f"{direction}_pure"] & (results[f"{direction}_support_code"] == 0) &
                        (abs(results[f"{direction}_Cv_t_jump"]) < 1e-5))) for direction in ("u", "tau")}}
    (output / "figure_presentation.json").write_text(json.dumps(presentation, indent=2) + "\n")
    for extension in ("pdf", "png"):
        figure.savefig(output / f"pure_freeeos_Cv_faces.{extension}", dpi=200)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "probe", "analyze", "plot"))
    parser.add_argument("--output", type=Path, default=EXPERIMENT)
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    {"prepare": prepare, "probe": probe_faces, "analyze": analyze, "plot": plot}[arguments.stage](output)


if __name__ == "__main__":
    main()
