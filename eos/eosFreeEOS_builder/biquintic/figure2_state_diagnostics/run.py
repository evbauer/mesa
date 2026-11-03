"""Figure 2 state and native-face diagnostics from the archived 036 fits."""

import argparse
import json
import os
from pathlib import Path
import sys

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"
sys.dont_write_bytecode = True

import numpy as np

CODE = Path(__file__).resolve().parents[1]
WORKSPACE = CODE.parents[3]
ARCHIVE = WORKSPACE / "freeeos_biquintic_progress/experiments/036_dense_derivatives"
REPLAY = WORKSPACE / "freeeos_biquintic_progress/experiments/068_figure2_consistency"
OUTPUT = WORKSPACE / "freeeos_biquintic_progress/experiments/072_figure2_state_diagnostics"
sys.path.insert(0, str(CODE))
sys.path.insert(0, str(CODE / "dense_slice_consistency"))
sys.path.insert(0, str(WORKSPACE / "freeeos_report"))
import audit_freeeos as legacy
import run_figure2 as figure2
from freeeos_potential import (CRAD, LN10, checksum, construct_twists,
                               evaluate, quintic_basis, read_table, source_constraints,
                               thermodynamic_partials, thermodynamics)
from smooth_patch import evaluate_fit

NAMES = ("biquintic", "native_1.0", "thermo_0.1", "legacy")
FIELDS = ("P", "E", "S")
FACE_FIELDS = ("P", "E", "S", "P_x", "E_x", "S_x", "P_t", "E_t", "S_t",
               "Cv", "Cv_x", "Cv_t")


def compatibility(values, slopes, curvatures, spacing):
    """Right minus left third derivative at the middle of three knots."""
    return (60 * (values[..., 2] - 2 * values[..., 1] + values[..., 0]) / spacing**3
            - 24 * (slopes[..., 2] - slopes[..., 0]) / spacing**2
            + 3 * (curvatures[..., 2] - 6 * curvatures[..., 1] + curvatures[..., 0]) / spacing)


def state_from_native(native, sample_q, sample_t, face=False):
    state = thermodynamic_partials(native, sample_q, sample_t) if face else thermodynamics(native, sample_q, sample_t)
    result = {field: state[field] for field in FIELDS}
    if face:
        density = state["rho"]
        temperature = state["T"]
        result.update({field + "_x": density * state[field + "_rho"] for field in FIELDS})
        result.update({field + "_t": temperature * state[field + "_T"] for field in FIELDS})
        result.update({field: state[field] for field in ("Cv", "Cv_x", "Cv_t")})
    return result


def legacy_state(cache, logq, logt, sample_q, sample_t):
    values = legacy.evaluate_spline(cache, logq, logt, sample_q, sample_t)
    temperature = 10.0**sample_t
    return {"P": np.exp(values[:, legacy.FIELDS["lnPgas"]]) + CRAD * temperature**4 / 3,
            "E": np.exp(values[:, legacy.FIELDS["lnE"]]),
            "S": np.exp(values[:, legacy.FIELDS["lnS"]])}


def scaled_jump(left, right, field):
    names = ("Cv", field)
    scale = np.maximum.reduce([abs(left[name]) for name in names] + [abs(right[name]) for name in names])
    if not np.all(np.isfinite(scale)) or np.any(scale == 0):
        if np.any(~np.isfinite(scale)):
            raise ValueError("Nonfinite jump scale")
    return np.divide(right[field] - left[field], scale,
                     out=np.zeros_like(scale), where=scale != 0), scale


def relative_summary(values):
    values = np.asarray(values)
    if not np.all(np.isfinite(values)):
        raise ValueError("Nonfinite diagnostic")
    return {"min": float(np.min(values)), "max": float(np.max(values)),
            "median_abs": float(np.median(abs(values))), "max_abs": float(np.max(abs(values)))}


def face_density_difference(nodes, complete, logq, logt, face_q, face_t, face_index,
                            temperature_index, side, step):
    offsets = np.array([-2, -1, 0] if side == "left" else [0, 1, 2])
    positions = face_q + step * offsets
    native = evaluate(nodes, complete, logq, logt, positions,
                      np.full(3, face_t), cell_indices=(
                          np.full(3, face_index - (side == "left")),
                          np.full(3, temperature_index)))
    cv = thermodynamics(native, positions, np.full(3, face_t))["Cv"]
    numerator = (3 * cv[2] - 4 * cv[1] + cv[0] if side == "left"
                 else -3 * cv[0] + 4 * cv[1] - cv[2])
    return numerator / (2 * step * LN10)


def plot(arrays, metadata):
    os.environ["MPLCONFIGDIR"] = str(OUTPUT / ".mplconfig")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"biquintic": "tab:blue", "native_1.0": "tab:green",
              "thermo_0.1": "tab:red", "legacy": "tab:orange"}
    labels = {"biquintic": "Biquintic", "native_1.0": "Native $C^4$ fit",
              "thermo_0.1": "Thermodynamic $C^4$ fit", "legacy": "Legacy bicubic"}
    source_labels = (f"Source row $\\log_{{10}}T={value:.2f}$" for value in metadata["bracketing_logT"])
    source_labels = tuple(source_labels)
    density = arrays["sample_logQ"] + 2 * arrays["sample_logT"] - 12
    face_density = arrays["face_logQ"] + 2 * arrays["sample_logT"][0] - 12
    units = {"P": r"dyn cm$^{-2}$", "E": r"erg g$^{-1}$", "S": r"erg g$^{-1}$ K$^{-1}$"}

    for kind in ("raw", "residual"):
        figure, axes = plt.subplots(3, 1, figsize=(9.1, 8.2), sharex=True, constrained_layout=True)
        for axis, field in zip(axes, FIELDS):
            reference = arrays[f"thermo_0.1_{field}"]
            if not np.all(reference > 0):
                raise ValueError(f"Reference {field} is not positive")
            for name in NAMES:
                values = arrays[f"{name}_{field}"]
                display = np.log10(values) if kind == "raw" else 1e6 * (values - reference) / reference
                axis.plot(density, display, color=colors[name],
                          linestyle="--" if name == "legacy" else "-", linewidth=1.25,
                          label=labels[name], zorder=4 if name == "legacy" else 2)
            for row_index, (marker, row_label) in enumerate(zip(("o", "s"), source_labels)):
                values = arrays[f"source_{row_index}_{field}"]
                row_reference = arrays[f"source_{row_index}_reference_{field}"]
                display = np.log10(values) if kind == "raw" else 1e6 * (values - row_reference) / row_reference
                axis.scatter(arrays[f"source_{row_index}_logRho"], display,
                             s=22, facecolors="none", edgecolors=("#333333", "#8a168a")[row_index],
                             marker=marker, label=row_label, zorder=5)
            for coordinate in face_density:
                axis.axvline(coordinate, color="0.7", alpha=.25, linewidth=.6)
            axis.set_ylabel((rf"$\log_{{10}} {field}$ [{units[field]}]" if kind == "raw"
                             else rf"$({field}-{field}_{{\rm ref}})/{field}_{{\rm ref}}$ [ppm]"))
        axes[0].legend(fontsize=8, ncol=3)
        axes[-1].set_xlabel(r"$\log_{10}\rho$ [g cm$^{-3}$]")
        figure.suptitle("Figure 2 FreeEOS-only state functions; default MESA selection is Skye-only\n"
                        + ("Actual source rows are shown at their own temperatures and densities"
                           if kind == "raw" else "Common reference: original red $C^4$ potential, evaluated at each point's own state"),
                        fontsize=10)
        for extension in ("pdf", "png"):
            figure.savefig(OUTPUT / f"state_{kind}.{extension}", dpi=180)
        plt.close(figure)

    figure, axes = plt.subplots(2, 1, figsize=(9.0, 6.8), sharex=True, constrained_layout=True)
    face_q = arrays["face_logQ"]
    row_labels = (f"Knot line $\\log_{{10}}T={value:.2f}$" for value in metadata["bracketing_logT"])
    for label, key, color, marker in (("Figure 2 between rows", "off", "tab:blue", "o"),
                                      (next(row_labels), "row_0", "#333333", "s"),
                                      (next(row_labels), "row_1", "#8a168a", "^")):
        axes[0].plot(face_q, arrays[f"{key}_Cv_x_jump"] / 1e8, marker=marker,
                     color=color, label=label, linewidth=1.2, markersize=4)
        axes[1].plot(face_q, arrays[f"{key}_Cv_x_scaled_jump"], marker=marker,
                     color=color, linewidth=1.2, markersize=4)
    for axis in axes:
        axis.axhline(0, color="0.55", linewidth=.7)
        for key in ("large_face", "quiet_face"):
            index = metadata[key]["index"]
            axis.axvline(face_q[index], color="0.5", linestyle=":" if key == "quiet_face" else "--", linewidth=.8)
    axes[0].legend(fontsize=8)
    axes[0].set_ylabel("Right $-$ left $C_{V,x}$ / $10^8$\n[erg g$^{-1}$ K$^{-1}$]")
    axes[1].set_ylabel("Jump / max($|C_V|$, $|C_{V,x}|$)\n[each face, both sides]")
    axes[1].set_xlabel(r"Native face $\log_{10}Q$ (row densities differ at fixed $Q$)")
    figure.suptitle("Native $Q$ face compatibility: exact source-knot lines and off-row Figure 2\n"
                    "Off-row jump also depends on estimated mixed jets; these are not source points", fontsize=10)
    for extension in ("pdf", "png"):
        figure.savefig(OUTPUT / f"face_compatibility.{extension}", dpi=180)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot-only", action="store_true")
    arguments = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    arrays_path = OUTPUT / "diagnostics.npz"
    results_path = OUTPUT / "results.json"
    if arguments.plot_only:
        metadata = json.loads(results_path.read_text())
        if checksum(arrays_path) != metadata["output_sha256"]["diagnostics.npz"]:
            raise ValueError("Saved 072 arrays fail checksum")
        with np.load(arrays_path) as saved:
            plot({key: saved[key] for key in saved.files}, metadata)
        metadata["output_sha256"].update({path.name: checksum(path) for path in OUTPUT.glob("*.pdf")})
        metadata["output_sha256"].update({path.name: checksum(path) for path in OUTPUT.glob("*.png")})
        metadata["script_sha256"] = checksum(__file__)
        results_path.write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n")
        return

    archived = json.loads((ARCHIVE / "results.json").read_text())
    replayed = json.loads((REPLAY / "results.json").read_text())
    inputs = {**archived["input_sha256"],
              "freeeos_biquintic_progress/experiments/036_dense_derivatives/curves_and_traces.npz":
                  archived["artifact_sha256"]["curves_and_traces.npz"],
              "freeeos_biquintic_progress/experiments/068_figure2_consistency/curves_and_traces.npz":
                  replayed["output_sha256"]["curves_and_traces.npz"]}
    for name, expected in inputs.items():
        if checksum(WORKSPACE / name) != expected:
            raise ValueError(f"Read-only input hash mismatch: {name}")
    with np.load(ARCHIVE / "curves_and_traces.npz") as saved:
        archive = {key: saved[key] for key in saved.files}
    sample_q, sample_t, face_q = (archive[key] for key in ("sample_logQ", "sample_logT", "face_logQ"))
    if len(sample_q) != 2001 or not np.array_equal(sample_q, np.linspace(2.25, 2.75, 2001)) or not np.all(sample_t == 6.5117106292763225):
        raise ValueError("Original 036 trajectory changed")
    table = WORKSPACE / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
    cache_path = WORKSPACE / "mesa/data/eosDT_data/cache/mesa-FreeEOS_02z70x.bin"
    _, logq, logt, rows = read_table(table)
    low, valid, source = source_constraints(rows, logq, logt)
    step_q, step_t = np.diff(logq)[0] * LN10, np.diff(logt)[0] * LN10
    nodes, complete, _ = construct_twists(low, valid, step_q, step_t)
    _, cache = legacy.read_cache(cache_path, len(logq), len(logt))
    fit_q, fit_t = archive["patch_logQ"] * LN10, archive["patch_logT"] * LN10
    arrays = {"sample_logQ": sample_q, "sample_logT": sample_t, "face_logQ": face_q}
    native = evaluate(nodes, complete, logq, logt, sample_q, sample_t)
    states = {"biquintic": state_from_native(native, sample_q, sample_t),
              "legacy": legacy_state(cache, logq, logt, sample_q, sample_t)}
    for name in ("native_1.0", "thermo_0.1"):
        native_fit = evaluate_fit(archive[f"{name}_coefficients"], fit_q, fit_t,
                                  sample_q * LN10, sample_t * LN10)
        states[name] = state_from_native(native_fit, sample_q, sample_t)
    for name, state in states.items():
        for field in FIELDS:
            arrays[f"{name}_{field}"] = state[field]
    reproduction = {}
    for name in NAMES:
        candidate = (figure2.legacy_state(cache, logq, logt, sample_q, sample_t) if name == "legacy"
                     else figure2.potential_state(
                         native if name == "biquintic" else evaluate_fit(
                             archive[f"{name}_coefficients"], fit_q, fit_t,
                             sample_q * LN10, sample_t * LN10), sample_q, sample_t))
        reproduction[name] = {}
        for field in ("Cv", "Cv_x", "Cv_t"):
            error = (candidate[field] - archive[f"{name}_{field}"]) / np.maximum(1, abs(archive[f"{name}_{field}"]))
            reproduction[name][field] = float(np.max(abs(error)))
            if reproduction[name][field] > 1e-11:
                raise AssertionError(f"036 replay failed: {name} {field}")

    face_indices = np.searchsorted(logq, face_q)
    if not np.array_equal(logq[face_indices], face_q):
        raise ValueError("Figure 2 faces no longer match source Q knots")
    lower_row = np.searchsorted(logt, sample_t[0]) - 1
    row_indices = (lower_row, lower_row + 1)
    face_t = np.full_like(face_q, sample_t[0])
    face_sides = {}
    for side, offset in (("left", -1), ("right", 0)):
        face_native = evaluate(nodes, complete, logq, logt, face_q, face_t,
                               cell_indices=(face_indices + offset, np.full(len(face_q), lower_row)))
        face_sides[side] = state_from_native(face_native, face_q, face_t, face=True)
        arrays[f"off_{side}_f_uuu"] = face_native[3, 0]
        for field in FACE_FIELDS:
            arrays[f"off_{side}_{field}"] = face_sides[side][field]
        for field in ("Cv", "Cv_x", "Cv_t"):
            error = (face_sides[side][field] - archive[f"{side}_{field}"]) / np.maximum(1, abs(archive[f"{side}_{field}"]))
            reproduction[f"face_{side}_{field}"] = float(np.max(abs(error)))
            if reproduction[f"face_{side}_{field}"] > 1e-11:
                raise AssertionError(f"036 face replay failed: {side} {field}")

    continuity = {}
    for field in FACE_FIELDS[:9]:
        left, right = face_sides["left"][field], face_sides["right"][field]
        continuity[field] = float(np.max(abs(right - left) / np.maximum(abs(left), abs(right))))
    for field in ("Cv_x", "Cv_t"):
        scaled, scale = scaled_jump(face_sides["left"], face_sides["right"], field)
        arrays[f"off_{field}_jump"] = face_sides["right"][field] - face_sides["left"][field]
        arrays[f"off_{field}_scaled_jump"] = scaled
        arrays[f"off_{field}_scale"] = scale

    q_triplets = np.stack([logq[face_indices - 1], logq[face_indices], logq[face_indices + 1]], axis=1)
    effective = evaluate(nodes, complete, logq, logt, q_triplets.ravel(),
                         np.full(q_triplets.size, sample_t[0]))
    effective_j = compatibility(*(effective[order, 0].reshape(len(face_q), 3)
                                  for order in range(3)), step_q)
    direct_j = -arrays["off_Cv_x_jump"] * 10.0**sample_t[0] / 4
    arrays["off_effective_J"] = effective_j
    arrays["off_direct_J"] = direct_j
    effective_error = float(np.max(abs(-4 * (effective_j - direct_j) /
                                       10.0**sample_t[0]) / arrays["off_Cv_x_scale"]))
    direct_native_error = float(np.max(abs(arrays["off_right_f_uuu"] - arrays["off_left_f_uuu"]
                                           - direct_j) * 4 / 10.0**sample_t[0]
                                       / arrays["off_Cv_x_scale"]))
    physical_shear_error = float(np.max(abs(arrays["off_Cv_t_jump"] + 2 * arrays["off_Cv_x_jump"])
                                        / arrays["off_Cv_x_scale"]))

    tangent_position = np.array([(sample_t[0] - logt[lower_row]) / np.diff(logt)[0]])
    tangent_basis = quintic_basis(tangent_position, step_t, 0)[0]
    supplied_contribution = np.zeros_like(face_q)
    estimated_contribution = np.zeros_like(face_q)
    for derivative_order in range(3):
        contribution = np.zeros_like(face_q)
        for endpoint, row_index in enumerate(row_indices):
            line_jets = [
                nodes[density_order, derivative_order,
                      face_indices[:, None] + np.array([-1, 0, 1]), row_index]
                for density_order in range(3)]
            tangent_line_j = compatibility(*line_jets, step_q)
            supplied_line_jets = [np.zeros_like(values) if (density_order, derivative_order)
                                  in ((2, 1), (1, 2), (2, 2)) else values
                                  for density_order, values in enumerate(line_jets)]
            supplied_line_j = compatibility(*supplied_line_jets, step_q)
            arrays[f"row_{endpoint}_J_tau_order_{derivative_order}"] = tangent_line_j
            weight = tangent_basis[endpoint, derivative_order]
            contribution += weight * tangent_line_j
            supplied_contribution += weight * supplied_line_j
            estimated_contribution += weight * (tangent_line_j - supplied_line_j)
        arrays[f"off_J_component_tau_order_{derivative_order}"] = contribution
    arrays["off_J_component_supplied"] = supplied_contribution
    arrays["off_J_component_estimated"] = estimated_contribution
    tangent_reconstruction_error = float(np.max(abs(
        sum(arrays[f"off_J_component_tau_order_{order}"] for order in range(3))
        - effective_j) * 4 / 10.0**sample_t[0] / arrays["off_Cv_x_scale"]))
    supplied_estimated_reconstruction_error = float(np.max(abs(
        supplied_contribution + estimated_contribution - effective_j)
        * 4 / 10.0**sample_t[0] / arrays["off_Cv_x_scale"]))

    source_row_checks = {}
    for row_number, row_index in enumerate(row_indices):
        row_t = logt[row_index]
        row_q_mask = (logq >= 2.24) & (logq <= 2.76)
        row_q = logq[row_q_mask]
        row_temperature = np.full_like(row_q, row_t)
        row_fit = evaluate_fit(archive["thermo_0.1_coefficients"], fit_q, fit_t,
                               row_q * LN10, row_temperature * LN10)
        reference = state_from_native(row_fit, row_q, row_temperature)
        arrays[f"source_{row_number}_logQ"] = row_q
        arrays[f"source_{row_number}_logT"] = row_temperature
        arrays[f"source_{row_number}_logRho"] = row_q + 2 * row_t - 12
        for field in FIELDS:
            arrays[f"source_{row_number}_{field}"] = source[field][row_q_mask, row_index]
            arrays[f"source_{row_number}_reference_{field}"] = reference[field]

        line_j = compatibility(*(low[order, 0, face_indices[:, None] + np.array([-1, 0, 1]), row_index]
                                 for order in range(3)), step_q)
        arrays[f"row_{row_number}_J"] = line_j
        row_sides = {}
        for side, offset in (("left", -1), ("right", 0)):
            row_native = evaluate(nodes, complete, logq, logt, face_q,
                                  np.full_like(face_q, row_t),
                                  cell_indices=(face_indices + offset, np.full(len(face_q), lower_row)))
            row_sides[side] = state_from_native(row_native, face_q, np.full_like(face_q, row_t), face=True)
            arrays[f"row_{row_number}_{side}_Cv"] = row_sides[side]["Cv"]
            arrays[f"row_{row_number}_{side}_Cv_x"] = row_sides[side]["Cv_x"]
        scaled, scale = scaled_jump(row_sides["left"], row_sides["right"], "Cv_x")
        arrays[f"row_{row_number}_Cv_x_jump"] = row_sides["right"]["Cv_x"] - row_sides["left"]["Cv_x"]
        arrays[f"row_{row_number}_Cv_x_scaled_jump"] = scaled
        arrays[f"row_{row_number}_Cv_x_scale"] = scale
        source_row_checks[str(row_t)] = {
            "compatibility_vs_direct_max_scaled": float(np.max(abs(-4 * line_j / 10.0**row_t
                                                                   - arrays[f"row_{row_number}_Cv_x_jump"]) / scale)),
            "source_points": int(len(row_q))}

    large_index = int(np.argmax(abs(arrays["off_Cv_x_scaled_jump"])))
    quiet_index = int(np.argmin(abs(arrays["off_Cv_x_scaled_jump"])))
    metadata = {
        "task": "072 Figure 2 state-function and native-face diagnostics",
        "placement": {"job": os.environ.get("SLURM_JOB_ID"), "host": os.uname().nodename,
                      "cpus_on_node": os.environ.get("SLURM_CPUS_ON_NODE"),
                      "job_end_epoch": os.environ.get("SLURM_JOB_END_TIME"), "numerical_threads": 1},
        "convention": "Original036 unsplit total Helmholtz potential; original saved green/red coefficients; no refit",
        "reference": "Original036 thermo_0.1 C4 total potential, evaluated at each point's own logQ/logT; residual=(candidate-reference)/reference; no fitted adjustment",
        "sampling": {"logT": float(sample_t[0]), "logQ": [float(sample_q[0]), float(sample_q[-1])],
                     "dense_queries": len(sample_q), "faces": len(face_q)},
        "bracketing_logT": [float(logt[index]) for index in row_indices],
        "twist_stencil_logT_span": [float(logt[lower_row - 4]), float(logt[lower_row + 5])],
        "source_row_coordinate": "Each source row uses its own log10rho=log10Q+2log10T-12; no source point lies on Figure2 T",
        "reproduction_max_relative": reproduction,
        "face_continuity_max_relative": continuity,
        "effective_line_compatibility_vs_direct_max": effective_error,
        "native_third_jump_vs_physical_Cv_x_max_scaled": direct_native_error,
        "physical_shear_jump_relation_max_scaled": physical_shear_error,
        "tangential_basis_reconstruction_max_scaled": tangent_reconstruction_error,
        "supplied_estimated_reconstruction_max_scaled": supplied_estimated_reconstruction_error,
        "tangential_decomposition": "Fixed original jets: order0 uses supplied line values/slopes/curvatures; order1 also uses estimated f_uu_tau; order2 also uses estimated f_u_tautau and f_uu_tautau. Terms can cancel and depend on reference convention; not physical error fractions.",
        "source_row_checks": source_row_checks,
        "jump_normalization": "right-minus-left partial / max(abs(Cv_left),abs(Cv_right),abs(partial_left),abs(partial_right)); zero only for finite all-zero inputs",
        "large_face": {"index": large_index, "logQ": float(face_q[large_index]),
                       "logRho": float(face_q[large_index] + 2 * sample_t[0] - 12)},
        "quiet_face": {"index": quiet_index, "logQ": float(face_q[quiet_index]),
                       "logRho": float(face_q[quiet_index] + 2 * sample_t[0] - 12)},
        "input_sha256": inputs,
        "versions": {"numpy": np.__version__},
        "new_fit_calls": 0,
    }
    metadata["residual_summary"] = {
        name: {field: relative_summary((arrays[f"{name}_{field}"] - arrays[f"thermo_0.1_{field}"])
                                       / arrays[f"thermo_0.1_{field}"]) for field in FIELDS}
        for name in NAMES if name != "thermo_0.1"}
    metadata["face_jump_summary"] = {field: relative_summary(arrays[f"off_{field}_scaled_jump"])
                                     for field in ("Cv_x", "Cv_t")}
    for key in ("large_face", "quiet_face"):
        index = metadata[key]["index"]
        metadata[key]["off_Cv_x_scaled_jump"] = float(arrays["off_Cv_x_scaled_jump"][index])
        metadata[key]["off_Cv_x_jump"] = float(arrays["off_Cv_x_jump"][index])
        metadata[key]["row_Cv_x_scaled_jumps"] = [float(arrays[f"row_{row_number}_Cv_x_scaled_jump"][index])
                                                     for row_number in range(2)]
        metadata[key]["off_Cv_x_jump_tangential_components"] = [float(-4 * arrays[f"off_J_component_tau_order_{order}"][index] /
                                                                   10.0**sample_t[0]) for order in range(3)]
        metadata[key]["off_Cv_x_jump_supplied_estimated_components"] = [
            float(-4 * arrays[f"off_J_component_{part}"][index] / 10.0**sample_t[0])
            for part in ("supplied", "estimated")]
    metadata["one_sided_Cv_x_finite_difference"] = {}
    for key in ("large_face", "quiet_face"):
        index = metadata[key]["index"]
        metadata["one_sided_Cv_x_finite_difference"][key] = {}
        for step in (1e-4, 5e-5, 1e-5):
            metadata["one_sided_Cv_x_finite_difference"][key][str(step)] = {
                side: float((face_density_difference(
                    nodes, complete, logq, logt, face_q[index], sample_t[0],
                    face_indices[index], lower_row, side, step)
                    - face_sides[side]["Cv_x"][index]) / arrays["off_Cv_x_scale"][index])
                for side in ("left", "right")}
    metadata["source_row_residual_summary"] = {
        str(logt[row_index]): {
            field: relative_summary((arrays[f"source_{row_number}_{field}"]
                                     - arrays[f"source_{row_number}_reference_{field}"])
                                    / arrays[f"source_{row_number}_reference_{field}"])
            for field in FIELDS}
        for row_number, row_index in enumerate(row_indices)}
    np.savez_compressed(arrays_path, **arrays)
    metadata["output_sha256"] = {"diagnostics.npz": checksum(arrays_path)}
    plot(arrays, metadata)
    metadata["output_sha256"].update({path.name: checksum(path) for path in OUTPUT.glob("*.pdf")})
    metadata["output_sha256"].update({path.name: checksum(path) for path in OUTPUT.glob("*.png")})
    metadata["script_sha256"] = checksum(__file__)
    for name, expected in inputs.items():
        if checksum(WORKSPACE / name) != expected:
            raise AssertionError(f"Read-only input changed: {name}")
    results_path.write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"reproduction": reproduction, "continuity": continuity,
                      "face_jump_summary": metadata["face_jump_summary"],
                      "large_face": metadata["large_face"], "quiet_face": metadata["quiet_face"],
                      "effective_error": effective_error, "source_row_checks": source_row_checks}, indent=2))


if __name__ == "__main__":
    main()
