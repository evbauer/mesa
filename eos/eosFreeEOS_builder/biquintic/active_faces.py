"""Classify exact component face jumps using actual default MESA EOS fractions."""

import argparse
import io
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

from freeeos_potential import LN10, checksum, construct_twists, evaluate, read_table, source_constraints, thermodynamic_partials, valid_cells
from run_audit import distribution, own_use
from validate_table import material_nodes


RESPONSE_FIELDS = ("Cv", "Cp", "chiRho", "chiT", "gamma1", "gamma3", "grad_ad")


def classifications(status, fraction):
    successful = (status == 0) & np.isfinite(fraction)
    return {"all_successful": successful,
            "active": successful & (fraction > 1e-12),
            "pure_FreeEOS": successful & (fraction >= 1 - 1e-12),
            "blend": successful & (fraction > 1e-12) & (fraction < 1 - 1e-12),
            "inactive": successful & (fraction <= 1e-12),
            "failed": ~successful}


def face_coordinates(cells, logq, logt, axis, generator, maximum):
    adjacent = cells[:-1, :] & cells[1:, :] if axis == 0 else cells[:, :-1] & cells[:, 1:]
    density_index, temperature_index = np.nonzero(adjacent)
    choice = generator.choice(len(density_index), min(maximum, len(density_index)), replace=False)
    density_index, temperature_index = density_index[choice], temperature_index[choice]
    fraction = generator.uniform(.1, .9, len(choice))
    if axis == 0:
        sample_q = logq[density_index + 1]
        sample_t = logt[temperature_index] + fraction * (logt[1] - logt[0])
        right = density_index + 1, temperature_index
    else:
        sample_q = logq[density_index] + fraction * (logq[1] - logq[0])
        sample_t = logt[temperature_index + 1]
        right = density_index, temperature_index + 1
    selected = own_use(sample_q, sample_t)
    return sample_q[selected], sample_t[selected], (density_index[selected], temperature_index[selected]), (right[0][selected], right[1][selected])


def jump_summary(jump, selected, sample_q, sample_t, fraction):
    if not np.any(selected):
        return {"count": 0}
    indices = np.flatnonzero(selected)
    worst = indices[np.argmax(abs(jump[selected]))]
    return {"count": len(indices), **distribution(jump[selected]),
            "worst": {"logQ": float(sample_q[worst]), "logT": float(sample_t[worst]),
                      "logRho": float(sample_q[worst] + 2 * sample_t[worst] - 12),
                      "fraction_FreeEOS": float(fraction[worst]), "signed_scaled_jump": float(jump[worst])}}


def save_figure(output, samples):
    os.environ.setdefault("MPLCONFIGDIR", str(output / ".mplconfig"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    figure, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    norm = TwoSlopeNorm(vmin=-7, vcenter=-3, vmax=1)
    for column, (direction, data) in enumerate(samples.items()):
        sample_q, sample_t = data["q"], data["t"]
        log_density = sample_q + 2 * sample_t - 12
        active = data["masks"]["active"]
        values = np.log10(np.maximum(abs(data["jumps"]["Cv_t"]), 1e-17))
        for row, selected in enumerate((data["masks"]["all_successful"], active)):
            axis = axes[row, column]
            if row == 1:
                inactive = data["masks"]["inactive"]
                axis.scatter(log_density[inactive], sample_t[inactive], c=".85", s=2, rasterized=True)
            image = axis.scatter(log_density[selected], sample_t[selected], c=values[selected],
                                 s=3, cmap="coolwarm", norm=norm, rasterized=True)
            axis.set_title(f"{direction} faces: {'all own-use' if row == 0 else 'actual FreeEOS fraction > 1e-12'} ({np.count_nonzero(selected):,})", fontsize=9)
            axis.set_xlabel("log10 rho [g cm^-3]")
            axis.set_ylabel("log10 T [K]")
    figure.colorbar(image, ax=axes, label="log10 |jump in dCv/dlnT| / max(1, |Cv|, |one-sided partials|)", extend="both")
    figure.suptitle("Exact one-sided potential cell faces, off-node tangential coordinates\nScaled-solar X=.70, Z=.02; gray = MESA FreeEOS inactive; component jumps, NOT final-EOS errors", fontsize=10)
    for extension in ("pdf", "png"):
        figure.savefig(output / f"active_Cv_faces.{extension}", dpi=180)
    plt.close(figure)


def main():
    workspace = Path(__file__).resolve().parents[4]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesa", type=Path, default=workspace / "mesa")
    parser.add_argument("--output", type=Path, default=workspace / "freeeos_biquintic_progress/experiments/014_active_domain")
    parser.add_argument("--max-faces", type=int, default=12000)
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    table = arguments.mesa / "data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
    _, logq, logt, rows = read_table(table)
    low, valid, _ = source_constraints(rows, logq, logt)
    nodes, complete, _ = construct_twists(low, valid, (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10)
    cells = valid_cells(complete) & valid_cells(material_nodes(rows, logq, logt))
    generator = np.random.default_rng(20261029)
    samples = {}
    for axis, direction in ((0, "u"), (1, "tau")):
        sample_q, sample_t, left_indices, right_indices = face_coordinates(cells, logq, logt, axis, generator, arguments.max_faces)
        left = thermodynamic_partials(evaluate(nodes, complete, logq, logt, sample_q, sample_t, cell_indices=left_indices), sample_q, sample_t)
        right = thermodynamic_partials(evaluate(nodes, complete, logq, logt, sample_q, sample_t, cell_indices=right_indices), sample_q, sample_t)
        jumps = {}
        for field in RESPONSE_FIELDS:
            for suffix in ("_x", "_t"):
                scale = np.maximum.reduce([abs(left[field]), abs(right[field]), abs(left[field + suffix]), abs(right[field + suffix]), np.ones(len(sample_q))])
                jumps[field + suffix] = (left[field + suffix] - right[field + suffix]) / scale
        samples[direction] = {"q": sample_q, "t": sample_t, "jumps": jumps}
    query_q = np.concatenate([data["q"] for data in samples.values()])
    query_t = np.concatenate([data["t"] for data in samples.values()])
    query_text = "".join(f"{density:.17g} {temperature:.17g}\n" for density, temperature in zip(query_q + 2 * query_t - 12, query_t))
    (output / "face_queries.txt").write_text(query_text)
    command = [sys.executable, "-B", str(Path(__file__).parent / "mesa_probe/run_probe.py"), "--mesa-dir", str(arguments.mesa.resolve())]
    completed = subprocess.run(command, input=query_text, capture_output=True, text=True, timeout=600)
    (output / "face_probe.tsv").write_text(completed.stdout)
    (output / "face_probe.stderr").write_text(completed.stderr)
    if completed.returncode:
        raise RuntimeError(f"Probe failed with status {completed.returncode}")
    probe = np.loadtxt(io.StringIO(completed.stdout), skiprows=1, ndmin=2)
    if probe.shape != (len(query_q), 14):
        raise ValueError(f"Unexpected probe output shape {probe.shape}")
    if not np.allclose(probe[:, :2], np.column_stack((query_q + 2 * query_t - 12, query_t)), rtol=0, atol=1e-12):
        raise ValueError("Probe reordered query coordinates")
    result = {"purpose": "Reclassify the identical 002 material own-use face queries by actual default MESA blend fractions",
              "not_measured": "Potential is not installed in MESA; these are component-polynomial Jacobian jumps, not final blended EOS errors or source accuracy",
              "sampling": "same RNG20261029 and max12000 perdirection as002; one exact face coordinate and .1-.9 offnode tangent; both cells valid material corners and own-use cuts",
              "scaling": "signed left-minus-right response derivative / max(1,abs(left/right response),abs(left/right derivative)); x,t mean physical natural-log rho,T",
              "provenance": {"source_sha256": checksum(table), "script_sha256": checksum(__file__), "preprocessing_sha256": checksum(Path(__file__).parent / "freeeos_potential.py"), "probe_command": command,
                             "composition": "X=.70,Y=.28,Z=.02, normalized mass_frac_solar.txt metals, one representative isotope per element; WDM_AV false"}, "faces": {}}
    offset = 0
    archives = {}
    for direction, data in samples.items():
        count = len(data["q"])
        local_probe = probe[offset:offset + count]
        offset += count
        fractions = local_probe[:, 5]
        masks = classifications(local_probe[:, 2], fractions)
        data["masks"] = masks
        data["fraction"] = fractions
        categories = {}
        for category, selected in masks.items():
            category_data = {"count": int(np.count_nonzero(selected))}
            if category != "failed":
                category_data["jumps"] = {field: jump_summary(jump, selected, data["q"], data["t"], fractions) for field, jump in data["jumps"].items()}
            categories[category] = category_data
        result["faces"][direction] = categories
        archives.update({f"{direction}_q": data["q"], f"{direction}_t": data["t"], f"{direction}_probe": local_probe})
        archives.update({f"{direction}_{field}_jump": jump for field, jump in data["jumps"].items()})
    np.savez_compressed(output / "active_faces.npz", **archives)
    (output / "face_results.json").write_text(json.dumps(result, indent=2) + "\n")
    save_figure(output, samples)
    print(json.dumps({direction: {category: {"count": values["count"], **({"Cv_t": values["jumps"]["Cv_t"], "chiT_t": values["jumps"]["chiT_t"]} if values.get("count") and category != "failed" else {})} for category, values in categories.items()} for direction, categories in result["faces"].items()}, indent=2))


if __name__ == "__main__":
    main()
