"""Table-only Jacobian jump of MESA's standard eps-grav thermodynamic form."""

import argparse
import json
from pathlib import Path

import numpy as np

from freeeos_potential import LN10, checksum, construct_twists, evaluate, query_cells, read_table, source_constraints, thermodynamic_partials
from run_audit import distribution


def gravity_residual(state, delta_x, delta_t):
    temperature, cp, gradient = state["T"], state["Cp"], state["grad_ad"]
    chi_rho, chi_t = state["chiRho"], state["chiT"]
    first = 1 - gradient * chi_t
    second = gradient * chi_rho
    rates = first * delta_t - second * delta_x
    residual = -temperature * cp * rates
    result = {"value": residual}
    for suffix in ("_x", "_t"):
        first_partial = -state["grad_ad" + suffix] * chi_t - gradient * state["chiT" + suffix]
        second_partial = state["grad_ad" + suffix] * chi_rho + gradient * state["chiRho" + suffix]
        rate_partial = first_partial * delta_t - second_partial * delta_x
        rate_partial += first if suffix == "_t" else -second
        result[suffix] = -temperature * (state["Cp" + suffix] * rates + cp * rate_partial)
        if suffix == "_t":
            result[suffix] += residual
    return result


def main():
    workspace = Path(__file__).resolve().parents[4]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=workspace / "freeeos_biquintic_progress/experiments/030_gravity_control")
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    archive = workspace / "freeeos_biquintic_progress/experiments/014_active_domain/active_faces.npz"
    saved = np.load(archive)
    table = workspace / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data"
    _, logq, logt, rows = read_table(table)
    low, valid, _ = source_constraints(rows, logq, logt)
    nodes, complete, _ = construct_twists(low, valid, (logq[1] - logq[0]) * LN10, (logt[1] - logt[0]) * LN10)
    result = {"purpose": "Simplified table-only EOS contribution to standard epsgrav value/Jacobian at exact component faces",
              "formula": "epsgrav=-T*Cp*((1-grad_ad*chiT)*delta_lnT-grad_ad*chiRho*delta_lnrho)/dt",
              "derivative": "EOS current x,t varied with start x,t held fixed, hence derivatives of increments are1; dt set1 and cancels normalization",
              "not_measured": "Not MESA residual invocation, Newton solve, time centering, composition/latentheat, finalblend, fullenergy-equation or stellar performance",
              "scaling": "jump / max(1,abs(left/right residual derivative),T*left/rightCv); an energy-Jacobian scale, NOT response-function scale",
              "provenance": {"source_sha256": checksum(table), "archive_sha256": checksum(archive), "script_sha256": checksum(__file__),
                             "preprocessing_sha256": checksum(Path(__file__).with_name("freeeos_potential.py")),
                             "MESA_consumer": "star/private/eps_grav.f90:163", "consumer_sha256": checksum(workspace / "mesa/star/private/eps_grav.f90")}, "faces": {}}
    arrays = {}
    for axis, direction in ((0, "u"), (1, "tau")):
        probe = saved[f"{direction}_probe"]
        pure = (probe[:, 2] == 0) & (probe[:, 5] >= 1 - 1e-12)
        sample_q, sample_t = saved[f"{direction}_q"][pure], saved[f"{direction}_t"][pure]
        if axis == 0:
            density_index = np.searchsorted(logq, sample_q, side="left") - 1
            temperature_index = query_cells(logq, logt, sample_q, sample_t)[1]
            right_indices = density_index + 1, temperature_index
        else:
            density_index = query_cells(logq, logt, sample_q, sample_t)[0]
            temperature_index = np.searchsorted(logt, sample_t, side="left") - 1
            right_indices = density_index, temperature_index + 1
        left = thermodynamic_partials(evaluate(nodes, complete, logq, logt, sample_q, sample_t,
                                     cell_indices=(density_index, temperature_index)), sample_q, sample_t)
        right = thermodynamic_partials(evaluate(nodes, complete, logq, logt, sample_q, sample_t,
                                      cell_indices=right_indices), sample_q, sample_t)
        face_result = {"count": len(sample_q), "categories": "actual pure FreeEOS centers only, fraction>=1-1e-12", "trials": {}}
        for amplitude in (0., 1e-4, 1e-3, 1e-2, .1):
            for ratio in (0., 1., -1.):
                delta_t, delta_x = amplitude, ratio * amplitude
                first, second = gravity_residual(left, delta_x, delta_t), gravity_residual(right, delta_x, delta_t)
                label = f"deltaT_{amplitude}_deltaRho_over_deltaT_{ratio}"
                trial = {"delta_lnT": delta_t, "delta_lnrho": delta_x, "jacobian_jumps": {}}
                for suffix in ("_x", "_t"):
                    scale = np.maximum.reduce([abs(first[suffix]), abs(second[suffix]), left["T"] * left["Cv"], right["T"] * right["Cv"], np.ones(len(sample_q))])
                    jump = (first[suffix] - second[suffix]) / scale
                    worst = int(np.argmax(abs(jump)))
                    trial["jacobian_jumps"][suffix] = {**distribution(jump), "worst": {"logQ": float(sample_q[worst]), "logT": float(sample_t[worst]),
                                                           "signed_scaled_jump": float(jump[worst]), "left": float(first[suffix][worst]), "right": float(second[suffix][worst])}}
                    arrays[f"{direction}_{label}{suffix}"] = jump
                face_result["trials"][label] = trial
        result["faces"][direction] = face_result
        arrays[f"{direction}_q"], arrays[f"{direction}_t"] = sample_q, sample_t
    np.savez_compressed(output / "jacobian_jumps.npz", **arrays)
    (output / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({direction: {label: trial["jacobian_jumps"]["_t"] for label, trial in data["trials"].items() if trial["delta_lnrho"] == 0} for direction, data in result["faces"].items()}, indent=2))


if __name__ == "__main__":
    main()
