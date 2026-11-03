"""Summarize coverage, mask causes, and saturation from archived arrays."""

import json

import numpy as np

from common import OUTPUT, digest, save_json


def range_of(mask, density, temperature, ye):
    positions = np.argwhere(mask)
    if not len(positions):
        return None
    log_density = np.log10(density[positions[:, 0]] / ye)
    log_temperature = np.log10(temperature[positions[:, 1]])
    return {"logRho": [float(log_density.min()), float(log_density.max())],
            "logT": [float(log_temperature.min()), float(log_temperature.max())]}


def main():
    ye = float(np.loadtxt(OUTPUT / "composition.tsv", skiprows=1)[5])
    report = {"HELM": {}, "FreeEOS": {}}
    for direction in ("T", "d"):
        filename = "helm_T_traces_corrected.npz" if direction == "T" else "helm_d_traces.npz"
        with np.load(OUTPUT / f"helm_{direction}_grid.npz") as grid, np.load(OUTPUT / filename) as data:
            external = data["route_success"] & data["fraction_good"] & ~data["pure"]
            masks = {"fully_FreeEOS_direct_probe": data["pure"],
                     "HELM_evaluation_failed_external": external & ~data["kernel_good"],
                     "routing_failed_or_bad_fractions": ~data["route_success"] | ~data["fraction_good"],
                     "HELM_Coulomb_branch_near": external & data["branch_near"],
                     "plotted_HELM": data["shown"]}
            if direction == "T":
                masks["HELM_precision_unresolved"] = external & data["kernel_good"] & ~data["resolved"]
            report["HELM"][direction] = {key: {"count": int(mask.sum()),
                                                "range": range_of(mask, grid["din"], grid["temperature"], ye)}
                                          for key, mask in masks.items()}
            report["HELM"][direction]["shape"] = list(data["shown"].shape)
            report["HELM"][direction]["native_face_density_count"] = len(grid["din"])
            report["HELM"][direction]["native_face_temperature_count"] = len(grid["temperature"])
            report["HELM"][direction]["source_sha256"] = digest(OUTPUT / filename)
    with np.load(OUTPUT / "freeeos_traces.npz") as data:
        for direction in ("u", "tau"):
            pure = data[f"{direction}_pure"]
            support = data[f"{direction}_support_code"]
            result = {"fully_FreeEOS": int(pure.sum()),
                      "supported": int((pure & (support == 0)).sum()),
                      "gap_by_060_cause": {str(code): int((pure & (support == code)).sum())
                                           for code in range(1, 5)}}
            report["FreeEOS"][direction] = result
    report["normalization"] = json.loads((OUTPUT / "freeeos_summary.json").read_text())
    save_json(OUTPUT / "coverage_diagnostics.json", report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
