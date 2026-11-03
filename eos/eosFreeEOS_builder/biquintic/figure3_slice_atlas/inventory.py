"""Write SHA-256 inventory for the completed task 071 deliverables."""

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
EXPERIMENT = ROOT / "freeeos_biquintic_progress/experiments/071_figure3_slice_atlas"
REPORT = ROOT / "freeeos_biquintic_progress/progress_reports"
ARCHIVE066 = ROOT / "freeeos_biquintic_progress/experiments/066_extended_figure3"
ARCHIVE060 = ROOT / "freeeos_biquintic_progress/experiments/060_pure_freeeos_faces"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    selected = [ARCHIVE066 / "curves.npz", ARCHIVE066 / "results.json",
                ARCHIVE060 / "face_manifest.npz", ARCHIVE060 / "face_results.npz",
                ROOT / "mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data",
                ROOT / "mesa/data/eosDT_data/cache/mesa-FreeEOS_02z70x.bin"]
    selected += sorted(path for path in HERE.glob("*.py"))
    selected += sorted(path for path in EXPERIMENT.iterdir() if path.is_file() and path.name != "artifact_manifest.json")
    selected += [REPORT / f"2026-11-02_figure3_slice_atlas.{extension}" for extension in ("md", "tex", "pdf")]
    result = {str(path.relative_to(ROOT)): {"sha256": digest(path), "bytes": path.stat().st_size}
              for path in selected}
    (EXPERIMENT / "artifact_manifest.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"Wrote SHA-256s for {len(result)} files")


if __name__ == "__main__":
    main()
