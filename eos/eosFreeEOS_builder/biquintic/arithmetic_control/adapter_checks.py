"""Reuse the unchanged047 checks with an explicit private053 output boundary."""

from pathlib import Path
import sys

DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(DIRECTORY.parent))
from runtime_adapter import run_tests


def main():
    output = DIRECTORY.parents[4] / "freeeos_biquintic_progress/experiments/053_arithmetic_control"
    original = run_tests.EXPERIMENT
    if (output / "adapter_validation").exists():
        raise RuntimeError("Preserve existing053 adapter checks")
    (output / "solar.identity").write_bytes((original / "solar.identity").read_bytes())
    run_tests.EXPERIMENT = output
    sys.argv = [str(Path(run_tests.__file__)), "--driver", str(output / "build/driver/runtime_driver"),
                "--output", str(output / "adapter_validation")]
    run_tests.main()


if __name__ == "__main__":
    main()
