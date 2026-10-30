"""Run unchanged047 adapter tests in the explicitly isolated059 directory."""

from pathlib import Path
import sys

DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(DIRECTORY.parent))
from runtime_adapter import run_tests


def main():
    output = DIRECTORY.parents[4] / "freeeos_biquintic_progress/experiments/059_arithmetic_ad"
    if (output / "adapter_validation").exists():
        raise RuntimeError("Preserve059 adapter validation")
    (output / "solar.identity").write_bytes((run_tests.EXPERIMENT / "solar.identity").read_bytes())
    run_tests.EXPERIMENT = output
    sys.argv = [str(Path(run_tests.__file__)), "--driver", str(output / "build/driver/runtime_driver"),
                "--output", str(output / "adapter_validation")]
    run_tests.main()


if __name__ == "__main__":
    main()
