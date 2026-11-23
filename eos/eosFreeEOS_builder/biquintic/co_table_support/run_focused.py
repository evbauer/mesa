from pathlib import Path
import os
import subprocess
import sys


WORKSPACE = Path(__file__).resolve().parents[5]
EXPERIMENT = WORKSPACE / "freeeos_biquintic_progress/experiments/075_co_table_support"
BUILD = EXPERIMENT / "build"
FAMILY = WORKSPACE / "freeeos_biquintic_progress/experiments/062_family_cache/family"
SDK = WORKSPACE / "mesasdk-x86_64-linux-26.6.1"


def run(name, executable, directory, mode):
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "4" if mode == "parallel" else "1"
    environment["OPENBLAS_NUM_THREADS"] = "1"
    environment["MKL_NUM_THREADS"] = "1"
    environment["LD_LIBRARY_PATH"] = f"{SDK}/lib:{SDK}/lib64:" + environment.get("LD_LIBRARY_PATH", "")
    command = [str(BUILD / executable), str(directory), mode]
    result = subprocess.run(command, env=environment, text=True, capture_output=True)
    (EXPERIMENT / f"{name}.log").write_text(
        f"command: {' '.join(command)}\nexit: {result.returncode}\n{result.stdout}{result.stderr}"
    )
    print(name, result.returncode, result.stdout.strip())
    if result.returncode:
        raise SystemExit(result.returncode)


def fixture(name):
    directory = EXPERIMENT / "fixtures" / name
    directory.mkdir(parents=True, exist_ok=True)
    lines = (FAMILY / "family.manifest").read_text().splitlines()
    target = next(position for position, line in enumerate(lines) if "mesa-FreeEOS_80z00x_CO.bq" in line)
    if name == "missing":
        lines.pop(target)
    elif name == "duplicate":
        lines.append(lines[target])
    elif name == "unknown":
        lines[target] = lines[target].replace("_CO.bq", "_XX.bq")
    elif name == "unavailable":
        tokens = lines[target].split()
        tokens[2] = "0"
        lines[target] = " ".join(tokens)
    elif name == "absent_payload":
        pass
    lines[1] = str(len(lines) - 2)
    (directory / "family.manifest").write_text("\n".join(lines) + "\n")
    for source in FAMILY.glob("*.bq"):
        if name == "absent_payload" and source.name == "mesa-FreeEOS_80z00x_CO.bq":
            continue
        link = directory / source.name
        if not link.exists():
            link.symlink_to(source)
    return directory


def main():
    EXPERIMENT.mkdir(parents=True, exist_ok=True)
    for mode in ("good", "parallel", "highrad", "nonideal", "duplicate", "badsha",
                 "unavailable", "missing", "invalid", "nonphysical"):
        directory = EXPERIMENT / "legacy_fixtures" / mode
        directory.mkdir(parents=True, exist_ok=True)
        run(f"legacy_{mode}", "test_family", directory, mode)
    for mode in ("", "parallel", "unsupported"):
        run(f"real_{mode or 'configured'}", "test_co", FAMILY, mode)
    for mode in ("missing", "duplicate", "unknown", "unavailable", "absent_payload"):
        run(f"variant_{mode}", "test_co", fixture(mode), mode)


if __name__ == "__main__":
    main()
