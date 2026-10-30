"""Private-cache and installed-library protection for API experiment drivers."""

import hashlib
import os
from pathlib import Path


CODE = Path(__file__).resolve().parent
WORKSPACE = CODE.parents[4]
MESA = WORKSPACE / "mesa"
EXPERIMENT = WORKSPACE / "freeeos_biquintic_progress/experiments/050_runtime_api"
SDK = WORKSPACE / "mesasdk-x86_64-linux-26.6.1"
MASSFILE = MESA / "eos/eosFreeEOS_builder/mass_frac_solar.txt"
INLIST = CODE.parent / "mesa_probe/inlist_probe"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def installed_inventory():
    cache = MESA / "data/eosDT_data/cache"
    return {"cache": {path.name: (path.stat().st_size, path.stat().st_mtime_ns)
                       for path in sorted(cache.glob("*.bin"))},
            "eos_library_sha256": sha256(MESA / "build/eos/lib/libeos.a")}


def private_environment(output):
    if os.environ.get("SLURM_JOB_NODELIST") != os.uname().nodename:
        raise RuntimeError("Must execute on the inherited allocated node")
    output = Path(output).resolve()
    cache_root, temp_root = output / "runtime_cache", output / "runtime_temp_cache"
    cache = cache_root / "eosDT_cache"
    cache.mkdir(parents=True, exist_ok=True)
    temp_root.mkdir(parents=True, exist_ok=True)
    for original in sorted((MESA / "data/eosDT_data/cache").glob("*.bin")):
        link = cache / original.name
        if not link.exists():
            link.symlink_to(original)
        if not link.is_symlink() or link.resolve() != original.resolve():
            raise RuntimeError("Private cache link replaced or misdirected: " + str(link))
    environment = dict(os.environ)
    environment.update(MESA_DIR=str(MESA), MESASDK_ROOT=str(SDK),
                       MESA_CACHES_DIR=str(cache_root), MESA_TEMP_CACHES_DIR=str(temp_root),
                       OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                       LD_LIBRARY_PATH=f"{SDK}/lib:{SDK}/lib64:" + environment.get("LD_LIBRARY_PATH", ""))
    return environment


def verify_private_cache(output):
    cache = Path(output) / "runtime_cache/eosDT_cache"
    replaced = [path.name for path in cache.glob("*.bin") if not path.is_symlink()]
    if replaced:
        raise RuntimeError("Private installed-cache links replaced: " + repr(replaced))
