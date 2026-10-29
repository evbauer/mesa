"""Run the default MESA mixture probe with exclusively experiment-owned caches."""

import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path


PROBE = Path(__file__).resolve().parent
WORKSPACE = PROBE.parents[4]
REPORT = WORKSPACE / 'freeeos_biquintic_progress/experiments/014_active_domain'
MASS_SHA = 'bd62922804e073a2a410480994b8a28e85b0bc53e469729d6ed5136ea9da10a6'


def cache_inventory(directory):
    return {path.name: {'size': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns}
            for path in sorted(directory.glob('*.bin'))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesa-dir', type=Path, default=WORKSPACE / 'mesa')
    parser.add_argument('--sdk-dir', type=Path, default=WORKSPACE / 'mesasdk-x86_64-linux-26.6.1')
    parser.add_argument('--inlist', type=Path, default=PROBE / 'inlist_probe')
    parser.add_argument('--massfile', type=Path, default=WORKSPACE / 'mesa/eos/eosFreeEOS_builder/mass_frac_solar.txt')
    args = parser.parse_args()
    mesa, sdk, inlist, massfile = (path.resolve() for path in
                                  [args.mesa_dir, args.sdk_dir, args.inlist, args.massfile])
    mass_sha = hashlib.sha256(massfile.read_bytes()).hexdigest()
    if mass_sha != MASS_SHA:
        raise ValueError('Mass pattern checksum differs from the declared scaled-solar pattern')
    source_cache = mesa / 'data/eosDT_data/cache'
    before = cache_inventory(source_cache)
    if not before:
        raise ValueError('No installed EOS binary caches to reuse')
    cache_root, temp_root = REPORT / 'runtime_cache', REPORT / 'runtime_temp_cache'
    cache_dir = cache_root / 'eosDT_cache'
    cache_dir.mkdir(parents=True, exist_ok=True)
    temp_root.mkdir(parents=True, exist_ok=True)
    for name in before:
        target = source_cache / name
        destination = cache_dir / name
        if destination.is_symlink():
            if destination.resolve() != target.resolve():
                raise ValueError(f'Unexpected cache link: {destination}')
        elif destination.exists():
            raise ValueError(f'Private cache regeneration detected: {destination}')
        else:
            destination.symlink_to(target)
    environment = dict(os.environ)
    environment.update(MESA_DIR=str(mesa), MESASDK_ROOT=str(sdk), OMP_NUM_THREADS='1',
                       MESA_CACHES_DIR=str(cache_root), MESA_TEMP_CACHES_DIR=str(temp_root),
                       LD_LIBRARY_PATH=f'{sdk}/lib:{sdk}/lib64:' + environment.get('LD_LIBRARY_PATH', ''))
    environment.pop('MESA_TEMP_CACHES_DISABLE', None)
    metadata = {'mesa_directory': str(mesa), 'sdk_directory': str(sdk), 'inlist': str(inlist),
                'inlist_sha256': hashlib.sha256(inlist.read_bytes()).hexdigest(), 'massfile': str(massfile),
                'massfile_sha256': mass_sha, 'input_X_Y_Z': [.7, .28, .02],
                'installed_cache_inventory_before': before, 'cache_root': str(cache_root),
                'temporary_cache_root': str(temp_root), 'OMP_NUM_THREADS': 1}
    metadata_path = REPORT / 'probe_metadata.json'
    metadata_path.write_text(json.dumps(metadata, indent=2) + '\n')
    command = [str(PROBE / 'build/mesa_probe'), str(mesa), str(inlist), str(massfile)]
    completed = subprocess.run(command, cwd=REPORT, env=environment)
    after = cache_inventory(source_cache)
    replaced = [path.name for path in cache_dir.glob('*.bin') if not path.is_symlink()]
    metadata.update(command=command, process_returncode=completed.returncode,
                    installed_cache_inventory_unchanged=before == after, private_replaced_cache_links=replaced)
    metadata_path.write_text(json.dumps(metadata, indent=2) + '\n')
    if before != after or replaced:
        raise RuntimeError('Cache safety check failed; see probe_metadata.json')
    raise SystemExit(completed.returncode)


if __name__ == '__main__':
    main()
