#!/usr/bin/env bash
set -euo pipefail

script_directory=$(cd "$(dirname "$0")" && pwd)
workspace=$(cd "$script_directory/../../../../.." && pwd)
mesa_directory="$workspace/mesa"
sdk_directory="$workspace/mesasdk-x86_64-linux-26.6.1"
output="$workspace/freeeos_biquintic_progress/experiments/050_runtime_api"
build_directory="../freeeos_biquintic_progress/experiments/050_runtime_api/build"

if [[ "${SLURM_JOB_NODELIST:-}" != "$(hostname -s)" ]]; then
    printf 'Run directly on the inherited allocated node.\n' >&2
    exit 1
fi
export MESASDK_ROOT="$sdk_directory"
export PATH="$sdk_directory/bin:$PATH"
export LD_LIBRARY_PATH="$sdk_directory/lib:$sdk_directory/lib64:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
dependency_paths=$(printf '%s:' "$mesa_directory"/build/*/lib/pkgconfig)
export PKG_CONFIG_PATH="$dependency_paths$sdk_directory/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
mkdir -p "$output"
mkdir -p "$output/build"
for package_directory in "$mesa_directory"/build/*; do
    package_name=${package_directory##*/}
    if [[ -d "$package_directory" && "$package_name" != eos && ! -e "$output/build/$package_name" ]]; then
        ln -s "$package_directory" "$output/build/$package_name"
    fi
done
make -C "$mesa_directory/eos" -j16 BUILD_DIR="$build_directory" all > "$output/build_eos.log" 2>&1
printf 'Private EOS build: %s/build/eos\n' "$output"
