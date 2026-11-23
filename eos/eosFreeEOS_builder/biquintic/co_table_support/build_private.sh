#!/usr/bin/env bash
set -euo pipefail

script_directory=$(cd "$(dirname "$0")" && pwd)
workspace=$(cd "$script_directory/../../../../.." && pwd)
mesa_directory="$workspace/mesa"
sdk_directory="$workspace/mesasdk-x86_64-linux-26.6.1"
output="$workspace/freeeos_biquintic_progress/experiments/075_co_table_support"
build_directory="../freeeos_biquintic_progress/experiments/075_co_table_support/native_build"

if [[ "${SLURM_JOB_NODELIST:-}" != "$(hostname -s)" ]]; then
    printf 'Run on the inherited single allocated node.\n' >&2
    exit 1
fi
export MESASDK_ROOT="$sdk_directory"
export PATH="$sdk_directory/bin:$PATH"
export LD_LIBRARY_PATH="$sdk_directory/lib:$sdk_directory/lib64:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$output/native_build/eos"
sha256sum "$mesa_directory"/build/*/lib/*.a > "$output/installed_before.sha256"
for package_directory in "$mesa_directory"/build/*; do
    [[ -d "$package_directory" ]] || continue
    package_name=${package_directory##*/}
    [[ "$package_name" == eos ]] && continue
    destination="$output/native_build/$package_name"
    if [[ ! -e "$destination" && ! -L "$destination" ]]; then
        ln -s "$package_directory" "$destination"
    fi
    if [[ ! -L "$destination" || "$(readlink -f "$destination")" != "$(readlink -f "$package_directory")" ]]; then
        printf 'Unexpected dependency link: %s\n' "$destination" >&2
        exit 1
    fi
done
dependency_paths=$(printf '%s:' "$output"/native_build/*/lib/pkgconfig)
export PKG_CONFIG_PATH="$dependency_paths$sdk_directory/lib/pkgconfig"
make -C "$mesa_directory/eos" -j12 BUILD_DIR="$build_directory" all > "$output/native_build_eos.log" 2>&1
sha256sum -c "$output/installed_before.sha256" > "$output/installed_postbuild_check.txt"
printf 'Private EOS build ready: %s/native_build/eos\n' "$output"
