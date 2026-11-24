#!/usr/bin/env bash
set -euo pipefail

script_directory=$(cd "$(dirname "$0")" && pwd)
workspace=$(cd "$script_directory/../../../../.." && pwd)
mesa_directory="$workspace/mesa"
sdk_directory="$workspace/mesasdk-x86_64-linux-26.6.1"
output="$workspace/freeeos_biquintic_progress/experiments/052_stellar_api"
private_eos="$workspace/freeeos_biquintic_progress/experiments/050_runtime_api/build/eos"
build_directory="../freeeos_biquintic_progress/experiments/052_stellar_api/build"

if [[ "${SLURM_JOB_NODELIST:-}" != "$(hostname -s)" ]]; then
    printf 'Run directly on the inherited allocated node.\n' >&2
    exit 1
fi
if (( ${SLURM_JOB_END_TIME:-0} - $(date +%s) < 1800 )); then
    printf 'At least thirty minutes of allocation time must remain.\n' >&2
    exit 1
fi
if [[ ! -f "$private_eos/lib/libeos.a" ]]; then
    printf 'Build and verify the private EOS first.\n' >&2
    exit 1
fi
export MESASDK_ROOT="$sdk_directory"
export PATH="$sdk_directory/bin:$PATH"
export LD_LIBRARY_PATH="$sdk_directory/lib:$sdk_directory/lib64:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$output/build"
sha256sum "$mesa_directory"/build/*/lib/*.a > "$output/installed_before.sha256"
for package_directory in "$mesa_directory"/build/*; do
    package_name=${package_directory##*/}
    [[ -d "$package_directory" ]] || continue
    destination="$output/build/$package_name"
    case "$package_name" in
        rates|net|kap|atm|star_data|star)
            if [[ -L "$destination" ]]; then
                printf 'Refusing mutable package symlink: %s\n' "$destination" >&2
                exit 1
            fi
            mkdir -p "$destination"
            ;;
        *)
            source="$package_directory"
            [[ "$package_name" != eos ]] || source="$private_eos"
            if [[ ! -e "$destination" && ! -L "$destination" ]]; then
                ln -s "$source" "$destination"
            fi
            if [[ ! -L "$destination" || "$(readlink -f "$destination")" != "$(readlink -f "$source")" ]]; then
                printf 'Unexpected dependency link: %s\n' "$destination" >&2
                exit 1
            fi
            ;;
    esac
done
dependency_paths=$(printf '%s:' "$output"/build/*/lib/pkgconfig)
export PKG_CONFIG_PATH="$dependency_paths$sdk_directory/lib/pkgconfig"
for package_name in rates kap net atm star_data star; do
    printf '%s build %s\n' "$(date --iso-8601=seconds)" "$package_name"
    make -C "$mesa_directory/$package_name" -j16 BUILD_DIR="$build_directory" all > "$output/build_$package_name.log" 2>&1
done
sha256sum -c "$output/installed_before.sha256" > "$output/installed_after_check.txt"
printf 'Private coherent EOS dependents built at %s/build\n' "$output"
