#!/usr/bin/env bash
set -euo pipefail

script_directory=$(cd "$(dirname "$0")" && pwd)
workspace=$(cd "$script_directory/../../../../.." && pwd)
mesa_directory="$workspace/mesa"
sdk_directory="$workspace/mesasdk-x86_64-linux-26.6.1"
output="$workspace/freeeos_biquintic_progress/experiments/065_native_integration"
build_directory="../freeeos_biquintic_progress/experiments/065_native_integration/build"
stage=${1:-all}

if [[ "${SLURM_JOB_NODELIST:-}" != "$(hostname -s)" ]]; then
    printf 'Run on the inherited single allocated node.\n' >&2
    exit 1
fi
if (( ${SLURM_JOB_END_TIME:-0} - $(date +%s) < 1800 )); then
    printf 'Reserve at least thirty minutes before starting a build.\n' >&2
    exit 1
fi
case "$stage" in eos|stellar|all) ;; *) exit 2 ;; esac
export MESASDK_ROOT="$sdk_directory"
export PATH="$sdk_directory/bin:$PATH"
export LD_LIBRARY_PATH="$sdk_directory/lib:$sdk_directory/lib64:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$output/build"
if [[ ! -f "$output/installed_before.sha256" ]]; then
    sha256sum "$mesa_directory"/build/*/lib/*.a > "$output/installed_before.sha256"
fi
sha256sum -c "$output/installed_before.sha256" > "$output/installed_prebuild_check.txt"
for package_directory in "$mesa_directory"/build/*; do
    [[ -d "$package_directory" ]] || continue
    package_name=${package_directory##*/}
    destination="$output/build/$package_name"
    case "$package_name" in
        eos|rates|kap|net|atm|star_data|star)
            if [[ -L "$destination" ]]; then
                printf 'Refusing mutable package symlink: %s\n' "$destination" >&2
                exit 1
            fi
            mkdir -p "$destination"
            ;;
        *)
            if [[ ! -e "$destination" && ! -L "$destination" ]]; then
                ln -s "$package_directory" "$destination"
            fi
            if [[ ! -L "$destination" || "$(readlink -f "$destination")" != "$(readlink -f "$package_directory")" ]]; then
                printf 'Unexpected dependency link: %s\n' "$destination" >&2
                exit 1
            fi
            ;;
    esac
done
dependency_paths=$(printf '%s:' "$output"/build/*/lib/pkgconfig)
export PKG_CONFIG_PATH="$dependency_paths$sdk_directory/lib/pkgconfig"
packages=()
if [[ "$stage" != stellar ]]; then packages+=(eos); fi
if [[ "$stage" != eos ]]; then packages+=(rates kap net atm star_data star); fi
for package_name in "${packages[@]}"; do
    timestamp=$(date +%Y%m%dT%H%M%S)
    logfile="$output/build_${package_name}_${timestamp}.log"
    printf '%s build %s log=%s\n' "$(date --iso-8601=seconds)" "$package_name" "$logfile"
    make -C "$mesa_directory/$package_name" -j16 BUILD_DIR="$build_directory" all > "$logfile" 2>&1
done
sha256sum -c "$output/installed_before.sha256" > "$output/installed_postbuild_check.txt"
printf 'Private build ready: %s/build\n' "$output"
