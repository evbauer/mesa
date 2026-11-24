#!/usr/bin/env bash
set -euo pipefail
script_directory=$(cd "$(dirname "$0")" && pwd)
workspace=$(cd "$script_directory/../../../../.." && pwd)
mesa_directory=${1:-"$workspace/mesa"}
sdk_directory=${2:-"$workspace/mesasdk-x86_64-linux-26.6.1"}
output=${3:?Supply a fresh private output directory; frozen experiments must not be overwritten}
stage=${4:-eos}
build_jobs=${MESA_BUILD_JOBS:-16}
if [[ ! "$build_jobs" =~ ^[1-9][0-9]*$ ]]; then
    printf 'MESA_BUILD_JOBS must be a positive integer.\n' >&2
    exit 2
fi
case "$stage" in eos|stellar) ;; *) exit 2 ;; esac
mesa_directory=$(realpath "$mesa_directory")
sdk_directory=$(realpath "$sdk_directory")
output=$(realpath -m "$output")
if [[ "${SLURMD_NODENAME:-}" != "$(hostname -s)" || -z "${SLURM_JOB_ID:-}" ]]; then
    printf 'Run inside the inherited allocation on its compute node.\n' >&2
    exit 1
fi
export MESASDK_ROOT="$sdk_directory"
export PATH="$sdk_directory/bin:$PATH"
export LD_LIBRARY_PATH="$sdk_directory/lib:$sdk_directory/lib64:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
mkdir -p "$output/native_build/eos"
if [[ ! -f "$output/installed_before.sha256" ]]; then
    sha256sum "$mesa_directory"/build/*/lib/*.a > "$output/installed_before.sha256"
fi
sha256sum -c "$output/installed_before.sha256" > "$output/installed_prebuild_check.txt"
for package_directory in "$mesa_directory"/build/*; do
    [[ -d "$package_directory" ]] || continue
    package_name=${package_directory##*/}
    [[ "$package_name" == eos ]] && continue
    destination="$output/native_build/$package_name"
    case "$package_name" in
        rates|kap|net|atm|star_data|star)
            if [[ "$stage" == stellar ]]; then
                if [[ -L "$destination" ]]; then
                    if [[ "$(readlink -f "$destination")" != "$(readlink -f "$package_directory")" ]]; then
                        printf 'Unexpected dependency link: %s\n' "$destination" >&2
                        exit 1
                    fi
                    unlink "$destination"
                fi
                mkdir -p "$destination"
                continue
            elif [[ -d "$destination" && ! -L "$destination" ]]; then
                continue
            fi
            ;;
    esac
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
build_directory=$(realpath --relative-to="$mesa_directory" "$output/native_build")
packages=(eos)
if [[ "$stage" == stellar ]]; then packages=(rates kap net atm star_data star); fi
for package_name in "${packages[@]}"; do
    make -C "$mesa_directory/$package_name" -j"$build_jobs" BUILD_DIR="$build_directory" all > "$output/native_build_$package_name.log" 2>&1
done
sha256sum -c "$output/installed_before.sha256" > "$output/installed_postbuild_check.txt"
printf 'Private %s build ready: %s/native_build\n' "$stage" "$output"
