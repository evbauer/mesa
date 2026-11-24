#!/usr/bin/env bash
set -euo pipefail
script_directory=$(cd "$(dirname "$0")" && pwd)
workspace=$(cd "$script_directory/../../../../.." && pwd)
mesa_directory=${1:-"$workspace/mesa"}
sdk_directory=${2:-"$workspace/mesasdk-x86_64-linux-26.6.1"}
output=${3:?Supply a fresh private output directory; frozen experiments must not be overwritten}
mesa_directory=$(realpath "$mesa_directory")
sdk_directory=$(realpath "$sdk_directory")
output=$(realpath -m "$output")
if [[ "${SLURMD_NODENAME:-}" != "$(hostname -s)" || -z "${SLURM_JOB_ID:-}" ]]; then
    printf 'Run on the inherited allocated compute node.\n' >&2
    exit 1
fi
export MESASDK_ROOT="$sdk_directory"
export LD_LIBRARY_PATH="$sdk_directory/lib:$sdk_directory/lib64:${LD_LIBRARY_PATH:-}"
export PATH="$sdk_directory/bin:$PATH"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
dependency_paths=$(printf '%s:' "$output"/native_build/*/lib/pkgconfig)
export PKG_CONFIG_PATH="$dependency_paths$sdk_directory/lib/pkgconfig"
mkdir -p "$output/tests/modules"
read -r -a includes <<< "$(pkg-config --define-variable=math_slot=crmath --cflags mesa-eos mesa-const mesa-math mesa-chem mesa-utils)"
read -r -a libraries <<< "$(pkg-config --define-variable=math_slot=crmath --libs --static mesa-eos mesa-const mesa-math mesa-chem mesa-utils)"
flags=(-O2 -g -fopenmp -ffree-line-length-none -fcheck=all -ffp-contract=off -I"$output/native_build/eos/modules" -I"$output/tests/modules" -J"$output/tests/modules")
link_flags=(-Wl,--enable-new-dtags -L"$sdk_directory/lib" -L"$sdk_directory/lib64" -Wl,-rpath,"$sdk_directory/lib" -Wl,-rpath,"$sdk_directory/lib64")
# The independent source reader is linked only to explicit historical tests.
for source in freeeos_source historical_source_test; do
    gfortran "${flags[@]}" "${includes[@]}" -c "$script_directory/$source.f90" -o "$output/tests/$source.o"
done
for test_name in test_reader test_tables test_loader; do
    extra=()
    if [[ "$test_name" != test_reader ]]; then extra=("$output/tests/freeeos_source.o" "$output/tests/historical_source_test.o"); fi
    gfortran "${flags[@]}" "${includes[@]}" "$script_directory/$test_name.f90" -o "$output/tests/$test_name" "${link_flags[@]}" "${extra[@]}" "${libraries[@]}"
done
source_directory="$script_directory/.."
for driver in native_driver dual_handle_driver; do
    driver_source="$source_directory/native_integration/$driver.f90"
    [[ "$driver" == dual_handle_driver ]] && driver_source="$source_directory/co_table_support/$driver.f90"
    gfortran "${flags[@]}" "${includes[@]}" "$source_directory/runtime_bridge/bridge_support.f90" "$driver_source" \
        -o "$output/tests/$driver" "${link_flags[@]}" "${libraries[@]}"
done
# Production adapter tests use real public EOS handles and the standard loader.
for driver in test_legacy_tables test_legacy_family test_control; do
    gfortran "${flags[@]}" "${includes[@]}" "$source_directory/runtime_bridge/bridge_support.f90" "$script_directory/$driver.f90" \
        -o "$output/tests/$driver" "${link_flags[@]}" "${libraries[@]}"
done
