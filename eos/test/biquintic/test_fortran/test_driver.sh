#!/usr/bin/env bash
set -euo pipefail

driver=${1:-./cell_driver}
scratch=${2:-build/driver_smoke}
cases=0
mkdir -p "$scratch"

fail() {
    printf 'FAIL driver smoke: %s\n' "$*" >&2
    exit 1
}

check_output() {
    local output=$1 expected_cells=$2 mode=$3
    awk -v cells="$expected_cells" -v mode="$mode" '
        BEGIN {
            header = "f f_u f_tau f_uu f_u_tau f_tau_tau f_uuu f_uu_tau f_u_tau_tau f_tau_tau_tau"
            scale = 2.5e12
        }
        NR % 3 == 1 {
            if ($0 != header) exit 1
            next
        }
        {
            if (NF != 10) exit 1
            cell = int((NR + 2) / 3)
            affine = mode == "affine" || (mode == "mixed" && cell % 2 == 0)
            for (component = 1; component <= 10; component++) {
                expected = 0
                if (component == 1) expected = affine ? -1.5e12 : scale
                if (affine && component == 2) expected = 5e12
                if (affine && component == 3) expected = -7.5e12
                if ($component !~ /^[+-]?[0-9]+[.][0-9]+E[+-][0-9]+$/) exit 1
                error = ($component - expected) / scale
                if (error < 0) error = -error
                if (error > 5e-12) exit 1
            }
        }
        END { if (NR != 3 * cells) exit 1 }
    ' "$output" || fail "unexpected numerical output: $output"
}

expect_success() {
    local name=$1 input=$2 expected_cells=$3 mode=$4
    "$driver" < "$input" > "$scratch/$name.stdout" 2> "$scratch/$name.stderr" || fail "$name failed"
    [[ ! -s "$scratch/$name.stderr" ]] || fail "$name wrote unexpected stderr"
    check_output "$scratch/$name.stdout" "$expected_cells" "$mode"
    cases=$((cases + 1))
}

expect_failure() {
    local name=$1 input=$2 expected_status=$3 expected_cells=${4:-0} diagnostic=${5:-}
    local actual_status
    if "$driver" < "$input" > "$scratch/$name.stdout" 2> "$scratch/$name.stderr"; then
        fail "$name unexpectedly succeeded"
    else
        actual_status=$?
    fi
    [[ "$actual_status" == "$expected_status" ]] || fail "$name returned $actual_status instead of $expected_status"
    if [[ "$expected_status" == 2 ]]; then
        rg -q 'Invalid or incomplete cell input' "$scratch/$name.stderr" || fail "$name lacks parse diagnostic"
    else
        rg -q "Cell rejected; cell=$((expected_cells + 1)); status=$diagnostic" "$scratch/$name.stderr" || fail "$name lacks status diagnostic"
    fi
    check_output "$scratch/$name.stdout" "$expected_cells" constant
    cases=$((cases + 1))
}

expect_success constant constant_cell.txt 1 constant
expect_success affine affine_cell.txt 1 affine
cat constant_cell.txt affine_cell.txt constant_cell.txt > "$scratch/mixed.txt"
expect_success mixed "$scratch/mixed.txt" 3 mixed
for ((cell = 0; cell < 128; cell++)); do
    cat constant_cell.txt affine_cell.txt
done > "$scratch/batch.txt"
expect_success batch "$scratch/batch.txt" 256 mixed
printf '' > "$scratch/empty.txt"
expect_success empty "$scratch/empty.txt" 0 constant
printf '%s' "$(cat constant_cell.txt)" > "$scratch/no_final_newline.txt"
expect_success no_final_newline "$scratch/no_final_newline.txt" 1 constant
sed 's/ /\t/g' constant_cell.txt > "$scratch/tabs.txt"
expect_success tabs "$scratch/tabs.txt" 1 constant
sed 's/$/\r/' constant_cell.txt > "$scratch/crlf.txt"
expect_success crlf "$scratch/crlf.txt" 1 constant
sed 's/e12/D12/g' constant_cell.txt > "$scratch/d_exponents.txt"
expect_success d_exponents "$scratch/d_exponents.txt" 1 constant

coproc STREAM_DRIVER { "$driver" 2> "$scratch/persistent.stderr"; }
stream_pid=$STREAM_DRIVER_PID
stream_read=${STREAM_DRIVER[0]}
stream_write=${STREAM_DRIVER[1]}
printf '' > "$scratch/persistent.stdout"
for input in constant_cell.txt affine_cell.txt; do
    cat "$input" >& "$stream_write"
    for ((line_number = 0; line_number < 3; line_number++)); do
        if ! IFS= read -r -t 5 -u "$stream_read" output_line; then
            kill "$stream_pid" 2>/dev/null || true
            wait "$stream_pid" 2>/dev/null || true
            fail 'persistent pipe did not flush three lines before EOF'
        fi
        printf '%s\n' "$output_line" >> "$scratch/persistent.stdout"
    done
done
exec {stream_write}>&-
wait "$stream_pid" || fail 'persistent pipe driver failed'
[[ ! -s "$scratch/persistent.stderr" ]] || fail 'persistent pipe wrote unexpected stderr'
check_output "$scratch/persistent.stdout" 2 mixed
cases=$((cases + 1))

for ((records = 1; records <= 5; records++)); do
    head -n "$records" constant_cell.txt > "$scratch/truncated_$records.txt"
    expect_failure "truncated_$records" "$scratch/truncated_$records.txt" 2
done
for ((record = 1; record <= 6; record++)); do
    sed "${record}s/.*/0/" constant_cell.txt > "$scratch/partial_$record.txt"
    expect_failure "partial_$record" "$scratch/partial_$record.txt" 2
done
sed '3s/.*/2.5e12 0 0 0 0 0 0 0 0 0/' constant_cell.txt > "$scratch/extra_token.txt"
expect_failure extra_token "$scratch/extra_token.txt" 2
sed '3s/.*/2.5e12 0 0 \/ 0 0 0 0 0/' constant_cell.txt > "$scratch/slash.txt"
expect_failure slash "$scratch/slash.txt" 2
sed '3s/.*/2.5e12 0 0 0 0 0 0 0 2*0/' constant_cell.txt > "$scratch/repetition.txt"
expect_failure repetition "$scratch/repetition.txt" 2
sed '3s/.*/2.5e12 0 0 0 0 0 0 0 ,/' constant_cell.txt > "$scratch/null.txt"
expect_failure null "$scratch/null.txt" 2
sed '1s/-0.7/not_a_number/' constant_cell.txt > "$scratch/malformed.txt"
expect_failure malformed "$scratch/malformed.txt" 2
sed '3s/.*/ /' constant_cell.txt > "$scratch/blank.txt"
expect_failure blank "$scratch/blank.txt" 2
awk 'BEGIN { for (position = 0; position < 5000; position++) printf "0 "; printf "\n" }' > "$scratch/oversized.txt"
expect_failure oversized "$scratch/oversized.txt" 2

queries=('-0.01 0.8' '1.01 0.8' '0.4 -0.01' '0.4 2.01')
for query in "${!queries[@]}"; do
    sed "2s/.*/${queries[$query]}/" affine_cell.txt > "$scratch/outside_$query.txt"
    expect_failure "outside_$query" "$scratch/outside_$query.txt" 1 0 3
done
for spacing in '0 0 0 2' '1 0 0 2' '0 1 0 0' '0 1 2 0'; do
    sed "1s/.*/$spacing/" affine_cell.txt > "$scratch/spacing.txt"
    expect_failure "spacing_$cases" "$scratch/spacing.txt" 1 0 2
done
for spacing in '0 1e-200 0 2' '0 1e200 0 2'; do
    sed "1s/.*/$spacing/" affine_cell.txt > "$scratch/unsafe_spacing.txt"
    expect_failure "unsafe_spacing_$cases" "$scratch/unsafe_spacing.txt" 1 0 2
done
for value in NaN Infinity -Infinity; do
    sed "3s/2.5e12/$value/" constant_cell.txt > "$scratch/nonfinite_corner.txt"
    expect_failure "nonfinite_corner_$cases" "$scratch/nonfinite_corner.txt" 1 0 1
    sed "2s/.*/$value 0.8/" constant_cell.txt > "$scratch/nonfinite_query.txt"
    expect_failure "nonfinite_query_$cases" "$scratch/nonfinite_query.txt" 1 0 1
    sed "1s/.*/$value 1 0 2/" constant_cell.txt > "$scratch/nonfinite_endpoint.txt"
    expect_failure "nonfinite_endpoint_$cases" "$scratch/nonfinite_endpoint.txt" 1 0 1
done
awk 'BEGIN {
    print "0 1 0 2"
    print "0.4 0.8"
    for (corner = 0; corner < 4; corner++) {
        for (component = 0; component < 9; component++) printf "1.7976931348623157e308 "
        printf "\n"
    }
}' > "$scratch/arithmetic_failure.txt"
expect_failure arithmetic_failure "$scratch/arithmetic_failure.txt" 1 0 4
cat constant_cell.txt "$scratch/outside_0.txt" > "$scratch/valid_then_outside.txt"
expect_failure valid_then_outside "$scratch/valid_then_outside.txt" 1 1 3
cat constant_cell.txt "$scratch/truncated_3.txt" > "$scratch/valid_then_truncated.txt"
expect_failure valid_then_truncated "$scratch/valid_then_truncated.txt" 2 1
{ cat constant_cell.txt; printf '0 1'; } > "$scratch/valid_then_partial_header.txt"
expect_failure valid_then_partial_header "$scratch/valid_then_partial_header.txt" 2 1
{ cat constant_cell.txt; printf 'garbage\n'; } > "$scratch/valid_then_garbage.txt"
expect_failure valid_then_garbage "$scratch/valid_then_garbage.txt" 2 1

printf 'PASS driver smoke cases: %d; includes a 256-cell stream, strict EOF/truncation and no extrapolation\n' "$cases"
