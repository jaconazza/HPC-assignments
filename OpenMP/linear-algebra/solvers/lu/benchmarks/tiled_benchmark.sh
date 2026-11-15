#!/usr/bin/env bash

EXERCISE="lu_optimized.c"
EXE="./lu_acc"
EXTRA_FLAGS="-fno-inline -DPOLYBENCH_TIME"

DATASETS=(MINI_DATASET SMALL_DATASET STANDARD_DATASET LARGE_DATASET EXTRALARGE_DATASET)
NTHREADS=(4)
BLOCKS=(8 16 32 64 128)

VALGRIND_MIN_BLOCK=64

OUT_CSV="results.csv"
PERF_EVENTS="task-clock,cycles,instructions,cache-references,cache-misses,branch-misses"

HEADER="exec_time_s,tool,dataset,nthreads,block,rep,wall_seconds,task_clock_ms,cycles,instructions,cache_references,cache_misses,branch_misses,vg_ir,vg_i1mr,vg_ilmr,vg_dr,vg_d1mr,vg_dlmr,vg_dw,vg_d1mw,vg_dlmw,exit_status"

perf_reps() {
    case "$1" in
        LARGE_DATASET|EXTRALARGE_DATASET) echo 3 ;;
        *) echo 5 ;;
    esac
}

use_valgrind() {
    case "$1" in
        MINI_DATASET|SMALL_DATASET|STANDARD_DATASET) return 0 ;;
        *) return 1 ;;
    esac
}

write_row() {
    local IFS=,
    echo "$*" >> "$OUT_CSV"
}

elapsed() {
    awk -v s="$1" -v e="$2" 'BEGIN { printf "%.3f", e - s }'
}

build() {
    local dataset="$1" threads="$2" block="$3"
    local flags="-D${dataset} -DNTHREADS=${threads} -DBLOCK=${block} ${EXTRA_FLAGS}"
    local out
    if ! out=$(make "EXERCISE=${EXERCISE}" "EXT_CFLAGS=${flags}" clean all 2>&1); then
        echo "$out"
        echo "build failed" >&2
        exit 1
    fi
}

polybench_time() {
    grep -oE '[0-9]+\.[0-9]+' <<< "$1" | tail -n 1
}

perf_value() {
    awk -F, -v k="$1" '$1 == k { print $2; exit }' <<< "$2"
}

run_perf() {
    local tmp start end stdout
    tmp=$(mktemp)
    start=$(date +%s.%N)
    stdout=$(LC_ALL=C perf stat -x, -e "$PERF_EVENTS" -o "$tmp" "$EXE" 2>/dev/null)
    PERF_RC=$?
    end=$(date +%s.%N)
    PERF_WALL=$(elapsed "$start" "$end")
    PERF_VALS=$(awk -F, '!/^#/ && NF >= 3 { split($3, a, ":"); print a[1] "," $1 }' "$tmp")
    PERF_TIME=$(polybench_time "$stdout")
    rm -f "$tmp"
}

parse_cachegrind() {
    awk '
    {
        sub(/^==[0-9]+==[ \t]*/, "")
        pos = index($0, ":")
        if (pos == 0) next
        key = substr($0, 1, pos - 1)
        rest = substr($0, pos + 1)
        gsub(/[ \t]+/, " ", key)
        sub(/^ /, "", key)
        sub(/ $/, "", key)
        gsub(/,/, "", rest)
        n = 0
        while (match(rest, /[0-9]+/)) {
            n++
            nums[n] = substr(rest, RSTART, RLENGTH)
            rest = substr(rest, RSTART + RLENGTH)
        }
        if (key == "I refs" && n >= 1) { ir = nums[1]; seen["ir"] = 1 }
        if (key == "I1 misses" && n >= 1) { i1mr = nums[1]; seen["i1mr"] = 1 }
        if (key == "LLi misses" && n >= 1) { ilmr = nums[1]; seen["ilmr"] = 1 }
        if (key == "D refs" && n >= 3) { dr = nums[2]; dw = nums[3]; seen["d"] = 1 }
        if (key == "D1 misses" && n >= 3) { d1r = nums[2]; d1w = nums[3]; seen["d1"] = 1 }
        if (key == "LLd misses" && n >= 3) { dlr = nums[2]; dlw = nums[3]; seen["dl"] = 1 }
    }
    END {
        if (seen["ir"] && seen["i1mr"] && seen["ilmr"] && seen["d"] && seen["d1"] && seen["dl"])
            print ir "," i1mr "," ilmr "," dr "," d1r "," dlr "," dw "," d1w "," dlw
        else
            print ",,,,,,,,"
    }'
}

run_valgrind() {
    local out err start end
    out=$(mktemp)
    err=$(mktemp)
    start=$(date +%s.%N)
    valgrind --tool=cachegrind --cache-sim=yes --cachegrind-out-file=/dev/null "$EXE" > "$out" 2> "$err"
    VG_RC=$?
    end=$(date +%s.%N)
    VG_WALL=$(elapsed "$start" "$end")
    VG_STATS=$(parse_cachegrind < "$err")
    VG_TIME=$(polybench_time "$(cat "$out")")
    rm -f "$out" "$err"
}

main() {
    echo "$HEADER" > "$OUT_CSV"

    local ds th bl rep reps
    for ds in "${DATASETS[@]}"; do
        for th in "${NTHREADS[@]}"; do
            for bl in "${BLOCKS[@]}"; do
                echo "[build] ${ds} threads=${th} block=${bl}"
                build "$ds" "$th" "$bl"

                reps=$(perf_reps "$ds")
                for ((rep = 1; rep <= reps; rep++)); do
                    run_perf
                    write_row "$PERF_TIME" perf "$ds" "$th" "$bl" "$rep" "$PERF_WALL" \
                        "$(perf_value task-clock "$PERF_VALS")" \
                        "$(perf_value cycles "$PERF_VALS")" \
                        "$(perf_value instructions "$PERF_VALS")" \
                        "$(perf_value cache-references "$PERF_VALS")" \
                        "$(perf_value cache-misses "$PERF_VALS")" \
                        "$(perf_value branch-misses "$PERF_VALS")" \
                        "" "" "" "" "" "" "" "" "" "$PERF_RC"
                done

                if use_valgrind "$ds" && [ "$bl" -ge "$VALGRIND_MIN_BLOCK" ]; then
                    run_valgrind
                    write_row "$VG_TIME" valgrind "$ds" "$th" "$bl" 1 "$VG_WALL" \
                        "" "" "" "" "" "" "$VG_STATS" "$VG_RC"
                fi
            done
        done
    done
}

main "$@"