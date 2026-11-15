#!/usr/bin/env bash
# Benchmark the LU sources with perf stat and Valgrind Cachegrind.
# Usage: ./perf_valgrind.sh ["DATASETS=MINI SMALL STANDARD LARGE EXTRALARGE"] \
#        [PERF_REPS=3] ["THREADS=4"] \
#        [OUT=lu_kernel_results.csv] [VALGRIND_REPS=3]
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
UTIL_DIR="$SCRIPT_DIR/../../../utilities"

DATASETS_TEXT="${1:-MINI SMALL STANDARD LARGE EXTRALARGE}"
PERF_REPS="${2:-3}"
THREADS_TEXT="${3:-4}"
OUT="${4:-lu_kernel_results.csv}"
VALGRIND_REPS="${5:-3}"
CC="${CC:-gcc}"
CFLAGS="${CFLAGS:--O2 -g}"
OMPFLAGS="${OMPFLAGS:--fopenmp}"
PERF_EVENTS="${PERF_EVENTS:-task-clock,cycles,instructions,cache-references,cache-misses,branch-misses}"
LU_VERSIONS="${LU_VERSIONS:-v0 v2 v2_simd v2_recip v2_hoist v2_rows v2_restrict v3 v4}"
LU_APPEND="${LU_APPEND:-0}"
CSV_HEADER="tool,version,n,threads,rep,wall_seconds,task_clock_ms,cycles,instructions,cache_references,cache_misses,branch_misses,vg_ir,vg_i1mr,vg_ilmr,vg_dr,vg_d1mr,vg_dlmr,vg_dw,vg_d1mw,vg_dlmw,exit_status,kernel_seconds"

export LC_ALL=C

if [[ ! "$PERF_REPS" =~ ^[1-9][0-9]*$ || ! "$VALGRIND_REPS" =~ ^(0|[1-9][0-9]*)$ ]]; then
  echo "Errore: PERF_REPS deve essere positivo e VALGRIND_REPS non negativo." >&2
  exit 2
fi
if [[ "$LU_APPEND" != 0 && "$LU_APPEND" != 1 ]]; then
  echo "Errore: LU_APPEND deve essere 0 oppure 1." >&2
  exit 2
fi

if ! command -v "$CC" >/dev/null 2>&1; then
  echo "Errore: compilatore '$CC' non trovato." >&2
  exit 2
fi
TOOLS=(perf awk)
if ((VALGRIND_REPS > 0)); then
  TOOLS+=(valgrind)
fi
for tool in "${TOOLS[@]}"; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "Errore: '$tool' non trovato nel PATH." >&2
    exit 2
  fi
done

read -r -a DATASETS <<< "$DATASETS_TEXT"
if (("${#DATASETS[@]}" == 0)); then
  echo "Errore: specificare almeno un dataset o una dimensione." >&2
  exit 2
fi

read -r -a THREADS <<< "$THREADS_TEXT"
if (("${#THREADS[@]}" == 0)); then
  echo "Errore: specificare almeno un numero di thread." >&2
  exit 2
fi
for thread_count in "${THREADS[@]}"; do
  if [[ ! "$thread_count" =~ ^[1-9][0-9]*$ ]]; then
    echo "Errore: numero di thread non valido: '$thread_count'." >&2
    exit 2
  fi
done

if [[ "$OUT" != /* ]]; then
  OUT="$PWD/$OUT"
fi
if [[ "$LU_APPEND" == 0 && -e "$OUT" ]]; then
  echo "Errore: il file di output esiste già: $OUT. Scegliere un altro nome o usare LU_APPEND=1." >&2
  exit 2
fi
mkdir -p -- "$(dirname -- "$OUT")"

TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/lu-perf-valgrind.XXXXXX")"
trap 'rm -rf -- "$TMP_DIR"' EXIT

read -r -a CFLAGS_ARRAY <<< "$CFLAGS"
read -r -a OMPFLAGS_ARRAY <<< "$OMPFLAGS"
SOURCES=(
  "v0:$SCRIPT_DIR/lu.c"
  "v2:$SCRIPT_DIR/lu_v2.c"
  "v2_simd:$SCRIPT_DIR/lu_v2_steps.c:1"
  "v2_recip:$SCRIPT_DIR/lu_v2_steps.c:2"
  "v2_hoist:$SCRIPT_DIR/lu_v2_steps.c:3"
  "v2_rows:$SCRIPT_DIR/lu_v2_steps.c:4"
  "v2_restrict:$SCRIPT_DIR/lu_v2_steps.c:5"
  "v3:$SCRIPT_DIR/lu_v3.c"
  "v4:$SCRIPT_DIR/lu_v4.c"
)
read -r -a REQUESTED_VERSIONS <<< "$LU_VERSIONS"
if (("${#REQUESTED_VERSIONS[@]}" == 0)); then
  echo "Errore: specificare almeno una versione in LU_VERSIONS." >&2
  exit 2
fi
SELECTED_SOURCES=()
for requested in "${REQUESTED_VERSIONS[@]}"; do
  found=0
  for entry in "${SOURCES[@]}"; do
    if [[ "${entry%%:*}" == "$requested" ]]; then
      found=1
      for selected in ${SELECTED_SOURCES[@]+"${SELECTED_SOURCES[@]}"}; do
        if [[ "${selected%%:*}" == "$requested" ]]; then
          echo "Errore: versione ripetuta: '$requested'." >&2
          exit 2
        fi
      done
      SELECTED_SOURCES+=( "$entry" )
      break
    fi
  done
  if ((found == 0)); then
    echo "Errore: versione sconosciuta: '$requested'." >&2
    exit 2
  fi
done
BINS=()
SIZES=()

for dataset in "${DATASETS[@]}"; do
  case "$dataset" in
    STANDARD) size_flag="-DSTANDARD_DATASET=" ;;
    MINI|SMALL|LARGE|EXTRALARGE) size_flag="-D${dataset}_DATASET" ;;
    *)
      if [[ ! "$dataset" =~ ^[1-9][0-9]*$ ]]; then
        echo "Errore: dataset o dimensione non validi: '$dataset'." >&2
        exit 2
      fi
      size_flag="-DN=$dataset"
      ;;
  esac
  N="$(printf '#include "lu.h"\nN\n' | \
    "$CC" "${CFLAGS_ARRAY[@]}" "${OMPFLAGS_ARRAY[@]}" "$size_flag" \
      -I"$SCRIPT_DIR" -E -P -x c -)"
  if [[ ! "$N" =~ ^[1-9][0-9]*$ ]]; then
    echo "Errore: dimensione non valida per '$dataset': '$N'." >&2
    exit 2
  fi
  for previous in ${SIZES[@]+"${SIZES[@]}"}; do
    if [[ "$N" == "$previous" ]]; then
      echo "Errore: dimensione N=$N ripetuta; selezionare dataset distinti." >&2
      exit 2
    fi
  done
  SIZES+=( "$N" )

  for entry in "${SELECTED_SOURCES[@]}"; do
    version="${entry%%:*}"
    source_info="${entry#*:}"
    source_file="${source_info%%:*}"
    step_flag=()
    if [[ "$source_info" == *:* ]]; then
      step_flag=("-DLU_STEP=${source_info##*:}")
    fi
    binary="$TMP_DIR/lu_${version}_n${N}"
    echo "Compilo $version, dataset=$dataset, N=$N..." >&2
    "$CC" "${CFLAGS_ARRAY[@]}" "${OMPFLAGS_ARRAY[@]}" \
      -UPOLYBENCH_GFLOPS -UPOLYBENCH_CYCLE_ACCURATE_TIMER -UPOLYBENCH_PAPI \
      -DPOLYBENCH_TIME -DPOLYBENCH_NO_FLUSH_CACHE "$size_flag" ${step_flag[@]+"${step_flag[@]}"} \
      -I"$UTIL_DIR" "$source_file" "$UTIL_DIR/polybench.c" -lm -o "$binary"
    BINS+=( "$version:$N:$binary" )
  done
done

if [[ "$LU_APPEND" == 1 ]]; then
  if [[ ! -f "$OUT" ]]; then
    echo "Errore: LU_APPEND=1 richiede un CSV esistente: $OUT" >&2
    exit 2
  fi
  plan="$TMP_DIR/planned.csv"
  for item in "${BINS[@]}"; do
    version="${item%%:*}"
    remaining="${item#*:}"
    N="${remaining%%:*}"
    for thread_count in "${THREADS[@]}"; do
      printf '%s,%s,%s\n' "$version" "$N" "$thread_count" >> "$plan"
    done
  done
  if ! awk -F, -v expected="$CSV_HEADER" '
    FNR == NR { planned[$1 "," $2 "," $3] = 1; next }
    FNR == 1 {
      header_seen = 1
      if ($0 != expected) {
        print "Errore: intestazione CSV incompatibile." > "/dev/stderr"
        bad = 1
      }
      next
    }
    {
      if (NF != 23) {
        print "Errore: riga CSV incompleta o non valida: " FNR > "/dev/stderr"
        bad = 1
      }
      if (planned[$2 "," $3 "," $4]) {
        print "Errore: combinazione gia presente: " $2 ", N=" $3 ", thread=" $4 > "/dev/stderr"
        bad = 1
        exit
      }
    }
    END {
      if (!header_seen) {
        print "Errore: CSV vuoto." > "/dev/stderr"
        bad = 1
      }
      exit bad
    }
  ' "$plan" "$OUT"; then
    exit 2
  fi
  last_byte="$(tail -c 1 "$OUT" | od -An -tu1 | tr -d '[:space:]')"
  if [[ "$last_byte" != 10 ]]; then
    printf '\n' >> "$OUT"
  fi
  echo "Aggiungo le nuove misure al CSV esistente: $OUT" >&2
else
  printf '%s\n' "$CSV_HEADER" > "$OUT"
fi

perf_values() {
  local perf_file="$1"
  awk -F, '
    {
      event = $3
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", event)
      value = $1
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
      gsub(/,/, "", value)
      if (value ~ /</ || value == "n/a") {
        print "Avviso: contatore perf non disponibile: " event " (" value ")." > "/dev/stderr"
        value = ""
      }
      if (event == "task-clock") task_clock = value
      else if (event == "cycles") cycles = value
      else if (event == "instructions") instructions = value
      else if (event == "cache-references") cache_references = value
      else if (event == "cache-misses") cache_misses = value
      else if (event == "branch-misses") branch_misses = value
    }
    END {
      printf "%s,%s,%s,%s,%s,%s\n", task_clock, cycles, instructions, \
        cache_references, cache_misses, branch_misses
    }
  ' "$perf_file"
}

cachegrind_values() {
  local cache_file="$1"
  awk '
    $1 == "events:" {
      for (i = 2; i <= NF; i++) event_index[$i] = i - 1
    }
    $1 == "summary:" {
      wanted[1] = "Ir"
      wanted[2] = "I1mr"
      wanted[3] = "ILmr"
      wanted[4] = "Dr"
      wanted[5] = "D1mr"
      wanted[6] = "DLmr"
      wanted[7] = "Dw"
      wanted[8] = "D1mw"
      wanted[9] = "DLmw"
      for (i = 1; i <= 9; i++) {
        event_col = event_index[wanted[i]]
        result[i] = event_col ? $(event_col + 1) : ""
      }
    }
    END {
      for (i = 1; i <= 9; i++) {
        gsub(/,/, "", result[i])
        printf "%s%s", result[i], (i == 9 ? "\n" : ",")
      }
    }
  ' "$cache_file"
}

run_timed() {
  local time_file="$1"
  local kernel_file="$2"
  local log_file="$3"
  shift 3
  local status
  set +e
  { TIMEFORMAT='%R'; time "$@" >"$kernel_file" 2>"$log_file"; } 2> "$time_file"
  status=$?
  set -e
  printf '%s' "$status"
}

kernel_time() {
  awk '
    NF {
      if (NF != 1 || $1 !~ /^[0-9]+([.][0-9]+)?$/) bad = 1
      value = $1
      count++
    }
    END {
      if (bad || count != 1) exit 1
      print value
    }
  ' "$1"
}

echo "Misure perf..."
for item in "${BINS[@]}"; do
  version="${item%%:*}"
  remaining="${item#*:}"
  N="${remaining%%:*}"
  binary="${remaining#*:}"
  for thread_count in "${THREADS[@]}"; do
    for ((rep = 1; rep <= PERF_REPS; rep++)); do
      perf_file="$TMP_DIR/perf_${version}_${N}_${thread_count}_${rep}.csv"
      time_file="$TMP_DIR/time_${version}_${N}_${thread_count}_${rep}.txt"
      kernel_file="$TMP_DIR/kernel_${version}_${N}_${thread_count}_${rep}.txt"
      log_file="$TMP_DIR/perf_${version}_${N}_${thread_count}_${rep}.log"
      echo "perf $version, N=$N, thread=$thread_count, ripetizione=$rep" >&2
      status="$(run_timed "$time_file" "$kernel_file" "$log_file" env OMP_NUM_THREADS="$thread_count" OMP_DYNAMIC=FALSE \
        perf stat -x, -o "$perf_file" -e "$PERF_EVENTS" -- "$binary")"
      if [[ "$status" != 0 ]]; then
        echo "Errore: perf non ha completato $version (exit $status)." >&2
        cat "$log_file" >&2
        [[ ! -f "$perf_file" ]] || cat "$perf_file" >&2
        exit 1
      fi
      if ! kernel_seconds="$(kernel_time "$kernel_file")"; then
        echo "Errore: tempo kernel non valido per $version, N=$N, ripetizione=$rep." >&2
        cat "$kernel_file" "$log_file" >&2
        exit 1
      fi
      wall_seconds="$(awk 'NF { print $1; exit }' "$time_file")"
      if [[ -f "$perf_file" ]]; then
        values="$(perf_values "$perf_file")"
      else
        echo "Errore: perf non ha prodotto i contatori per $version." >&2
        exit 1
      fi
      printf 'perf,%s,%s,%s,%s,%s,%s,,,,,,,,,,%s,%s\n' \
        "$version" "$N" "$thread_count" "$rep" "$wall_seconds" "$values" "$status" "$kernel_seconds" >> "$OUT"
    done
  done
done

echo "Misure Valgrind Cachegrind..."
for item in "${BINS[@]}"; do
  version="${item%%:*}"
  remaining="${item#*:}"
  N="${remaining%%:*}"
  binary="${remaining#*:}"
  for thread_count in "${THREADS[@]}"; do
    for ((rep = 1; rep <= VALGRIND_REPS; rep++)); do
      cache_file="$TMP_DIR/cache_${version}_${N}_${thread_count}_${rep}.out"
      log_file="$TMP_DIR/valgrind_${version}_${N}_${thread_count}_${rep}.log"
      time_file="$TMP_DIR/time_vg_${version}_${N}_${thread_count}_${rep}.txt"
      kernel_file="$TMP_DIR/kernel_vg_${version}_${N}_${thread_count}_${rep}.txt"
      stderr_file="$TMP_DIR/stderr_vg_${version}_${N}_${thread_count}_${rep}.log"
      echo "valgrind $version, N=$N, thread=$thread_count, ripetizione=$rep" >&2
      status="$(run_timed "$time_file" "$kernel_file" "$stderr_file" env OMP_NUM_THREADS="$thread_count" OMP_DYNAMIC=FALSE \
        valgrind --tool=cachegrind --cachegrind-out-file="$cache_file" \
        --log-file="$log_file" "$binary")"
      if [[ "$status" != 0 ]]; then
        echo "Errore: Valgrind non ha completato $version (exit $status)." >&2
        cat "$stderr_file" >&2
        [[ ! -f "$log_file" ]] || cat "$log_file" >&2
        exit 1
      fi
      if ! kernel_seconds="$(kernel_time "$kernel_file")"; then
        echo "Errore: tempo kernel Cachegrind non valido per $version, N=$N, ripetizione=$rep." >&2
        cat "$kernel_file" "$stderr_file" >&2
        exit 1
      fi
      wall_seconds="$(awk 'NF { print $1; exit }' "$time_file")"
      if [[ -f "$cache_file" ]]; then
        values="$(cachegrind_values "$cache_file")"
      else
        echo "Errore: Cachegrind non ha prodotto i contatori per $version." >&2
        [[ ! -f "$log_file" ]] || cat "$log_file" >&2
        exit 1
      fi
      printf 'valgrind,%s,%s,%s,%s,%s,,,,,,,%s,%s,%s\n' \
        "$version" "$N" "$thread_count" "$rep" "$wall_seconds" "$values" "$status" "$kernel_seconds" >> "$OUT"
    done
  done
done

echo "Risultati salvati in: $OUT"


#USED on Nano1:
#OMP_DYNAMIC=FALSE OMP_PROC_BIND=TRUE OMP_PLACES=cores \
#  ./perf_valgrind.sh

# plot_lu_results.py uses wall_seconds, not kernel_seconds.