# Analysis of LU Factorization: Implementation, Optimizations, and Results

## 1. Objective and summary

The program implements LU factorization without pivoting and compares a
serial version with several OpenMP variants. This report describes how the
kernel works, the incremental optimizations, and the results in
[`lu_kernel_results.csv`](./lu_kernel_results.csv). The main measurements
are **kernel** execution times, not total process times.

For `N=4000`, the median of the original serial version is **60.087601 s**,
that of the OpenMP pragma version is **38.538919 s**, that of the SIMD
version is **37.820548 s**, and those of `static` and `static,1` are
**37.536679 s** and **37.481836 s**, respectively. The ratio between the
original time and the `static` time is **1.601×**. The differences between
the last two versions are much smaller than the improvement over the
original.

**Important limitation:** the matrix used by the program is singular. The
measurements describe the behavior of these runs, but do not certify the
numerical correctness of the factorization or its performance on a generic
non-singular matrix.

## 2. Program structure

The program operates on a square `N × N` matrix `A`, normally of type
`double`. Dataset sizes are defined in [`lu.h`](./lu.h).

| Phase | What happens |
|---|---|
| Allocation | Matrix `A` is allocated using the PolyBench macros. |
| Initialization | `init_array` sets `A[i][j] = ((i+1)(j+1))/N`. |
| Timing | The PolyBench macros start the timer immediately before `kernel_lu`. |
| Computation | `kernel_lu` modifies `A` in place, without allocating separate `L` and `U` matrices. |
| Completion | The timer is stopped and printed; printing the matrix to prevent dead-code elimination depends on the PolyBench macros and arguments; memory is then freed. |

### 2.1 What `kernel_lu` does

The serial code is in [`lu.c`](./lu.c). The idea is to process one pivot `k`
at a time, normalize the part of the pivot row to its right, and update the
submatrix below it:

```text
for k = 0, ..., N-1:
    for j = k+1, ..., N-1:
        A[k][j] = A[k][j] / A[k][k]

    for i = k+1, ..., N-1:
        for j = k+1, ..., N-1:
            A[i][j] = A[i][j] - A[i][k] * A[k][j]
```

The diagonal value `A[k][k]` is the pivot. Positions above the diagonal
contain the normalized coefficients of `U`; the diagonal of `U` is
implicitly equal to one. The lower part, including the diagonal, contains
`L`: this is an *in-place* storage convention for `A = L U`. No alternative
pivots are searched for and no rows are swapped.

The loop over `k` is sequential: the next pivot depends on the preceding
updates. For a given `k`, however, the columns of the pivot row are
independent; once normalization is complete, updates to different rows can
also run in parallel. A barrier is needed between normalization and updating,
and a second barrier before the next pivot. The kernel performs **O(N³)**
arithmetic operations and the matrix occupies **O(N²)** memory; the
submatrix shrinks as `k` increases.

## 3. Implemented optimizations

The intermediate versions are implemented in [`lu_v2_steps.c`](./lu_v2_steps.c),
by selecting `LU_STEP` from 1 to 5. Each step retains the previous ones.

| Order | Name in CSV | Change | Goal and limitation |
|---:|---|---|---|
| 0 | `v0` | Original serial version | Reference; the compiler may still apply its own optimizations. |
| 1 | `v2` | One persistent `omp parallel` region and two `omp for` loops per pivot | Distribute columns and rows without recreating the team at every pivot. Implicit barriers preserve dependencies, but incur waiting and synchronization. |
| 2 | `v2_simd` | `omp for simd` for normalization and `omp simd` on the inner update loop | Expose element-level parallelism to vector instructions. The directive alone does not prove that the compiler generated SIMD instructions. |
| 3 | `v2_recip` | One reciprocal `1/A[k][k]` per thread and pivot, followed by multiplications | Reduce the number of divisions; this can change rounding and does not guarantee a speedup if the compiler already performs the transformation. |
| 4 | `v2_hoist` | Load `A[i][k]` outside the inner loop over `j` | Avoid repeatedly loading the same value for each row. |
| 5 | `v2_rows` | Pointers `Ai` and `Ak` to the two rows | Make contiguous accesses explicit and potentially simplify address calculation. |
| 6 | `v2_restrict` | `restrict` on the row pointers | Communicate the aliasing constraint to the compiler and enable optimizations. The pivot and update rows are distinct because `i > k`. |
| 7 | `v3`, **static** | Explicit `static` scheduling, with different index handling | Divide the work into approximately equal blocks; compared with `v2_restrict`, the scheduling clause is not the only change. |
| 8 | `v4`, **static,1** | `schedule(static,1)` on the update rows | Assign rows cyclically; normalization still uses `static`. |

`static` without an explicit chunk distributes approximately equal blocks;
`static,1` assigns one row at a time in round-robin order. For example, with
four threads, the latter assigns rows 0, 4, 8, … to the first thread and
rows 1, 5, 9, … to the second. Block assignment can improve locality, while
cyclic assignment can help with irregular workloads. Here, for a fixed
pivot, the updated rows have loops of the same length, so a load-balancing
benefit is not guaranteed. The sources for the two versions are
[`lu_v3.c`](./lu_v3.c) and [`lu_v4.c`](./lu_v4.c).

## 4. Evaluation methodology

| Aspect | Available measurements or method |
|---|---|
| Dataset | `N=32`, 128, 1024, 2000, 4000. |
| Threads | `threads=4` in every row; `v0` remains serial. The actual number of threads in the team is not recorded. |
| Repetitions | Three for each version and size in the `perf` data. The median and, where useful, the minimum and maximum are reported. |
| Coverage | 135 `perf` rows = 9 versions × 5 sizes × 3 runs. There are also 81 Cachegrind rows, only for `N=32`, 128, and 1024. All reported `exit_status` values are 0. |
| Kernel time | `kernel_seconds` measures the `kernel_lu` call, including OpenMP team creation and synchronization, but excluding initialization, printing, and deallocation. |
| Process time | `wall_seconds` covers the entire invocation and is not used for the speedups reported here. |
| Counters | `perf` counts the entire process, including the OpenMP runtime; Cachegrind simulates events for the entire process and alters execution. Its times are not native execution times. |

The current script [`perf_valgrind.sh`](./perf_valgrind.sh) enables the
PolyBench timer and disables PolyBench's preliminary cache flush. The timer
uses `gettimeofday` and prints six decimal places: it is not a monotonic
clock, and its resolution limits interpretation, especially for small
datasets. The kernel measurements in the `perf` rows are nevertheless
collected **under `perf`**, whose overhead may have an effect. The CSV does
not record the CPU, actual compiler and flags, affinity, frequency, system
load, date, run order, actual number of threads, or multiplexing; these
conditions cannot be reconstructed or verified from the CSV alone.

We define the speedup for each size as:

```text
speedup(version) = median(original kernel_seconds) /
                   median(version kernel_seconds)
```

A value greater than 1 indicates a lower median than the original; it is a
**ratio of medians**, not the median of the per-run speedups. The min–max
values describe only three samples, not confidence intervals. Even when the
ranges do not overlap, this alone does not constitute statistical evidence.

## 5. Times and speedups

### 5.1 Medians for all versions

Times in seconds, `perf`, 4 requested threads, and three repetitions. Values
are reported to the timer's precision:

| Version | N=32 | N=128 | N=1024 | N=2000 | N=4000 |
|---|---:|---:|---:|---:|---:|
| Original `v0` | 0.000031 | 0.001747 | 1.076367 | 7.419219 | 60.087601 |
| OpenMP `v2` | 0.000294 | 0.000939 | 0.597172 | 4.641734 | 38.538919 |
| SIMD | 0.000288 | 0.000884 | 0.582767 | 4.595963 | 37.820548 |
| + reciprocal | 0.000292 | 0.000879 | 0.583056 | 4.574943 | 37.821681 |
| + hoisting | 0.000290 | 0.000838 | 0.568133 | 4.549707 | 37.554381 |
| + pointers | 0.000291 | 0.000830 | 0.591544 | 4.545261 | 37.537228 |
| + `restrict` | 0.000287 | 0.000803 | 0.564595 | 4.544747 | 37.531454 |
| `static` | 0.000286 | 0.000819 | 0.575045 | 4.527864 | 37.536679 |
| `static,1` | 0.000299 | 0.001219 | 0.582838 | 4.522652 | 37.481836 |

For `N=32`, the original takes about 31 microseconds and the parallel
versions about 286–299 microseconds: the workload is too small to amortize
team creation and synchronization. At `N=128`, the original is slower than
the parallel versions shown, but these times are still under a few
milliseconds: small changes and timer precision have a large impact.

### 5.2 Comparison of the main versions

| N | OpenMP speedup | SIMD speedup | static speedup | static,1 speedup |
|---:|---:|---:|---:|---:|
| 32 | 0.105× | 0.108× | 0.108× | 0.104× |
| 128 | 1.860× | 1.976× | 2.133× | 1.433× |
| 1024 | 1.802× | 1.847× | 1.872× | 1.847× |
| 2000 | 1.598× | 1.614× | 1.639× | 1.640× |
| 4000 | 1.559× | 1.589× | 1.601× | 1.603× |

Each value is calculated as **original time / version time**.

![Speedup by size relative to the original](./grafici_relazione_kernel/speedup_per_dataset.svg)

For `N=4000`, moving from the original to OpenMP reduces the median from
60.087601 to 38.538919 s. SIMD is faster than the pragma-only version in the
reported campaign: 37.820548 versus 38.538919 s, about **1.86%** lower.
Compared with SIMD, `static` has a median about **0.75%** lower; this
difference aggregates multiple changes and does not identify the
contribution of any single optimization.

### 5.3 Successive steps and variability at N=4000

![Median times for the successive steps, N=4000](./grafici_relazione_kernel/tempi_passaggi_n4000.svg)

| Version | Median kernel (s) | Min–max (s) | Speedup over original |
|---|---:|---:|---:|
| Original | 60.087601 | 60.007218–60.116156 | 1.000× |
| OpenMP | 38.538919 | 38.535083–38.546236 | 1.559× |
| SIMD | 37.820548 | 37.809507–37.846906 | 1.589× |
| + reciprocal | 37.821681 | 37.806470–38.231522 | 1.589× |
| + hoisting | 37.554381 | 37.544441–37.559699 | 1.600× |
| + pointers | 37.537228 | 37.532494–37.551365 | 1.601× |
| + `restrict` | 37.531454 | 37.524102–37.534687 | 1.601× |
| `static` | 37.536679 | 37.529436–37.556152 | 1.601× |
| `static,1` | 37.481836 | 37.480024–37.490305 | 1.603× |

The reciprocal and SIMD versions have nearly identical medians; one of the
reciprocal runs is slower. The subsequent hoisting step has a lower median,
while pointers and `restrict` produce even smaller differences. These are
comparisons between cumulative versions: the compiler may optimize
automatically, and variability and run order do not make it possible to
attribute each difference with certainty to the single change indicated.

## 6. `static` versus `static,1`

The percentage difference is `(T_static,1 / T_static − 1) × 100`: a negative
value indicates a lower median for `static,1`.

| N | `static` median [min–max] s | `static,1` median [min–max] s | Difference |
|---:|---:|---:|---:|
| 32 | 0.000286 [0.000286–0.000287] | 0.000299 [0.000296–0.000299] | +4.55% |
| 128 | 0.000819 [0.000811–0.000899] | 0.001219 [0.001208–0.001243] | +48.84% |
| 1024 | 0.575045 [0.553745–0.579151] | 0.582838 [0.561424–0.593323] | +1.36% |
| 2000 | 4.527864 [4.524003–4.533600] | 4.522652 [4.505185–4.523440] | −0.12% |
| 4000 | 37.536679 [37.529436–37.556152] | 37.481836 [37.480024–37.490305] | −0.15% |

![Percentage difference between static,1 and static by dataset](./grafici_relazione_kernel/static_vs_static1_dataset.svg)

The largest observed difference is at `N=128`, but it concerns very short
times and does not show that the same trend applies to other inputs. For
`N=1024`, the ranges overlap; for `N=4000`, the three available `static,1`
runs are lower than those for `static`, but the difference between medians
is only 0.054843 s. With three runs and no additional controls for the
machine, run order, affinity, and input, it is not possible to establish
general superiority.

## 7. Instructions, cycles, and IPC

IPC (*instructions per cycle*) is calculated for each run as
`instructions / cycles`; the table and chart report the median of the three
ratios. The instruction and cycle columns are separate medians, expressed in
billions. These counters cover the **entire process**, not just the kernel,
and also include thread and runtime activity.

![IPC for the various versions, N=4000](./grafici_relazione_kernel/ipc_passaggi_n4000.svg)

| Version, N=4000 | Median instructions (billions) | Median cycles (billions) | Median IPC | Median kernel time (s) |
|---|---:|---:|---:|---:|
| Original | 170.956 | 89.358 | 1.913 | 60.087601 |
| OpenMP | 174.394 | 226.936 | 0.769 | 38.538919 |
| SIMD | 110.533 | 222.516 | 0.497 | 37.820548 |
| + reciprocal | 110.801 | 222.614 | 0.498 | 37.821681 |
| + hoisting | 89.650 | 220.970 | 0.406 | 37.554381 |
| + pointers | 89.346 | 220.937 | 0.404 | 37.537228 |
| + `restrict` | 89.264 | 220.886 | 0.404 | 37.531454 |
| `static` | 89.379 | 220.843 | 0.405 | 37.536679 |
| `static,1` | 89.227 | 220.623 | 0.404 | 37.481836 |

The original has a higher IPC but a longer elapsed time: **IPC is not a
speedup**. In parallel runs, cycles are aggregated across threads and also
include runtime and synchronization; multiple threads can make progress
simultaneously. In addition, a single vector instruction can perform
multiple useful operations.

Between SIMD and hoisting, the median instruction count drops from 110.533
to 89.650 billion, about **−18.89%**, while cycles drop from 222.516 to
220.970 billion. IPC decreases from 0.497 to 0.406, even though kernel time
decreases slightly. The data are consistent with fewer instructions over a
similar number of cycles, but do not isolate the cause or prove that
vectorization actually occurred. For `N=4000`, the median process cache
misses are about 4.635 billion for `static` and 4.700 billion for
`static,1`: this is not enough to explain a time difference or, by itself,
determine which scheduling policy has better locality.

## 8. Correctness and final considerations

The initialized matrix is the outer product of the vectors `i+1` and `j+1`,
scaled by `1/N`: for `N>1`, it has rank one and is **singular**. Without
pivoting, the kernel may encounter zero or numerically unstable pivots and
produce `NaN` or `Inf`. An `exit_status` of zero only indicates that the
process terminated without a reported error: **it does not verify that the
factorization is correct**. Rounding differences between division and
reciprocal multiplication may also matter.

In summary:

1. For the larger datasets, the OpenMP versions have lower median times than
   the original; at `N=4000`, the observed speedup of `static` is 1.601×.
2. The micro-optimizations following SIMD show small and non-uniform
   improvements in time, despite an observed reduction in instruction count.
3. `static` and `static,1` are very close for `N=2000` and 4000; `static,1`
   has slightly lower medians in these measurements. The available data do
   not establish a universally better choice.
4. IPC, counters, and variation on small datasets should be interpreted
   together with kernel time and the limitations of the measurement, not in
   isolation.
