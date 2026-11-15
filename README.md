# HPC Assignments: Team 6

This GitHub repo contains the assignments of **Team 6** for the **HPC course at Unimore 2026/2027**.

## Repository structure

```text
├── OpenMP
│   ├── dataset_specs.md
│   ├── linear-algebra
│   │   └── solvers
│   │       └── lu
│   │           ├── Makefile
│   │           ├── README.md
│   │           ├── benchmarks
│   │           │   ├── lu_v2_steps.c
│   │           │   ├── perf_valgrind.sh
│   │           │   └── tiled_benchmark.sh
│   │           ├── docs
│   │           │   ├── direct_optimization_analysis.md
│   │           │   └── tiled_version_analysis.md
│   │           ├── lu.c
│   │           ├── lu.h
│   │           ├── lu_optimized.c
│   │           ├── lu_v3.c
│   │           ├── lu_v4.c
│   │           ├── plots_final_analysis
│   │           │   ├── confronto_speedup_lu.svg
│   │           │   └── confronto_tempi_kernel_lu.svg
│   │           ├── presentation
│   │           │   ├── HPC OpenMP Lu Team 6.pdf
│   │           │   └── HPC OpenMP Lu Team 6.pptx
│   │           └── results
│   │               ├── lu_kernel_results.csv
│   │               ├── lu_results.csv
│   │               └── lu_tile_results.csv
│   └── utilities
│       ├── benchmark_list
│       ├── c.mk
│       ├── c2.mk
│       ├── common.mk
│       ├── create_cpped_version.sh
│       ├── papi_counters.list
│       ├── polybench.c
│       ├── polybench.h
│       ├── template-for-new-benchmark.c
│       ├── template-for-new-benchmark.h
│       └── time_benchmark.sh
├── README.md
└── common
    ├── c.mk
    ├── polybench.c
    ├── polybench.h
    └── polybenchUtilFuncts.h
```
where:

- **`OpenMP/`**: assignment on the LU decomposition in OpenMP (source code and analysis)
- **`common/`**: PolyBench, used to measure the execution time and for part of the profiling operations

## Contributors

- Nardone Massimiliano
- Nazzaro Jacopo
- Vecchio Ester
