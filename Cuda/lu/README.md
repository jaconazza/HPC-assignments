# OpenMP LU decomposition in CUDA

This folder contains the implementation, analysis and results of the optimization of the **LU decomposition using `OpenMP`**.

## Contents

- **Source code** of the different optimized versions (direct and tiled)
- **Results** in the [`results`](results/) folder
- **Plots** of the speedup comparison in the [`plot_final_analysis`](plot_final_analysis/) folder
- **Documentation** for both the direct and the tiled optimization in [`docs`](docs/) folder
- **Presentation** of the work (`.pdf` and `.pptx`) in [`presentation`](presentation/) folder
- **Benchmarks** scripts for fast profiling [`benchmarks`](benchmarks/) folder

## Build

To build the project, use the following command:

```bash
make EXERCISE=<program_name>.c EXT_CFLAGS="-D<dataset>_DATASET -DPOLYBENCH_TIME -g -Wall" clean all run
```
