# Analysis of Tiled LU decomposition

This document describe how the tiled version of LU decomposition works, explaining the main advantages and downsides of using this version of the alogrithm.
[`lu_optimized.c`](./lu_optimized.c) contains a C implementation, meanwhile [`lu_tile_results.csv`](./lu_tile_results.csv) contains the profiling results.

## 1. Idea and implementation

### 1.1 Idea

The LU decomposition is an algebric algorithm that is able to decompose a square `NxN` matrix in a product of a lower square matrix `L` and an upper square matrix `U`.

To achieve it, we can directly apply the **Gaussian pivot algorithm** on the original dataset,  calculating the upper triangular matrix by nullifying the elements under the pivot (diagonal of the matrix at every iteration) and consequently obtaining the lower triangular one from the multipliers used in every row operations in the methods.

Despite its simplicity, each pivot must be **calculated sequentially**, underlying an un-deletable sequential dependencies.
Moreover if the matrix has large dimensions, the row update process could not be fitted entirely in the cache, not gaining advantages from dataset’s local and temporal locality.

To increase performance, we have implemented a **parallel tiled decomposition**, which consist of **subdividing the entire original matrix in `tiles`** (also referred as `“blocks”`), each of which can be encapsulated in an `omp task` computed by a **parallel thread**.

### 1.2. Algorithm and implementation

The algorithm is similar to the Gaussian pivot method, where every iteration consists of 3 main steps:

1. **Calculate the Sequential LU (Gaussian method)** :
It consist of using the classical LU gaussian method on tile in position `(0,0)`.
The first block contains the `pivots` of the matrix, which are crucial for the next steps of the process.
This first step is implemented by `sequentialLU` function, that takes the block as an argument and compute the Gaussian pivoting sequentially.
$$A_{kk} = L_{kk}\,U_{kk}$$
Element-wise, for each pivot $p = 1, \dots, B$ and each row $i > p$:

$$l_{ip} = \frac{a_{ip}}{a_{pp}}, \qquad a_{ij} \leftarrow a_{ij} - l_{ip}\,a_{pj} \quad (j > p)$$

2. **Calculate the Triangular Forward substitution** :
Each block on the right or below the pivot, can be seen as a set of **triangular linear systems** that can be directly resolve using `triangular forward substitution method` where the variable vector to resolve is, iteratively, every `columns` (for the right-side tiles) or `row` (for the down-side tiles) of each block.
This second step is implemented by `matrixTriangularUpperSubsitution` and `matrixTriangularLowerSubsitution` which compute the TFS on the right (upper) and lower tiles.
Right-side tiles ($j > k$): solve a lower triangular system for each column of the tile.

$$L_{kk}\,U_{kj} = A_{kj} \quad\Longrightarrow\quad U_{kj} = L_{kk}^{-1}\,A_{kj}$$

Down-side tiles ($i > k$): solve an upper triangular system for each row of the tile.

$$L_{ik}\,U_{kk} = A_{ik} \quad\Longrightarrow\quad L_{ik} = A_{ik}\,U_{kk}^{-1}$$

3. **Calculate the Schur complements on the remaining sub-matrix**:
Like the classical Gaussian method, this step takes the result from the triangular systems (previously executed) and update the remaining `tiles` by subtracting to each original `block` the `GEMM` of the **upper and lower blocks on the same row and columns**.
As a result, it is obtained a **new sub-matrix** on which the entire process can be **recursively called until convergence**.
Update of every remaining tile ($i > k$, $j > k$):

$$A_{ij} \leftarrow A_{ij} - L_{ik}\,U_{kj}$$

It is important to remind that this method **does not remove dependencies**.
In fact, these **3 steps must be computed sequentially** to obtain a correct result and **not falling into `race-conditions`**.
More specifically the `TFS omp tasks` **depends on the result of the `sequential LU`** on the pivot tile, meanwhile the `Schur complement` **is bounded to both right and down side `TFS`**.
The main advantage is that **each block from the same step is completely independent to each other**, parallelizing them using omp parallel tasks.
To implement this feature, each `pragma omp task` has the clause `depends()` where are defined the dependencies of each algorithm steps that must be completed before execute that specific `task` by a `thread`.

## 2. Results Analysis

The advantages grow with the dataset size, until the speedup saturates.
Working on `smaller / medium size tile`, make the **data fit** **into the cache**, optimizing the spatial and temporal locality.
In fact, with higher dataset’s size, an entire row may not entirely fit in the cache during the row-update phase, forcing the process to evict datas even if they will be used in the next steps of the execution.
Meanwhile, using small blocks, the elaborated row is entirely fitted in the cache, possibly re-using it in the next steps of the computation without requesting it to the `DRAM`.

We can evidence the advantages of tiled parallel version by looking at the following reports:

### 2.1 Execution Time & Speedup

Mean execution time in milliseconds. Sequential uses 1 thread, tiled versions use 4 threads.

`Mean execution time in seconds` comparison between **tiled** version and **sequential** one on different datasets and blocksize.

| Dataset | Sequential | 8x8 | 16x16 | 32x32 | 64x64 | 128x128 |
|---|---:|---:|---:|---:|---:|---:|
| MINI | 0.000032 | 0.000182 | 0.000132 | 0.000098 | 0.000097 | 0.000095 |
| SMALL | 0.001782 | 0.003217 | 0.000962 | 0.000898 | 0.002203 | 0.001998 |
| STANDARD | 1.074850 | 1.412324 | 0.457003 | 0.361645 | 0.341256 | 0.421099 |
| LARGE | 7.427436 | 10.755739 | 3.006443 | 1.908307 | 1.706959 | 1.673793 |
| EXTRALARGE | 60.060153 | 81.415253 | 24.178431 | 15.294228 | 13.450494 | 12.922725 |

`Speedup` of the **tiled** version (32x32 tiles, 4 threads) over the **sequential** one.

| Dataset | Sequential (s) | Tiled 32x32 (s) | Speedup |
|---|---:|---:|---:|
| MINI | 0.000032 | 0.000098 | 0.32x |
| SMALL | 0.001782 | 0.000898 | 1.98x |
| STANDARD | 1.074850 | 0.361645 | 2.97x |
| LARGE | 7.427436 | 1.908307 | 3.89x |
| EXTRALARGE | 60.060153 | 15.294228 | 3.93x |

The **parallel tiled** version is `4.5` **times faster** than the *sequential* method on `large datasets`, meanwhile on `small datasets` it reaches only `0.7x` of the sequential speed (i.e. it is **slower**).
In `smaller dataset`, if the tile dimension is too high, the entire matrix could possibly fit in **only one tile**, which is the pivot one computed sequentially, so the execution time is equal to the sequential program with a bias due to the **task management of `openMP`** (even if the other threads are not used).
Even if we chose a tile dimension extremely small to divide the dataset we do not gain any advantages, due to high amount of thread and task **overheads**, making the tiled version **slower**.

Instead, on `large dataset`, the advantages are evident, gaining a **speedup of `4.5x`** on the large datasets.
It is important to underline that, after reaching a speedup of `4.5x` on the `large dataset`, the value **remains the same on the extra large** one, represented in the graph as a `plateau`, due to the high saturation of the `4 threads`.

Considering an original `NxN` matrix, the cost of the entire process sequentially is `(2/3)*N^3`
If we subdivide the matrix in tile of `BxB` we obtain the following costs at each iteration:

- **Sequential LU**:
   
  $$\left(\tfrac{2}{3}B^3\right)\cdot\frac{N}{B} = \tfrac{2}{3}\,N B^2$$ 
  
  where `N/B` is the number of tile.
  
- **TFS**: 
  
  $$B^3 \cdot 2 \cdot \frac{N}{B} = 2\,N B^2$$ 
  
  tiles for every row / column.

- **Schur complement**: 
  
  $$2 \cdot B^3 \cdot \frac{N}{B} = 2\,N B^2$$ 
  
  to complete the entire calculations.

Considering that the `TFS` and `GEMM` can be parallelized, with (only in theory) optimal infinite thread the cost is reduced to the only `sequential LU` for every pivot block, gaining a $$N^2/B^2$$ of **maximum speedup**.
Doing the maths:
$$S_{\max} = \frac{T_1}{T_\infty} = \frac{\tfrac{2}{3}N^3}{\tfrac{2}{3}N B^2} = \frac{N^2}{B^2} = \left(\frac{N}{B}\right)^2$$

### 2.2 Cache Analysis

Cache miss rate = `cache-misses / cache-references` (perf).

On `MINI Dataset`:

| Version | Tile size | Miss rate |
|---|---|---:|
| Sequential | - | 7.68% |
| Tiled | 8x8 | 7.47% |
| Tiled | 16x16 | 7.39% |
| Tiled | 32x32 | 7.59% |
| Tiled | 64x64 | 7.36% |
| Tiled | 128x128 | 7.62% |

On `EXTRALARGE dataset`:

| Version | Tile size | Miss rate |
|---|---|---:|
| Sequential | - | 5.45% |
| Tiled | 8x8 | 2.78% |
| Tiled | 16x16 | 1.20% |
| Tiled | 32x32 | 0.51% |
| Tiled | 64x64 | 3.18% |
| Tiled | 128x128 | 5.07% |

If we look at the **first table**, the behavior on `small dataset` is approximately the same respect the **sequential** program.
This due to the fact that, for `small dataset`, even the sequential program can **fit the matrix in the cache**, so there is not a real advantage in using tiled decomposition.
Moreover, **small blocks just raise the amount of overhead**, not gaining any kind of optimization cache-side.

Meanwhile, for `large dataset`, where the **original matrix (or the row) does not fit entirely in the cache** and the tiled version is much better, having a `miss rate` way more lower than the **sequential** one in `medium size tiles`.

In fact we can see that on `medium blocks (like 32-64 dimension)` the tile is **fine-grained and also maximize the local and temporal reuse of the data**.
For `small blocks`, instead, the **tiling is redundant**, like in the previous example, forcing deconstruction even if its not necessary.
Yet in `higher dimension` tile the advantage is not too high, due to the fact that **big tiles of data does not fit in entirely the cache and does not optimize the data reuse and locality**, making the block (or its row) not fit in the cache and needing a `refill from lower caches` or `DRAM`.

Doing the maths on `L1` cache:

Assuming `double` elements (8 bytes each) and a 32 KB L1 cache the usaege of cache is the following:

- **MINI** `(32x32)`: the whole matrix takes 32 × 32 × 8 = `8 KB`, so it fits entirely in `L1`.
- **SMALL** `(128x128)`: the whole matrix takes 128 × 128 × 8 = `128 KB`, which is **larger** than the 32 KB of L1, so it **does not fit**. However,**the two rows used by the algorithm at each step** (the pivot row and the row being updated) take only `2 × 1 KB`, so **they still fit completely in `L1`**.
In fact, in these examples, the `tiling` does not provide any advantage in cache miss rate.
- **LARGE and EXTRA LARGE** `(2000x2000 and 4000x4000)`: the **original matrix does not fit in the cache**, and **just two rows barely fit**, leaving no room for other data.
For the `large` dataset each row is `2000*8=16k`,really limited.  
For the `extra` large `4000*8 = 32K` only one row fits, and for compute two rows at time, only half of them is loaded.
As a consequence, **half of the pivot row will be evicted and re-loaded in future steps, not optimizing data reuse**.
Meanwhile `tiling`, **makes the entire `medium size block` fit in the cache, optimizing reuses**.
Is critical the `tile size`, in fact with `128x128` size the block does not entirely fit in the cache, making the cache miss-rate increase.
The main gain comes from reusing each tile across consecutive steps, instead of reloading the whole sub-matrix at every step.

## 3. Conclusions

In conclusion we can evince that `tiling` is optimizing the `LU decomposition` gaining an high speedup thanks to `openMP` `task` and `thread` systems.
However the **hyper-parameters** of the algorithm, like `tile size`, are critical to actually gain a substantial optimization and speedup, and not generate to much overhead from the `task manager`.
