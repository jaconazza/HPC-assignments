#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
/* Default data type is double, default size is 1024. */
#include "lu.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))

// /* Array initialization. */
// static void init_array(int n, DATA_TYPE POLYBENCH_2D(A, N, N, n, n)) {
//   int i, j;

//   for (i = 0; i < n; i++)
//     for (j = 0; j < n; j++)
//       A[i][j] = ((DATA_TYPE)(i + 1) * (j + 1)) / n;
// }

/* Array init with det!=1 by getting elements between 0-15 and dividing them
 * by 16.
 * the diagonal is then setted with values from 1 to N, making the amtrix
 * diagonal dominant
 * */
static void init_array(int n, DATA_TYPE POLYBENCH_2D(A, N, N, n, n)) {
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      A[i][j] = (DATA_TYPE)((i + 2 * j) % 16) / 16; // values in [0,1)

  for (i = 0; i < n; i++)
    A[i][i] = (DATA_TYPE)n; // diagonal dominant
}

/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static void print_array(int n, DATA_TYPE POLYBENCH_2D(A, N, N, n, n)) {
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      fprintf(stderr, DATA_PRINTF_MODIFIER, A[i][j]);
      if ((i * n + j) % 20 == 0)
        fprintf(stderr, "\n");
    }
  fprintf(stderr, "\n");
}

// verify the result from the sequential and parrallel decomposition
static int8_t verifyResult(int n, DATA_TYPE POLYBENCH_2D(A, N, N, n, n),
                           DATA_TYPE POLYBENCH_2D(Aref, N, N, n, n)) {
  double maxErr = 0.0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      double e = fabs((double)A[i][j] - (double)Aref[i][j]);
      if (e > maxErr)
        maxErr = e;
    }
  fprintf(stderr, "Max error: %g\n", maxErr);
  return maxErr < 1e-9 ? 1 : -1;
}

// ---- LU Decomposition ----

/**
 * Calculate the LU decomposition using the classical pivoting algorithm.
 * The indexes are used to apply LU decomposition in eventual sub-matrix
 * contained in the given one
 *
 * @param A matrix to decompose
 * @param startIdx starting index
 * @param endIdx ending index
 */
static void sequentialLU(DATA_TYPE POLYBENCH_2D(A, N, N, n, n),
                         const int startIdx, const int endIdx) {
  // iterate on diagonal pivots
  for (size_t k = startIdx; k <= endIdx; k++) {

    // given the pivot, normalize the lower elements to get the multipliers
    for (size_t j = k + 1; j <= endIdx; j++)
      A[k][j] = A[k][j] / A[k][k];

    /**
     * update every row as Ri - mi * Rk
     * where:
     * Ri is the row to update
     * Rk is the pivot row
     * mi is the multiplier for the i-row
     * */
    for (size_t i = k + 1; i <= endIdx; i++)
      for (size_t j = k + 1; j <= endIdx; j++)
        A[i][j] = A[i][j] - A[i][k] * A[k][j];
  }
}

/**
 * Update one block of the trailing submatrix using the Schur complement.
 *
 * Computes Aij -= A(i,k) * A(k,j), where A(i,k) holds the block of L below the
 * pivot block and A(k,j) holds the block of U to the right of the pivot block.
 *
 * @param A matrix, updated in place
 * @param kStart first row/column index of the pivot block
 * @param dimK size of the pivot block k
 * @param rowStartIdx starting row index to update
 * @param rowDim row length
 * @param colStartIdx starting column index to update
 * @param colDim col length
 *
 * @details
 * After splitting the matrix in square blocks of size BxB and
 * resolving the linear traingular systems, (side-blocks U and below-blocks L of
 * the pivot one), the remainging (n-B)x(n-B) block is computer using Schur
 * complement.
 *
 * Example with 9x9 matrix divided in 3x3 blocks A[0-2][0-2]:
 *
 * [A00  A01 A02]    [L00  0    0 ]   [U00 U01 U02]
 * [A10  A11 A12] =  [L10  L11  0 ] · [ 0  U11 U12]
 * [A20  A21 A22]    [L20  L21  L22]  [ 0   0  U22]
 *
 * the entire [A11 A12] sub-matrix has to be computed using Schur complement
 *            [A21 A22]
 *
 * This function is called for every block in the (n-B)x(n-B) sub-matrix and is
 * calculated using GEMM between the L blocks and U blocks on the same row and
 * columns, then subtracting their product to the block examined.
 * The results is a new block (sub-matrix) Aij' which will be decomposed and
 * managed in the next iterations of the main algorithm.
 *
 * Algebric example for 22:
 *
 * A22 = L20*U02 + L21*U12 + L22*U22
 * -> A22 - L20*U02 - L21*U12 = L22*U22 = A22'
 *
 * The decomposition of A22' into L22 and U22 whill be done in next iterations
 */
static void calculateSchurComplement(DATA_TYPE POLYBENCH_2D(A, N, N, n, n),
                                     const int kStart, const int dimK,
                                     const int rowStartIdx, const int rowDim,
                                     const int colStartIdx, const int colDim) {
  // Matrix product between L and U matrix on the same row and column of the
  // analyzed block
  for (size_t i = rowStartIdx; i < rowStartIdx + rowDim; i++)
    for (size_t k = kStart; k < kStart + dimK; k++) {
      DATA_TYPE lik = A[i][k];

      // updaete of every element subtracting the product to Aij
      for (size_t j = colStartIdx; j < colStartIdx + colDim; j++)
        A[i][j] -= lik * A[k][j];
    }
}

/**
 * Compute the triangular forward substitution iterativly on columns of
 * a side-block respect the pivot one.
 *
 * @param A matrix, updated in place
 * @param columnStart starting column index used to repeat the Triangular FS on
 * different columns
 * @param columnDim number of columns
 * @param startL startin index of the Low triangular matrix L (on A)
 * @param dimL dimension of the L matrix
 *
 * @details
 * After resolving the pivoting block, the blocks on the left-side are computed
 * using triangular forward subsitution method.
 * Example:
 *
 * [A00  A01 A02]    [L00  0    0 ]   [U00 U01 U02]
 * [A10  A11 A12] =  [L10  L11  0 ] · [ 0  U11 U12]
 * [A20  A21 A22]    [L20  L21  L22]  [ 0   0  U22]
 *
 * A01 = L00 * U01.
 * L00 is given from the LU decomposition of A00, we can get the value of U01
 * resolving iterativly the traingular linear system on the columns A01 and L00
 *
 * for every column i:
 *    L00 * U01[i] = A01[i]
 *
 * this function is called for every side-block of the pivot one and in-place
 * updating the matrix A
 */
static void
matrixTriangularUpperSubsitutio(DATA_TYPE POLYBENCH_2D(A, N, N, n, n),
                                const int columnStart, const int columnDim,
                                const int startL, const int dimL) {
  // iterate on the columns of A to resolve iterativly the TFS
  for (size_t k = columnStart; k < columnStart + columnDim; k++) {

    // Resolve the TFS
    for (size_t row = startL; row < startL + dimL; row++) {
      DATA_TYPE sum = 0;

      for (size_t col = startL; col < row; col++) {
        sum += A[row][col] * A[col][k];
      }

      A[row][k] = (A[row][k] - sum) / A[row][row];
    }
  }
}

/**
 * Compute the Triangular Forward Substitution iterativly on rows of a
 * lower-block respects the pivot one.
 *
 * @param A matrix, updated in place
 * @param rowstart starting row index used to repeat the Triangular FS on
 * different rows
 * @param rowDim number of rows
 * @param startU startin index of the side triangular matrix U (on A)
 * @param dimU dimension of the U matrix
 *
 * @details
 * After resolving the pivoting block, the lower-blocks are computed
 * using triangular forward subsitution method.
 * It is forward and not backward due to an algebric transpose:
 *
 * [A00  A01 A02]    [L00  0    0 ]   [U00 U01 U02]
 * [A10  A11 A12] =  [L10  L11  0 ] · [ 0  U11 U12]
 * [A20  A21 A22]    [L20  L21  L22]  [ 0   0  U22]
 *
 * A10 = L10 * U00 where U00 is given and upper-triangular.
 * We can then resolve iterativly on every row of A to get the L matrix.
 *
 * L10[i] * U00 = A10[i], even if U is upper-triangular, L10[i] and A10[i] are
 * orizontal vector, we can then transpose both the terms of the equations and
 * getting:
 *
 * (L10[i] * U00)^T = A10[i]^T -> U00^T * L10[i]^T = A10[i]^T
 * which can be resolved using a forward substitution.
 * This is done for every column in A10.
 *
 */
static void
matrixTriangularLowerSubsitution(DATA_TYPE POLYBENCH_2D(A, N, N, n, n),
                                 const int rowStart, const int rowDim,
                                 const int startU, const int dimU) {
  // iterate on the row of A and resolve the TFS
  for (size_t row = rowStart; row < rowStart + rowDim; row++) {

    // Resolve the TFS
    for (size_t col = startU; col < startU + dimU; col++) {
      DATA_TYPE sum = 0;

      for (size_t k = startU; k < col; k++) {
        sum += A[row][k] * A[k][col];
      }

      A[row][col] = A[row][col] - sum;
    }
  }
}

/**
 * Compute the parallel LU decomposition
 *
 * @param n Dimension of the matrix
 * @param A matrix, decomposed in place
 * @param blockSize dimension of the sub-matrix used to split the matrix in
 * tiles (or blocks)
 *
 * @details
 * The entire A matrix is splitted into blocks of dimension blockSizeXblockSize
 * and compute the LU decomposition using algebric tricks.
 * More specifically:
 *
 * - Pivot block (bxb) -> computed using classical LU
 *
 * - Right-side blocks ( bx(n-b) subdivided in square blocks of bxb) -> computed
 * using Triangular Forward substitution
 *
 * - Below-side blocks ( (n-b)xb subdivided in square blocks of bxb) -> computed
 * using Triangular Forward substitution
 *
 * - Rest of the Sub-matrix ( (n-b)x(n-b) subdivided in square blocks of bxb )
 * -> updated using Schur complement and iterativly decomposed by the main
 * iterations of the algortm
 *
 * (if the split it not precise, border blocks are partials of (n-K)xblockSize
 * or blockSizex(n-k) or (n-k)x(n-k) size).
 *
 * Example:
 * Suppose we have an original 9x9 matrix, splitted in 3x3 blocks:
 *
 * [A00  A01 A02]    [L00   0    0 ]   [U00 U01 U02]
 * [A10  A11 A12] =  [L10  L11   0 ] · [ 0  U11 U12]
 * [A20  A21 A22]    [L20  L21  L22]   [ 0   0  U22]
 *
 * The first block in 00 contains the pivots and can be directly decomposed with
 * a classical sequential LU decomposition, obtaining:
 * A00 = L00 * U00
 *
 * The right-side blocks (respect to the pivot one) will contains elements of
 * the U matrix. Following GEMM rules, these elements can be computed as:
 *
 * A01 = L00 * U01
 * A02 = L00 * U01
 * -> knowing L00, resolvable using TFS on every comun of each the block
 *
 *
 * The below-side blocks will contains elements of the L matrix.
 * Following GEMM rules, these elements can be computed as:
 *
 * A10 = L10 * U00
 * A20 = L20 * U00
 * -> knowing U00, resolvable using TFS on every row of each the block
 *
 * The rest of the Submatrix can be updated using the Schur complement and will
 * be decomposed in next iterations of the main algorithm.
 * A11 - L10·U01 = L11·U11 = A11'
 * A12 - L10·U02 = L11·U12 = A12'
 * A21 - L20·U01 = L21·U11 = A21'
 * A22 - L20·U02 - L21·U12 = L22·U22 = A22'
 *
 * After the first iterations:
 *
 * [LU00  U01 U02 ]
 * [L10  A11' A12']
 * [L20  A21' A22']
 *
 * And the algorithm will iterativlt contnue the composition on the submatrix:
 * [A11' A12']
 * [A21' A22']
 */

// NOTA: visto che lavoro a blocchi e un elemento viene toccato solo quando
// tocco il suo blocco, prendo come riferimento dei depends il primo elemento
// del blocco.

// clang-format off
static void kernelBlockLU(const int n, DATA_TYPE POLYBENCH_2D(A, N, N, n, n),
                          const int blockSize) {
  // iterate on every pivot block
  for (size_t k = 0; k < _PB_N; k += blockSize) {

    // compute the LU decomposition on the main base/pivot block
    // manage the dimension if  blockSize > remaingin col/rows
    int kDim = MIN(_PB_N - k, blockSize);

    #pragma omp task shared(A) firstprivate(k, kDim) depend(inout : A[k][k])
    sequentialLU(A, k, k + kDim - 1);

    // compute the TFS on each of the right blocks
    for (size_t i = k + blockSize; i < _PB_N; i += blockSize) {
      #pragma omp task shared(A) firstprivate(i, blockSize, k, kDim) depend(in: A[k][k]) depend(inout: A[k][i])
      {
        int iDim = MIN(_PB_N - i, blockSize);
        matrixTriangularUpperSubsitutio(A, i, iDim, k, kDim);
      }
    }

    // comppute the TFS on each of the below blocks
    for (size_t i = k + blockSize; i < _PB_N; i += blockSize) {
      #pragma omp task shared(A) firstprivate(i, blockSize, k, kDim) depend(in: A[k][k]) depend(inout: A[i][k])
      {
        int iDim = MIN(_PB_N - i, blockSize);
        matrixTriangularLowerSubsitution(A, i, iDim, k, kDim);
      }
    }

    // Schur complement calculations on the rest blocks of the
    // (n-blockSize)x(n-blockSize) submatrix
    for (size_t i = k + blockSize; i < _PB_N; i += blockSize) {
      int iDim = MIN(_PB_N - i, blockSize);

      for (size_t j = k + blockSize; j < _PB_N; j += blockSize) {
        int jDim = MIN(_PB_N - j, blockSize);
        
        // calculate Bij using Schur Complement
        #pragma omp task shared(A) firstprivate(i, j, iDim, jDim) depend(in: A[k][j], A[i][k]) depend(inout: A[i][j])
        calculateSchurComplement(A, k, kDim, i, iDim, j, jDim);
      }
    }
  }
}
// clang-format on

// -----------------------

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static void kernel_lu(int n, DATA_TYPE POLYBENCH_2D(A, N, N, n, n)) {
  int i, j, k;

  for (k = 0; k < _PB_N; k++) {
    for (j = k + 1; j < _PB_N; j++)
      A[k][j] = A[k][j] / A[k][k];
    for (i = k + 1; i < _PB_N; i++)
      for (j = k + 1; j < _PB_N; j++)
        A[i][j] = A[i][j] - A[i][k] * A[k][j];
  }
}

int main(int argc, char **argv) {
  /* Retrieve problem size. */
  int n = N;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(B, DATA_TYPE, N, N, n, n);

  /* Initialize array(s). */
  init_array(n, POLYBENCH_ARRAY(A));
  init_array(n, POLYBENCH_ARRAY(B));

  kernel_lu(n, POLYBENCH_ARRAY(A));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
#pragma omp parallel num_threads(NTHREADS)
#pragma omp single
  kernelBlockLU(n, POLYBENCH_ARRAY(B), BLOCK);

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */

  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(A)));

  int8_t checkEqual = verifyResult(n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));
  fprintf(stderr, "Check Result equality: %d\n", checkEqual);
  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);

  return 0;
}
