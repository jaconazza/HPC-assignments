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

/* Array initialization. */
static void init_array(int n, DATA_TYPE POLYBENCH_2D(A, N, N, n, n)) {
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      A[i][j] = ((DATA_TYPE)(i + 1) * (j + 1)) / n;
}

/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static void print_array(int n, DATA_TYPE POLYBENCH_2D(A, N, N, n, n))

{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      fprintf(stderr, DATA_PRINTF_MODIFIER, A[i][j]);
      if ((i * n + j) % 20 == 0)
        fprintf(stderr, "\n");
    }
  fprintf(stderr, "\n");
}

// --------------------------

static void sequentialLU(DATA_TYPE **A, const int startIdx, const int endIdx) {

  for (size_t k = startIdx; k <= endIdx; k++) {
    for (size_t j = k + 1; j <= endIdx; j++)
      A[k][j] = A[k][j] / A[k][k];

    for (size_t i = k + 1; i <= endIdx; i++)
      for (size_t j = k + 1; j <= endIdx; j++)
        A[i][j] = A[i][j] - A[i][k] * A[k][j];
  }
}

// GEMM : A11 = L10*U01 + L11*U11 -> A11 - A10*A01 = B11 next iteration will be
// This GEMM and Substraction is done using 3 loops
static void calculateSchurComplement(DATA_TYPE **A, const int kStart,
                                     const int dimK, const int rowStartIdx,
                                     const int rowDim, const int colStartIdx,
                                     const int colDim) {

  for (size_t i = rowStartIdx; i < rowStartIdx + rowDim; i++)
    for (size_t k = kStart; k < kStart + dimK; k++) {
      DATA_TYPE lik = A[i][k];
      for (size_t j = colStartIdx; j < colStartIdx + colDim; j++)
        A[i][j] -= lik * A[k][j];
    }
}

/*

Compute the matrix resolution of the traingular forward subsitution

L traingular
X = matrix
A = given matrix

the normal triangular substitution is done with vector:

Lx = b

So in the matrix form its simply iterativly done for every column of the matrix
A and X

for i in columns:
  resolve LX[i] = A[i]

calling the normal traingular resolution
*/
static void matrixTriangularUpperSubsitutio(DATA_TYPE **A,
                                            const int columnStart,
                                            const int columnDim,
                                            const int startL, const int dimL) {

  for (size_t k = columnStart; k < columnStart + columnDim; k++) {

    for (size_t row = startL; row < startL + dimL; row++) {
      DATA_TYPE sum = 0;

      for (size_t col = startL; col < row; col++) {
        sum += A[row][col] * A[col][k];
      }

      A[row][k] = (A[row][k] - sum) / A[row][row];
    }
  }
}

// both are traingular forward because:
/*

this function computes:
  for every row i
    x[i] * U = A[i]

  but x[i] and A[i] are orizontal vector, like 1xb
  so the algorithm is a normal forward, but changing the access index:

  (x[i] * U)^T = A[i]^T
  -> U^T * x[i]*T = A[i]^T

  where x[i]^T and A[i]^T are column vector and U^T is a lower traingle matrix
  So -> its a forward triangular subsitution
*/
static void matrixTriangularLowerSubsitution(DATA_TYPE **A, const int rowStart,
                                             const int rowDim, const int startU,
                                             const int dimU) {

  for (size_t row = rowStart; row < rowStart + rowDim; row++) {

    for (size_t col = startU; col < startU + dimU; col++) {
      DATA_TYPE sum = 0;

      for (size_t k = startU; k < col; k++) {
        sum += A[row][k] * A[k][col];
      }

      A[row][col] = A[row][col] - sum;
    }
  }
}

/*

[A00  A01 A02]    [L00  0    0 ]   [U00 U01 U02]
[A10  A11 A12] =  [L10  L11  0 ] · [ 0  U11 U12]
[A20  A21 A22]    [L20  L21  L22]  [ 0   0  U22]

A00 = L00 * U00 -> sequential LU

A01 = L00 * U01 -> resolve linear system Triangular forward subsitution ->
L00 * X = A01

A02 = L00 * U01 -> triangular forward subsitution (iterativly for every
column of A)

A10 = L10 * U00 -> resolve linear system using Triangular Backward
subsitution -> X * U00 = A10 A20 = L20 * U00 -> Triangular backward
substitution for every column of A

A11 = L10·U01 + L11·U11
A12 = L10·U02 + L11·U12
A21 = L20·U01 + L21·U11
A22 = L20·U02 + L21·U12 + L22·U22




or in an alternative form:
the blocks A01 and A02 are fused into a unique rectangular block
the blocks A10 and A20 are fused into a unique rectangular block

so:

A = [A00 A01]
    [A10 A11]

where:

A00 is bxb
A01 is bx(n-b)
A10 is (n-b)xb
A11 is (n-b)x(n-b)
*/
static void kernelBlockLU(const int n, DATA_TYPE POLYBENCH_2D(A, N, N, n, n),
                          const int blockSize) {

  for (size_t k = 0; k < _PB_N; k += blockSize) {
    // compute the LU decomposition on the main base/pivot block
    int kDim = MIN(_PB_N - k, blockSize);
    sequentialLU(A, k, k + kDim - 1);

    for (size_t i = k + blockSize; i < _PB_N; i += blockSize) {
      int iDim = MIN(_PB_N - i, blockSize);
      matrixTriangularUpperSubsitutio(A, i, iDim, k, kDim);
    }

    for (size_t i = k + blockSize; i < _PB_N; i += blockSize) {
      int iDim = MIN(_PB_N - i, blockSize);
      matrixTriangularLowerSubsitution(A, i, iDim, k, kDim);
    }

    // Schur complement calculations
    for (size_t i = k + blockSize; i < _PB_N; i += blockSize) {
      int iDim = MIN(_PB_N - i, blockSize);
      for (size_t j = k + blockSize; j < _PB_N; j += blockSize) {
        int jDim = MIN(_PB_N - j, blockSize);
        // calculate Bij using Schur Complement
        calculateSchurComplement(A, k, kDim, i, iDim, j, jDim);
      }
    }
  }
}

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

  /* Initialize array(s). */
  init_array(n, POLYBENCH_ARRAY(A));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_lu(n, POLYBENCH_ARRAY(A));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(A)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);

  return 0;
}
