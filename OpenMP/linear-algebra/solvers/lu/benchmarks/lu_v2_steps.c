#include <stdio.h>
#include <string.h>

#include <polybench.h>
#include "lu.h"

#if !defined(LU_STEP) || LU_STEP < 1 || LU_STEP > 5
#error "Compile with -DLU_STEP=1, 2, 3, 4 or 5"
#endif

static void init_array(int n, DATA_TYPE POLYBENCH_2D(A, N, N, n, n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      A[i][j] = ((DATA_TYPE)(i + 1) * (j + 1)) / n;
}

static void print_array(int n, DATA_TYPE POLYBENCH_2D(A, N, N, n, n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      fprintf(stderr, DATA_PRINTF_MODIFIER, A[i][j]);
      if ((i * n + j) % 20 == 0)
        fprintf(stderr, "\n");
    }
  fprintf(stderr, "\n");
}

static void kernel_lu(int n, DATA_TYPE POLYBENCH_2D(A, N, N, n, n))
{
  int i, j, k;

  #pragma omp parallel private(i, j, k)
  {
    for (k = 0; k < _PB_N; k++)
    {
#if LU_STEP >= 2
      const DATA_TYPE inv = (DATA_TYPE)1 / A[k][k];
#endif

      #pragma omp for simd
      for (j = k + 1; j < _PB_N; j++)
#if LU_STEP >= 2
        A[k][j] *= inv;
#else
        A[k][j] = A[k][j] / A[k][k];
#endif

      #pragma omp for
      for (i = k + 1; i < _PB_N; i++)
      {
#if LU_STEP >= 3
        const DATA_TYPE aik = A[i][k];
#endif
#if LU_STEP == 4
        DATA_TYPE *Ai = &A[i][0];
        const DATA_TYPE *Ak = &A[k][0];
#elif LU_STEP >= 5
        DATA_TYPE *restrict Ai = &A[i][0];
        const DATA_TYPE *restrict Ak = &A[k][0];
#endif

        #pragma omp simd
        for (j = k + 1; j < _PB_N; j++)
#if LU_STEP >= 4
          Ai[j] -= aik * Ak[j];
#elif LU_STEP >= 3
          A[i][j] -= aik * A[k][j];
#else
          A[i][j] = A[i][j] - A[i][k] * A[k][j];
#endif
      }
    }
  }
}

int main(int argc, char **argv)
{
  int n = N;
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);

  init_array(n, POLYBENCH_ARRAY(A));
  polybench_start_instruments;
  kernel_lu(n, POLYBENCH_ARRAY(A));
  polybench_stop_instruments;
  polybench_print_instruments;
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(A)));
  POLYBENCH_FREE_ARRAY(A);

  return 0;
}
