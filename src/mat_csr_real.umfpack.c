#include <bf/mat_csr_real.h>

#include <suitesparse/umfpack.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/rand.h>

void bfMatCsrRealLU(BfMat const *mat, BfMat **L, BfMat **U, BfPerm *P, BfPerm *Q) {
  BF_ERROR_BEGIN();

  BfMatCsrReal const *matCsrReal = bfMatConstToMatCsrRealConst(mat);

  BfSize n = bfMatGetNumRows(mat);
  if (n != bfMatGetNumCols(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* Can't exceed range of int since SuperLU uses ints for indices. */
  if (n > INT_MAX)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize nnz = matCsrReal->rowptr[n];

  BfReal *data = matCsrReal->data;

  /** Convert CSR indices from BfSize to int: */

  int *colind = bfMemAlloc(nnz, sizeof(int));
  if (colind == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0; i < nnz; ++i)
    colind[i] = matCsrReal->colind[i];

  int *rowptr = bfMemAlloc(n + 1, sizeof(int));
  if (rowptr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0; i <= n; ++i)
    rowptr[i] = matCsrReal->rowptr[i];

  /* Use UMFPACK to compute the LU decomposition of `mat`. Since
   * UMFPACK operates on CSC sparse matrices, and since `mat` is a CSR
   * sparse matrix, we actually compute the LU decomposition A^T =
   * U^T*L*T. */

  int status;

  void *symbolic;
  status = umfpack_di_symbolic(n, n, rowptr, colind, data, &symbolic, NULL, NULL);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  void *numeric;
  status = umfpack_di_numeric(rowptr, colind, data, symbolic, &numeric, NULL, NULL);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  int Ut_nnz, Lt_nnz, n_row, n_col, Ut_diag_nnz;
  status = umfpack_di_get_lunz(&Ut_nnz, &Lt_nnz, &n_row, &n_col, &Ut_diag_nnz, numeric);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  if ((BfSize)n_row != n || (BfSize)n_col != n)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  int *Ut_rowptr = bfMemAlloc(n + 1, sizeof(int));
  if (Ut_rowptr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  int *Ut_colind = bfMemAlloc(Ut_nnz, sizeof(int));
  if (Ut_colind == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfReal *Ut_data = bfMemAlloc(Ut_nnz, sizeof(BfReal));
  if (Ut_data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  int *Lt_rowptr = bfMemAlloc(n + 1, sizeof(int));
  if (Lt_rowptr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  int *Lt_colind = bfMemAlloc(Lt_nnz, sizeof(int));
  if (Lt_colind == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfReal *Lt_data = bfMemAlloc(Lt_nnz, sizeof(BfReal));
  if (Lt_data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  int *P_data = bfMemAlloc(n, sizeof(int));
  if (P_data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  int *Q_data = bfMemAlloc(n, sizeof(int));
  if (Q_data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfReal *Rs = bfMemAlloc(n, sizeof(BfReal));
  if (Rs == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  int do_recip = 0;
  status = umfpack_di_get_numeric(Ut_rowptr, Ut_colind, Ut_data,
                                  Lt_rowptr, Lt_colind, Lt_data,
                                  P_data, Q_data, NULL, &do_recip, Rs, numeric);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfReal *b = bfMemAlloc(n, sizeof(BfReal));
  bfSeed(0);
  bfRealRandn(n, b);

  BfReal *x = bfMemAlloc(n, sizeof(BfReal));

  // TEST SOLVE
  status = umfpack_di_solve(UMFPACK_At, rowptr, colind, data, x, b, numeric, NULL, NULL);
  if (status)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  FILE *fp;

  fp = fopen("Ut_rowptr.bin", "w");
  fwrite(Ut_rowptr, n + 1, sizeof(int), fp);
  fclose(fp);

  fp = fopen("Ut_colind.bin", "w");
  fwrite(Ut_colind, Ut_nnz, sizeof(int), fp);
  fclose(fp);

  fp = fopen("Ut_data.bin", "w");
  fwrite(Ut_data, Ut_nnz, sizeof(BfReal), fp);
  fclose(fp);

  fp = fopen("Lt_rowptr.bin", "w");
  fwrite(Lt_rowptr, n + 1, sizeof(int), fp);
  fclose(fp);

  fp = fopen("Lt_colind.bin", "w");
  fwrite(Lt_colind, Lt_nnz, sizeof(int), fp);
  fclose(fp);

  fp = fopen("Lt_data.bin", "w");
  fwrite(Lt_data, Lt_nnz, sizeof(BfReal), fp);
  fclose(fp);

  fp = fopen("P.bin", "w");
  fwrite(P_data, n, sizeof(int), fp);
  fclose(fp);

  fp = fopen("Q.bin", "w");
  fwrite(Q_data, n, sizeof(int), fp);
  fclose(fp);

  fp = fopen("Rs.bin", "w");
  fwrite(Rs, n, sizeof(BfReal), fp);
  fclose(fp);

  fp = fopen("x.bin", "w");
  fwrite(x, n, sizeof(BfReal), fp);
  fclose(fp);

  fp = fopen("b.bin", "w");
  fwrite(b, n, sizeof(BfReal), fp);
  fclose(fp);

  (void)L;
  (void)U;
  (void)P;
  (void)Q;
  BF_ASSERT(false);

  BF_ERROR_END() {}

  umfpack_di_free_symbolic(&symbolic);
  umfpack_di_free_numeric(&numeric);
}
