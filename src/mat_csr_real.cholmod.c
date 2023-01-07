#include <bf/mat_csr_real.h>

#include <string.h>

#include <bf/cholmod.h>
#include <bf/error.h>
#include <bf/error_macros.h>

BfMat *bfMatCsrRealCholesky(BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfMatCsrReal const *matCsrReal = bfMatConstToMatCsrRealConst(mat);
  HANDLE_ERROR();

  BfSize n = matCsrReal->super.numRows;

  if (n != matCsrReal->super.numCols)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize nnz = matCsrReal->rowptr[n];

  if (!bfCholmodIsStarted())
    bfCholmodStart();

  cholmod_common *c = bfCholmodGetCommon();

  /* NOTE: Cholmod operates on sparse matrices in CSC format, but this
   * is a CSR format matrix. This means we just compute an upper
   * triangular Cholesky factor instead of a lower triangular one. */

  cholmod_sparse A = {
    .nrow = n,
    .ncol = n,
    .nzmax = nnz,
    .p = matCsrReal->rowptr,
    .i = matCsrReal->colind,
    .nz = NULL,
    .x = matCsrReal->data,
    .z = NULL,
    .stype = 0,
    .itype = CHOLMOD_LONG,
    .xtype = CHOLMOD_REAL,
    .dtype = CHOLMOD_DOUBLE,
    .sorted = true,
    .packed = true
  };

  cholmod_factor *R_factor = cholmod_analyze(&A, c);
  if (R_factor == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  cholmod_factorize(&A, R_factor, c);

  if (R_factor->minor != n)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  assert(R_factor->is_ll && !R_factor->is_super && R_factor->is_monotonic);

  cholmod_sparse *R_sparse = cholmod_factor_to_sparse(R_factor, c);

  BfSize *R_rowptr = malloc((n + 1)*sizeof(BfSize));
  memcpy(R_rowptr, R_sparse->p, n*sizeof(BfSize));
  R_rowptr[n] = R_sparse->nzmax;

  BfMatCsrReal *R = bfMatCsrRealNew();
  bfMatCsrRealInit(R, n, n, R_rowptr, R_sparse->i, R_sparse->x);

  END_ERROR_HANDLING()
    bfMatCsrRealDeinitAndDealloc(&R);

  cholmod_free_factor(&R_factor, c);
  cholmod_free_sparse(&R_sparse, c);
  free(R_rowptr);

  return bfMatCsrRealToMat(R);
}
