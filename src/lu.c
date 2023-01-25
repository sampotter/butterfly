#include <bf/lu.h>

#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_csr_real.h>
#include <bf/vec_real.h>

#include <suitesparse/umfpack.h>

struct BfLu {
  int *rowptr;
  int *colind;
  BfReal *data;

  void *symbolic;
  void *numeric;
};

BfLu *bfLuNew() {
  BEGIN_ERROR_HANDLING();

  BfLu *lu = malloc(sizeof(BfLu));
  if (lu == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return lu;
}

void bfLuInit(BfLu *lu, BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfMatCsrReal const *matCsrReal = bfMatConstToMatCsrRealConst(mat);
  HANDLE_ERROR();

  BfSize n = bfMatGetNumRows(mat);
  if (n != bfMatGetNumCols(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* Can't exceed range of int since SuperLU uses ints for indices. */
  if (n > INT_MAX)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize nnz = matCsrReal->rowptr[n];

  lu->rowptr = malloc((n + 1)*sizeof(int));
  if (lu->rowptr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0; i <= n; ++i)
    lu->rowptr[i] = matCsrReal->rowptr[i];

  lu->colind = malloc(nnz*sizeof(int));
  if (lu->colind == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0; i < nnz; ++i)
    lu->colind[i] = matCsrReal->colind[i];

  lu->data = malloc(nnz*sizeof(BfReal));
  if (lu->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  memcpy(lu->data, matCsrReal->data, nnz*sizeof(BfReal));

  /* Use UMFPACK to compute the LU decomposition of `mat`. Since
   * UMFPACK operates on CSC sparse matrices, and since `mat` is a CSR
   * sparse matrix, we actually compute the LU decomposition A^T =
   * U^T*L*T. */

  int status;

  status = umfpack_di_symbolic(
    n, n, lu->rowptr, lu->colind, lu->data, &lu->symbolic, NULL, NULL);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  status = umfpack_di_numeric(
    lu->rowptr, lu->colind, lu->data, lu->symbolic, &lu->numeric, NULL, NULL);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  END_ERROR_HANDLING() {
    bfLuDeinit(lu);
  }
}

void bfLuDeinit(BfLu *lu) {
  free(lu->rowptr);
  free(lu->colind);
  free(lu->data);

  lu->rowptr = NULL;
  lu->colind = NULL;
  lu->data = NULL;

  umfpack_di_free_numeric(&lu->numeric);
  lu->numeric = NULL;

  umfpack_di_free_symbolic(&lu->symbolic);
  lu->symbolic = NULL;
}

void bfLuDealloc(BfLu **lu) {
  free(*lu);
  *lu = NULL;
}

void bfLuDeinitAndDealloc(BfLu **lu) {
  bfLuDeinit(*lu);
  bfLuDealloc(lu);
}

static BfVec *solve_vecReal(BfLu const *lu, BfVecReal const *b) {
  BEGIN_ERROR_HANDLING();

  BfVecReal *x = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(x, b->super.size);
  HANDLE_ERROR();

  int status = umfpack_di_solve(
    UMFPACK_At,
    lu->rowptr,
    lu->colind,
    lu->data,
    x->data,
    b->data,
    lu->numeric,
    NULL,
    NULL);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  END_ERROR_HANDLING() {
    bfVecRealDeinitAndDealloc(&x);
  }

  return bfVecRealToVec(x);
}

BfVec *bfLuSolve(BfLu const *lu, BfVec const *b) {
  switch (bfVecGetType(b)) {
  case BF_TYPE_VEC_REAL:
    return solve_vecReal(lu, bfVecConstToVecRealConst(b));
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

static BfVec *solveLowerVec_vecReal(BfLu const *lu, BfVecReal const *b, bool permute) {
  BEGIN_ERROR_HANDLING();

  BfVecReal *x = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(x, b->super.size);
  HANDLE_ERROR();

  int status = umfpack_di_solve(
    permute ? UMFPACK_Q_Ut : UMFPACK_Ut,
    lu->rowptr,
    lu->colind,
    lu->data,
    x->data,
    b->data,
    lu->numeric,
    NULL,
    NULL);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  END_ERROR_HANDLING() {
    bfVecRealDeinitAndDealloc(&x);
  }

  return bfVecRealToVec(x);
}

BfVec *bfLuSolveLowerVec(BfLu const *lu, BfVec const *b, bool permute) {
  switch (bfVecGetType(b)) {
  case BF_TYPE_VEC_REAL:
    return solveLowerVec_vecReal(lu, bfVecConstToVecRealConst(b), permute);
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

static BfVec *solveUpperVec_vecReal(BfLu const *lu, BfVecReal const *b, bool permute) {
  BEGIN_ERROR_HANDLING();

  BfVecReal *x = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(x, b->super.size);
  HANDLE_ERROR();

  int status = umfpack_di_solve(
    permute ? UMFPACK_Lt_P : UMFPACK_Lt,
    lu->rowptr,
    lu->colind,
    lu->data,
    x->data,
    b->data,
    lu->numeric,
    NULL,
    NULL);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  END_ERROR_HANDLING() {
    bfVecRealDeinitAndDealloc(&x);
  }

  return bfVecRealToVec(x);
}

BfVec *bfLuSolveUpperVec(BfLu const *lu, BfVec const *b, bool permute) {
  switch (bfVecGetType(b)) {
  case BF_TYPE_VEC_REAL:
    return solveUpperVec_vecReal(lu, bfVecConstToVecRealConst(b), permute);
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

static BfVec *scale_vecReal(BfLu const *lu, BfVecReal const *b) {
  BEGIN_ERROR_HANDLING();

  BfVecReal *x = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(x, b->super.size);
  HANDLE_ERROR();

  int status = umfpack_di_scale(x->data, b->data, lu->numeric);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  END_ERROR_HANDLING() {
    bfVecRealDeinitAndDealloc(&x);
  }

  return bfVecRealToVec(x);
}

BfVec *bfLuScale(BfLu const *lu, BfVec const *b) {
  switch (bfVecGetType(b)) {
  case BF_TYPE_VEC_REAL:
    return scale_vecReal(lu, bfVecConstToVecRealConst(b));
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

void bfLuDump(BfLu const *lu) {
  BEGIN_ERROR_HANDLING();

  int status = 0;

  int lnz, unz, n_row, n_col, nz_udiag;
  status = umfpack_di_get_lunz(&lnz, &unz, &n_row, &n_col, &nz_udiag, lu->numeric);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  assert(n_row == n_col);
  int n = n_row;

  int *Lp = malloc((n + 1)*sizeof(int));
  if (Lp == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  int *Lj = malloc(lnz*sizeof(int));
  if (Lj == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  double *Lx = malloc(lnz*sizeof(double));
  if (Lx == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  int *Up = malloc((n + 1)*sizeof(int));
  if (Up == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  int *Ui = malloc(unz*sizeof(int));
  if (Ui == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  double *Ux = malloc(unz*sizeof(double));
  if (Ux == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  int *P = malloc(n*sizeof(int));
  if (P == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  int *Q = malloc(n*sizeof(int));
  if (Q == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  double *Dx = malloc(n*sizeof(double));
  if (Dx == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  double *Rs = malloc(n*sizeof(double));
  if (Rs == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  int do_recip;
  status = umfpack_di_get_numeric(
    Lp, Lj, Lx, Up, Ui, Ux, P, Q, Dx, &do_recip, Rs, lu->numeric);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  FILE *fp = NULL;

  fp = fopen("Lp.bin", "w");
  fwrite(Lp, sizeof(int), n + 1, fp);
  fclose(fp);

  fp = fopen("Lj.bin", "w");
  fwrite(Lj, sizeof(int), lnz, fp);
  fclose(fp);

  fp = fopen("Lx.bin", "w");
  fwrite(Lx, sizeof(double), lnz, fp);
  fclose(fp);

  fp = fopen("Up.bin", "w");
  fwrite(Up, sizeof(int), n + 1, fp);
  fclose(fp);

  fp = fopen("Ui.bin", "w");
  fwrite(Ui, sizeof(int), unz, fp);
  fclose(fp);

  fp = fopen("Ux.bin", "w");
  fwrite(Ux, sizeof(double), unz, fp);
  fclose(fp);

  fp = fopen("P.bin", "w");
  fwrite(P, sizeof(int), n, fp);
  fclose(fp);

  fp = fopen("Q.bin", "w");
  fwrite(Q, sizeof(int), n, fp);
  fclose(fp);

  fp = fopen("Dx.bin", "w");
  fwrite(Dx, sizeof(double), n, fp);
  fclose(fp);

  fp = fopen("Rs.bin", "w");
  fwrite(Rs, sizeof(double), n, fp);
  fclose(fp);

  END_ERROR_HANDLING() {}

  free(Lp);
  free(Lj);
  free(Lx);
  free(Up);
  free(Ui);
  free(Ux);
  free(P);
  free(Q);
  free(Dx);
  free(Rs);
}
