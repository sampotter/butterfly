#include <bf/lu_csr_real.h>

#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_csr_real.h>
#include <bf/vec_real.h>

#include <suitesparse/umfpack.h>

struct BfLuCsrReal {
  int *rowptr;
  int *colind;
  BfReal *data;

  void *symbolic;
  void *numeric;
};

BfLuCsrReal *bfLuCsrRealNew() {
  BEGIN_ERROR_HANDLING();

  BfLuCsrReal *luCsrReal = malloc(sizeof(BfLuCsrReal));
  if (luCsrReal == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return luCsrReal;
}

void bfLuCsrRealInit(BfLuCsrReal *luCsrReal, BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfMatCsrReal const *matCsrReal = bfMatConstToMatCsrRealConst(mat);
  HANDLE_ERROR();

  BfSize n = bfMatGetNumRows(mat);
  if (n != bfMatGetNumCols(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* Can't exceed range of int since UMFPACK uses ints for indices. */
  if (n > INT_MAX)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize nnz = matCsrReal->rowptr[n];

  luCsrReal->rowptr = malloc((n + 1)*sizeof(int));
  if (luCsrReal->rowptr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0; i <= n; ++i)
    luCsrReal->rowptr[i] = matCsrReal->rowptr[i];

  luCsrReal->colind = malloc(nnz*sizeof(int));
  if (luCsrReal->colind == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0; i < nnz; ++i)
    luCsrReal->colind[i] = matCsrReal->colind[i];

  luCsrReal->data = malloc(nnz*sizeof(BfReal));
  if (luCsrReal->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  memcpy(luCsrReal->data, matCsrReal->data, nnz*sizeof(BfReal));

  /* Use UMFPACK to compute the LU decomposition of `mat`. Since
   * UMFPACK operates on CSC sparse matrices, and since `mat` is a CSR
   * sparse matrix, we actually compute the LU decomposition A^T =
   * U^T*L*T. */

  int status;

  status = umfpack_di_symbolic(
    n, n, luCsrReal->rowptr, luCsrReal->colind, luCsrReal->data,
    &luCsrReal->symbolic, NULL, NULL);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  status = umfpack_di_numeric(
    luCsrReal->rowptr, luCsrReal->colind, luCsrReal->data,
    luCsrReal->symbolic, &luCsrReal->numeric, NULL, NULL);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  END_ERROR_HANDLING() {
    bfLuCsrRealDeinit(luCsrReal);
  }
}

void bfLuCsrRealDeinit(BfLuCsrReal *luCsrReal) {
  free(luCsrReal->rowptr);
  free(luCsrReal->colind);
  free(luCsrReal->data);

  luCsrReal->rowptr = NULL;
  luCsrReal->colind = NULL;
  luCsrReal->data = NULL;

  umfpack_di_free_numeric(&luCsrReal->numeric);
  luCsrReal->numeric = NULL;

  umfpack_di_free_symbolic(&luCsrReal->symbolic);
  luCsrReal->symbolic = NULL;
}

void bfLuCsrRealDealloc(BfLuCsrReal **luCsrReal) {
  free(*luCsrReal);
  *luCsrReal = NULL;
}

void bfLuCsrRealDeinitAndDealloc(BfLuCsrReal **luCsrReal) {
  bfLuCsrRealDeinit(*luCsrReal);
  bfLuCsrRealDealloc(luCsrReal);
}

static BfVec *solve_vecReal(BfLuCsrReal const *luCsrReal, BfVecReal const *b) {
  BEGIN_ERROR_HANDLING();

  BfVecReal *x = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(x, b->super.size);
  HANDLE_ERROR();

  int status = umfpack_di_solve(
    UMFPACK_At,
    luCsrReal->rowptr,
    luCsrReal->colind,
    luCsrReal->data,
    x->data,
    b->data,
    luCsrReal->numeric,
    NULL,
    NULL);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  END_ERROR_HANDLING() {
    bfVecRealDeinitAndDealloc(&x);
  }

  return bfVecRealToVec(x);
}

BfVec *bfLuCsrRealSolveVec(BfLuCsrReal const *luCsrReal, BfVec const *b) {
  switch (bfVecGetType(b)) {
  case BF_TYPE_VEC_REAL:
    return solve_vecReal(luCsrReal, bfVecConstToVecRealConst(b));
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

static BfVec *solveLowerVec_vecReal(BfLuCsrReal const *luCsrReal, BfVecReal const *b, bool permute) {
  BEGIN_ERROR_HANDLING();

  BfVecReal *x = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(x, b->super.size);
  HANDLE_ERROR();

  int status = umfpack_di_solve(
    permute ? UMFPACK_Q_Ut : UMFPACK_Ut,
    luCsrReal->rowptr,
    luCsrReal->colind,
    luCsrReal->data,
    x->data,
    b->data,
    luCsrReal->numeric,
    NULL,
    NULL);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  END_ERROR_HANDLING() {
    bfVecRealDeinitAndDealloc(&x);
  }

  return bfVecRealToVec(x);
}

BfVec *bfLuCsrRealSolveLowerVec(BfLuCsrReal const *luCsrReal, BfVec const *b, bool permute) {
  switch (bfVecGetType(b)) {
  case BF_TYPE_VEC_REAL:
    return solveLowerVec_vecReal(luCsrReal, bfVecConstToVecRealConst(b), permute);
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

static BfVec *solveUpperVec_vecReal(BfLuCsrReal const *luCsrReal, BfVecReal const *b, bool permute) {
  BEGIN_ERROR_HANDLING();

  BfVecReal *x = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(x, b->super.size);
  HANDLE_ERROR();

  int status = umfpack_di_solve(
    permute ? UMFPACK_Lt_P : UMFPACK_Lt,
    luCsrReal->rowptr,
    luCsrReal->colind,
    luCsrReal->data,
    x->data,
    b->data,
    luCsrReal->numeric,
    NULL,
    NULL);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  END_ERROR_HANDLING() {
    bfVecRealDeinitAndDealloc(&x);
  }

  return bfVecRealToVec(x);
}

BfVec *bfLuCsrRealSolveUpperVec(BfLuCsrReal const *luCsrReal, BfVec const *b, bool permute) {
  switch (bfVecGetType(b)) {
  case BF_TYPE_VEC_REAL:
    return solveUpperVec_vecReal(luCsrReal, bfVecConstToVecRealConst(b), permute);
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

static BfVec *scale_vecReal(BfLuCsrReal const *luCsrReal, BfVecReal const *b) {
  BEGIN_ERROR_HANDLING();

  BfVecReal *x = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(x, b->super.size);
  HANDLE_ERROR();

  int status = umfpack_di_scale(x->data, b->data, luCsrReal->numeric);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  END_ERROR_HANDLING() {
    bfVecRealDeinitAndDealloc(&x);
  }

  return bfVecRealToVec(x);
}

BfVec *bfLuCsrRealScaleVec(BfLuCsrReal const *luCsrReal, BfVec const *b) {
  switch (bfVecGetType(b)) {
  case BF_TYPE_VEC_REAL:
    return scale_vecReal(luCsrReal, bfVecConstToVecRealConst(b));
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

void bfLuCsrRealDump(BfLuCsrReal const *luCsrReal) {
  BEGIN_ERROR_HANDLING();

  int status = 0;

  int lnz, unz, n_row, n_col, nz_udiag;
  status = umfpack_di_get_lunz(&lnz, &unz, &n_row, &n_col, &nz_udiag, luCsrReal->numeric);
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
    Lp, Lj, Lx, Up, Ui, Ux, P, Q, Dx, &do_recip, Rs, luCsrReal->numeric);
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
