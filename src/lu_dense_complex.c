#include <bf/lu_dense_complex.h>

#include <assert.h>
#include <string.h>

#include <bf/blas.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_dense_complex.h>
#include <bf/mat_perm.h>
#include <bf/mat_product.h>

struct BfLuDenseComplexImpl {
  lapack_int m;
  lapack_int n;
  lapack_complex_double *data;
  lapack_int *ipiv;
};

/** Interface: Lu */

static BfLuVtable LU_VTABLE = {
  .Solve = (__typeof__(&bfLuSolve))bfLuDenseComplexSolve,
  .SolveLower = (__typeof__(&bfLuSolveLower))bfLuDenseComplexSolveLower,
  .SolveUpper = (__typeof__(&bfLuSolveUpper))bfLuDenseComplexSolveUpper,
  .Scale = (__typeof__(&bfLuScale))bfLuDenseComplexScale,
  .SolveVec = (__typeof__(&bfLuSolveVec))bfLuDenseComplexSolveVec,
  .SolveLowerVec = (__typeof__(&bfLuSolveLowerVec))bfLuDenseComplexSolveLowerVec,
  .SolveUpperVec = (__typeof__(&bfLuSolveUpperVec))bfLuDenseComplexSolveUpperVec,
  .ScaleVec = (__typeof__(&bfLuScaleVec))bfLuDenseComplexScaleVec,
  .GetMatView = (__typeof__(&bfLuGetMatView))bfLuDenseComplexGetMatView,
};

BfMat *bfLuDenseComplexSolve(BfLuDenseComplex const *luDenseComplex, BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex const *matDenseComplex = bfMatConstToMatDenseComplexConst(mat);
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(mat);
  if (luDenseComplex->impl->m != m)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize n = bfMatGetNumCols(mat);

  lapack_complex_double const *a = luDenseComplex->impl->data;
  lapack_int lda = luDenseComplex->impl->n;
  lapack_int const *ipiv = luDenseComplex->impl->ipiv;

  BfMatDenseComplex *solution = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(solution, m, n);
  HANDLE_ERROR();

  lapack_complex_double *b = solution->data;

  memcpy(b, matDenseComplex->data, m*n*sizeof(lapack_complex_double));

  lapack_int ldb = n;

  lapack_int info = LAPACKE_zgetrs(LAPACK_ROW_MAJOR, 'N', m, n, a, lda, ipiv, b, ldb);
  if (info != 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  END_ERROR_HANDLING() {
    assert(false);
  }

  return bfMatDenseComplexToMat(solution);
}

BfMat *bfLuDenseComplexSolveLower(BfLuDenseComplex const *luDenseComplex, BfMat const *B, bool permute) {
  (void)luDenseComplex;
  (void)B;
  (void)permute;
  assert(false);
}

BfMat *bfLuDenseComplexSolveUpper(BfLuDenseComplex const *luDenseComplex, BfMat const *B, bool permute) {
  (void)luDenseComplex;
  (void)B;
  (void)permute;
  assert(false);
}

BfMat *bfLuDenseComplexScale(BfLuDenseComplex const *luDenseComplex, BfMat const *B) {
  (void)luDenseComplex;
  (void)B;
  assert(false);
}

BfVec *bfLuDenseComplexSolveVec(BfLuDenseComplex const *luDenseComplex, BfVec const *b) {
  (void)luDenseComplex;
  (void)b;
  assert(false);
}

BfVec *bfLuDenseComplexSolveLowerVec(BfLuDenseComplex const *luDenseComplex, BfVec const *b, bool permute) {
  (void)luDenseComplex;
  (void)b;
  (void)permute;
  assert(false);
}

BfVec *bfLuDenseComplexSolveUpperVec(BfLuDenseComplex const *luDenseComplex, BfVec const *b, bool permute) {
  (void)luDenseComplex;
  (void)b;
  (void)permute;
  assert(false);
}

BfVec *bfLuDenseComplexScaleVec(BfLuDenseComplex const *luDenseComplex, BfVec const *b) {
  (void)luDenseComplex;
  (void)b;
  assert(false);
}

BfMat *bfLuDenseComplexGetMatView(BfLuDenseComplex *luDenseComplex) {
  BEGIN_ERROR_HANDLING();

  BfMatProduct *PLU = bfMatProductNew();
  HANDLE_ERROR();

  bfMatProductInit(PLU);
  HANDLE_ERROR();

  BfMat *P = bfMatPermToMat(
    bfMatPermNewViewFromLapackPivots(
      luDenseComplex->impl->m, luDenseComplex->impl->ipiv));
  HANDLE_ERROR();

  bfMatProductPostMultiply(PLU, P);
  HANDLE_ERROR();

  BfSize n = luDenseComplex->impl->m;
  if (n != luDenseComplex->impl->n)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfComplex *ptr = luDenseComplex->impl->data;

  BfMat *L = bfMatDenseComplexToMat(bfMatDenseComplexNewViewFromPtr(n, n, ptr));
  HANDLE_ERROR();

  L->props |= BF_MAT_PROPS_LOWER_TRI;
  L->props |= BF_MAT_PROPS_UNIT;

  bfMatProductPostMultiply(PLU, L);
  HANDLE_ERROR();

  BfMat *U = bfMatDenseComplexToMat(bfMatDenseComplexNewViewFromPtr(n, n, ptr));
  HANDLE_ERROR();

  U->props |= BF_MAT_PROPS_UPPER_TRI;

  bfMatProductPostMultiply(PLU, U);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    assert(false);
  }

  return bfMatProductToMat(PLU);
}

/** Upcasting: LuDenseComplex -> Lu */

BfLu *bfLuDenseComplexToLu(BfLuDenseComplex *luDenseComplex) {
  return &luDenseComplex->super;
}

/** Downcasting: Lu -> LuDenseComplex */

/** Implementation: LuDenseComplex */

BfLuDenseComplex *bfLuDenseComplexNew() {
  BEGIN_ERROR_HANDLING();

  BfLuDenseComplex *luDenseComplex = malloc(sizeof(BfLuDenseComplex));
  if (luDenseComplex == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return luDenseComplex;
}

void bfLuDenseComplexInit(BfLuDenseComplex *luDenseComplex, BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  bfLuInit(&luDenseComplex->super, &LU_VTABLE);

  BfMatDenseComplex const *matDenseComplex = bfMatConstToMatDenseComplexConst(mat);
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);

  lapack_complex_double *a = malloc(m*n*sizeof(lapack_complex_double));
  if (a == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  memcpy(a, matDenseComplex->data, m*n*sizeof(lapack_complex_double));

  BfSize lda = n;

  lapack_int *ipiv = malloc(m*sizeof(lapack_int));
  if (ipiv == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  lapack_int info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, m, n, a, lda, ipiv);
  if (info)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  luDenseComplex->impl = malloc(sizeof(BfLuDenseComplexImpl));
  if (luDenseComplex->impl == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  luDenseComplex->impl->m = m;
  luDenseComplex->impl->n = n;
  luDenseComplex->impl->data = a;
  luDenseComplex->impl->ipiv = ipiv;

  END_ERROR_HANDLING() {
    assert(false);
  }
}

void bfLuDenseComplexDeinit(BfLuDenseComplex *luDenseComplex) {
  luDenseComplex->impl->m = 0;
  luDenseComplex->impl->n = 0;

  free(luDenseComplex->impl->data);
  luDenseComplex->impl->data = NULL;

  free(luDenseComplex->impl->ipiv);
  luDenseComplex->impl->ipiv = NULL;
}

void bfLuDenseComplexDealloc(BfLuDenseComplex **luDenseComplex) {
  free(*luDenseComplex);
  *luDenseComplex = NULL;
}

void bfLuDenseComplexDeinitAndDealloc(BfLuDenseComplex **luDenseComplex) {
  bfLuDenseComplexDeinit(*luDenseComplex);
  bfLuDenseComplexDealloc(luDenseComplex);
}
