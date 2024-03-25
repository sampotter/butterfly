#include <bf/chol_csr_real.h>

#include <limits.h>

#include <bf/assert.h>
#include <bf/cholmod.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_csr_real.h>
#include <bf/mem.h>
#include <bf/vec_real.h>

struct BfCholCsrRealImpl {
  cholmod_sparse A;
  cholmod_factor *RFactor;
};

/** Interface: Chol */

static BfCholVtable CHOL_VTABLE = {
  .Solve = (__typeof__(&bfCholSolve))bfCholCsrRealSolve,
  .SolveVec = (__typeof__(&bfCholSolveVec))bfCholCsrRealSolveVec,
};

BfMat *bfCholCsrRealSolve(BfCholCsrReal const *cholCsrReal, BfMat const *B) {
  (void)cholCsrReal;
  (void)B;
  BF_DIE();
}

static BfVec *solve_vecReal(BfCholCsrReal const *cholCsrReal, BfVecReal const *b) {
  BF_ERROR_BEGIN();

  BfVecReal *x = bfVecRealNew();
  HANDLE_ERROR();

  BfSize n = bfVecRealGetSize(b);

  bfVecRealInit(x, n);
  HANDLE_ERROR();

  cholmod_dense b_ = {
    .nrow = n,
    .ncol = 1,
    .nzmax = n,
    .d = n,
    .x = b->data,
    .z = NULL,
    .xtype = CHOLMOD_REAL,
    .dtype = CHOLMOD_DOUBLE,
  };

  cholmod_common *c = bfCholmodGetCommon();

  // TODO: we can probably optimize this a bit by using cholmod_solve2
  // instead of cholmod_solve.
  cholmod_dense *x_ = cholmod_solve(CHOLMOD_A, cholCsrReal->impl->RFactor, &b_, c);

  bfMemCopy(x_->x, n, sizeof(BfReal), x->data);
  cholmod_free_dense(&x_, c);

  BF_DIE();

  BF_ERROR_END() {
    BF_DIE();
  }

  return bfVecRealToVec(x);
}

BfVec *bfCholCsrRealSolveVec(BfCholCsrReal const *cholCsrReal, BfVec const *b) {
  switch (bfVecGetType(b)) {
  case BF_TYPE_VEC_REAL:
    return solve_vecReal(cholCsrReal, bfVecConstToVecRealConst(b));
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

/** Upcasting: CholCsrReal -> Chol */

/** Downcasting: Chol -> CholCsrReal */

/** Implementation: CholCsrReal */

BfCholCsrReal *bfCholCsrRealNew() {
  BF_ERROR_BEGIN();

  BfCholCsrReal *cholCsrReal = bfMemAlloc(1, sizeof(BfCholCsrReal));
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return cholCsrReal;
}

void bfCholCsrRealInit(BfCholCsrReal *cholCsrReal, BfMat const *mat) {
  BF_ERROR_BEGIN();

  bfCholInit(&cholCsrReal->super, &CHOL_VTABLE);

  BfMatCsrReal const *matCsrReal = bfMatConstToMatCsrRealConst(mat);
  HANDLE_ERROR();

  BfSize n = bfMatGetNumRows(mat);
  if (n != bfMatGetNumCols(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* Can't exceed range of int since UMFPACK uses ints for indices. */
  if (n > INT_MAX)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize nnz = matCsrReal->rowptr[n];

  if (!bfCholmodIsStarted())
    bfCholmodStart();

  cholCsrReal->impl->A = (cholmod_sparse) {
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

  cholmod_common *c = bfCholmodGetCommon();

  cholCsrReal->impl->RFactor = cholmod_analyze(&cholCsrReal->impl->A, c);
  if (cholCsrReal->impl->RFactor == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  cholmod_factorize(&cholCsrReal->impl->A, cholCsrReal->impl->RFactor, c);

  if (cholCsrReal->impl->RFactor->minor != n)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BF_ASSERT(cholCsrReal->impl->RFactor->is_ll);
  BF_ASSERT(!cholCsrReal->impl->RFactor->is_super);
  BF_ASSERT(cholCsrReal->impl->RFactor->is_monotonic);

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfCholCsrRealDeinit(BfCholCsrReal *cholCsrReal) {
  BF_DIE(); // TODO: implement
}

void bfCholCsrRealDealloc(BfCholCsrReal **cholCsrReal) {
  BF_DIE(); // TODO: implement
}

void bfCholCsrRealDeinitAndDealloc(BfCholCsrReal **cholCsrReal) {
  BF_DIE(); // TODO: implement
}
