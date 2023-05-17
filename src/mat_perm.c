#include <bf/mat_perm.h>

#include <bf/assert.h>
#include <bf/blas.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_dense_complex.h>
#include <bf/mem.h>

typedef enum {
  PERM_TYPE_BF,
  PERM_TYPE_LAPACK,
} PermType;

struct BfMatPermImpl {
  PermType permType;
  union {
    BfPerm *perm;
    struct {
      BfSize size;
      BfLapackInt *ipiv;
    } permLapack;
  };
};

static BfSize getSize(BfMatPerm const *matPerm) {
  switch (matPerm->impl->permType) {
  case PERM_TYPE_BF:
    return matPerm->impl->perm->size;
  case PERM_TYPE_LAPACK:
    return matPerm->impl->permLapack.size;
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return BF_SIZE_BAD_VALUE;
  }
}

static BfSize numBytes(BfMatPerm const *matPerm) {
  BfSize numBytes = sizeof(BfMatPermImpl);
  switch (matPerm->impl->permType) {
  case PERM_TYPE_BF:
    numBytes += getSize(matPerm)*sizeof(BfSize);
    break;
  case PERM_TYPE_LAPACK:
    numBytes += getSize(matPerm)*sizeof(BfLapackInt);
    break;
  }
  return numBytes;
}

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .GetType = (__typeof__(&bfMatGetType))bfMatPermGetType,
  .NumBytes = (__typeof__(&bfMatNumBytes))bfMatPermNumBytes,
  .GetNumRows = (__typeof__(&bfMatGetNumRows))bfMatPermGetNumRows,
  .GetNumCols = (__typeof__(&bfMatGetNumCols))bfMatPermGetNumCols,
  .Solve = (__typeof__(&bfMatSolve))bfMatPermSolve,
  .GetInverse = (__typeof__(&bfMatGetInverse))bfMatPermGetInverse
};

BfType bfMatPermGetType(BfMatPerm const *matPerm) {
  (void)matPerm;
  return BF_TYPE_MAT_PERM;
}

BfSize bfMatPermNumBytes(BfMatPerm const *matPerm) {
  return numBytes(matPerm);
}

BfSize bfMatPermGetNumRows(BfMatPerm const *matPerm) {
  return matPerm->super.numRows;
}

BfSize bfMatPermGetNumCols(BfMatPerm const *matPerm) {
  return matPerm->super.numCols;
}

static BfMat *solve_matDenseComplex(BfMatPerm const *matPerm, BfMat const *otherMat) {
  BF_ERROR_BEGIN();

  BfSize m = bfMatGetNumRows(otherMat);
  BfSize n = bfMatGetNumCols(otherMat);

  if (getSize(matPerm) != m)
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfMatDenseComplex const *otherMatDenseComplex = bfMatConstToMatDenseComplexConst(otherMat);

  BfMatDenseComplex *result = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(result, m, n);
  HANDLE_ERROR();

  if (matPerm->impl->permType == PERM_TYPE_BF) {
    for (BfSize i = 0; i < m; ++i) {
      BfComplex *writePtr = result->data + matPerm->impl->perm->index[i]*result->rowStride;
      BfComplex const *readPtr = otherMatDenseComplex->data + i*otherMatDenseComplex->rowStride;
      bfMemCopy(readPtr, n, sizeof(BfComplex), writePtr);
    }
  }

  else if (matPerm->impl->permType == PERM_TYPE_LAPACK) {
    /* Copy over data first: */
    for (BfSize i = 0; i < m; ++i) {
      BfComplex *writePtr = result->data + i*result->rowStride;
      BfComplex const *readPtr = otherMatDenseComplex->data + i*otherMatDenseComplex->rowStride;
      bfMemCopy(readPtr, n, sizeof(BfComplex), writePtr);
    }

    /* Apply row pivots using LAPACKE: */
    BfSize k1 = 1;
    BfSize k2 = matPerm->impl->permLapack.size;
    BfLapackInt const *ipiv = matPerm->impl->permLapack.ipiv;
    LAPACKE_zlaswp(LAPACK_ROW_MAJOR, n, result->data, n, k1, k2, ipiv, 1);
  }

  else RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  END_ERROR_HANDLING() {
    BF_ASSERT(false);
  }

  return bfMatDenseComplexToMat(result);
}

BfMat *bfMatPermSolve(BfMatPerm const *matPerm, BfMat const *mat) {\
  switch (bfMatGetType(mat)) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    return solve_matDenseComplex(matPerm, mat);
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    BF_ASSERT(false);
  }
}

BfMat *bfMatPermGetInverse(BfMatPerm const *matPerm) {
  BF_ERROR_BEGIN();

  if (matPerm->impl->permType != PERM_TYPE_BF)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfPerm permInverse = bfPermGetReversePerm(matPerm->impl->perm);
  HANDLE_ERROR();

  BfMatPerm *matPermInverse = bfMatPermNewFromPerm(&permInverse);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    BF_ASSERT(false);
  }

  bfPermDeinit(&permInverse);

  return bfMatPermToMat(matPermInverse);
}

/** Upcasting: MatPerm -> Mat */

BfMat *bfMatPermToMat(BfMatPerm *matPerm) {
  return &matPerm->super;
}

/** Downcasting: Mat -> MatPerm */

/** Implementation: MatPerm */

BfMatPerm *bfMatPermNew() {
  BF_ERROR_BEGIN();

  BfMatPerm *matPerm = bfMemAlloc(1, sizeof(BfMatPerm));
  if (matPerm == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return matPerm;
}

BfMatPerm *bfMatPermNewFromPerm(BfPerm const *perm) {
  BF_ERROR_BEGIN();

  BfMatPerm *matPerm = bfMatPermNew();
  HANDLE_ERROR();

  bfMatPermInitFromPerm(matPerm, perm);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    BF_ASSERT(false);
  }

  return matPerm;
}

BfMatPerm *bfMatPermNewViewFromLapackPivots(BfSize size, BfLapackInt *ipiv) {
  BF_ERROR_BEGIN();

  BfMatPerm *matPerm = bfMatPermNew();
  HANDLE_ERROR();

  bfMatPermInitViewFromLapackPivots(matPerm, size, ipiv);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    BF_ASSERT(false);
  }

  return matPerm;
}

void bfMatPermInitFromPerm(BfMatPerm *matPerm, BfPerm const *perm) {
  BF_ERROR_BEGIN();

  bfMatInit(&matPerm->super, &MAT_VTABLE, perm->size, perm->size);
  HANDLE_ERROR();

  matPerm->impl = bfMemAlloc(perm->size, sizeof(BfMatPermImpl));
  matPerm->impl->permType = PERM_TYPE_BF;
  matPerm->impl->perm = bfPermCopy(perm);

  END_ERROR_HANDLING() {
    BF_ASSERT(false);
  }
}

void bfMatPermInitViewFromLapackPivots(BfMatPerm *matPerm, BfSize size, BfLapackInt *ipiv) {
  BF_ERROR_BEGIN();

  bfMatInit(&matPerm->super, &MAT_VTABLE, size, size);
  HANDLE_ERROR();

  matPerm->super.props |= BF_MAT_PROPS_VIEW;

  matPerm->impl = bfMemAlloc(1, sizeof(BfMatPermImpl));
  matPerm->impl->permType = PERM_TYPE_LAPACK;
  matPerm->impl->permLapack.size = size;
  matPerm->impl->permLapack.ipiv = ipiv;

  END_ERROR_HANDLING() {
    BF_ASSERT(false);
  }
}

void bfMatPermDeinit(BfMatPerm *matPerm) {
  BF_ASSERT(matPerm->impl->permType == PERM_TYPE_LAPACK);

  if (!(matPerm->super.props &= BF_MAT_PROPS_VIEW))
    free(matPerm->impl->permLapack.ipiv);

  free(matPerm->impl);
  matPerm->impl = NULL;
}

void bfMatPermDealloc(BfMatPerm **matPerm) {
  free(*matPerm);
  *matPerm = NULL;
}

void bfMatPermDelete(BfMatPerm **matPerm) {
  bfMatPermDeinit(*matPerm);
  bfMatPermDealloc(matPerm);
}
