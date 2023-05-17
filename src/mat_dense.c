#include <bf/mat_dense.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>

/** Upcasting: */

BfMat *bfMatDenseToMat(BfMatDense *matDense) {
  return &matDense->super;
}

BfMat const *bfMatDenseConstToMatConst(BfMatDense const *matDense) {
  return &matDense->super;
}

/** Downcasting: */

BfMatDense *bfMatToMatDense(BfMat *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_DENSE_REAL)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatDense *)mat;
  }
}

BfMatDense const *bfMatConstToMatDenseConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_DENSE_REAL)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatDense const *)mat;
  }
}

/** Implementation: MatDense */

BfMatDense *bfMatDenseNew() {
  BF_ERROR_BEGIN();

  BfMatDense *mat = bfMemAlloc(1, sizeof(BfMatDense));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return mat;
}

void bfMatDenseInit(BfMatDense *matDense, BfMatVtable *matVtable,
                    BfMatDenseVtable *matDenseVtable,
                    BfSize numRows, BfSize numCols,
                    BfSize rowStride, BfSize colStride) {
  BF_ERROR_BEGIN();

  BfMat *mat = &matDense->super;

  bfMatInit(mat, matVtable, numRows, numCols);
  HANDLE_ERROR();

  matDense->vtable = matDenseVtable;
  matDense->rowStride = rowStride;
  matDense->colStride = colStride;

  END_ERROR_HANDLING()
    bfMatDeinit(mat);
}

void bfMatDenseDeinit(BfMatDense *matDense) {
  bfMatDeinit(&matDense->super);
}

BfSize bfMatDenseGetRowStride(BfMatDense const *matDense) {
  return matDense->rowStride;
}

BfSize bfMatDenseGetColStride(BfMatDense const *matDense) {
  return matDense->colStride;
}
