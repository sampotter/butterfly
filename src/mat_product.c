#include <bf/mat_product.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>

/** Interface: MatProduct */

static BfMatVtable MAT_VTABLE = {
  .Copy = (__typeof__(&bfMatProductCopy))bfMatProductCopy,
  .Delete = (__typeof__(&bfMatProductDelete))bfMatProductDelete,
  .GetType = (__typeof__(&bfMatProductGetType))bfMatProductGetType,
  .GetNumRows = (__typeof__(&bfMatProductGetNumRows))bfMatProductGetNumRows,
  .GetNumCols = (__typeof__(&bfMatProductGetNumCols))bfMatProductGetNumCols,
  .ScaleCols = (__typeof__(&bfMatProductScaleCols))bfMatProductScaleCols,
  .Mul = (__typeof__(&bfMatProductMul))bfMatProductMul,
  .Solve = (__typeof__(&bfMatSolve))bfMatProductSolve,
};

BfMat *bfMatProductCopy(BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfMatProduct const *matProduct = NULL;
  BfMatProduct *matProductCopy = NULL;
  BfMat const *factor = NULL;
  BfMat *factorCopy = NULL;

  matProduct = bfMatConstToMatProductConst(mat);
  HANDLE_ERROR();

  matProductCopy = bfMatProductNew();
  HANDLE_ERROR();

  bfMatProductInit(matProductCopy);
  HANDLE_ERROR();

  BfSize numFactors = bfMatProductNumFactors(matProduct);

  for (BfSize i = 0; i < numFactors; ++i) {
    factor = bfMatProductGetFactorConst(matProduct, i);
    factorCopy = bfMatCopy(factor);
    HANDLE_ERROR();
    bfMatProductPostMultiply(matProductCopy, factorCopy);
  }

  END_ERROR_HANDLING() {
    bfMatDelete(&factorCopy);
    bfMatProductDeinitAndDealloc(&matProductCopy);
  }

  return bfMatProductToMat(matProductCopy);
}

void bfMatProductDelete(BfMat **mat) {
  bfMatProductDeinitAndDealloc((BfMatProduct **)mat);
}

BfType bfMatProductGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_PRODUCT;
}

bool bfMatProductInstanceOf(BfMat const *mat, BfType type) {
  return bfTypeDerivedFrom(bfMatGetType(mat), type);
}

BfSize bfMatProductGetNumRows(BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfMat const *factor = NULL;
  BfSize numRows = BF_SIZE_BAD_VALUE;

  if (bfMatIsTransposed(mat))
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  factor = bfMatProductGetFactor((BfMatProduct *)mat, 0);
  HANDLE_ERROR();

  numRows = bfMatGetNumRows(factor);

  END_ERROR_HANDLING() {}

  return numRows;
}

BfSize bfMatProductGetNumCols(BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfMatProduct const *matProduct = NULL;
  BfMat const *factor = NULL;
  BfSize numCols = BF_SIZE_BAD_VALUE;

  matProduct = bfMatConstToMatProductConst(mat);
  HANDLE_ERROR();

  if (bfMatIsTransposed(mat))
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfSize numFactors = bfMatProductNumFactors(matProduct);

  factor = bfMatProductGetFactorConst(matProduct, numFactors - 1);
  HANDLE_ERROR();

  numCols = bfMatGetNumCols(factor);

  END_ERROR_HANDLING() {}

  return numCols;
}

void bfMatProductScaleCols(BfMat *mat, BfVec const *vec) {
  BEGIN_ERROR_HANDLING();

  BfMatProduct *matProduct = bfMatToMatProduct(mat);
  HANDLE_ERROR();

  BfSize numFactors = bfMatProductNumFactors(matProduct);
  if (numFactors == 0)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMat *factor = bfMatProductGetFactor(matProduct, numFactors - 1);

  bfMatScaleCols(factor, vec);

  END_ERROR_HANDLING() {}
}

BfMat *bfMatProductMul(BfMat const *mat, BfMat const *otherMat) {
  BEGIN_ERROR_HANDLING();

  BfMatProduct const *matProduct = NULL;
  BfMat const *factor = NULL;
  BfMat *prev = NULL;
  BfMat *result = NULL;

  matProduct = bfMatConstToMatProductConst(mat);
  HANDLE_ERROR();

  BfSize numFactors = bfMatProductNumFactors(matProduct);

  BfSize i = numFactors - 1;

  factor = bfMatProductGetFactorConst(matProduct, i);

  prev = bfMatMul(factor, otherMat);
  HANDLE_ERROR();

  while (i > 0) {
    factor = bfMatProductGetFactorConst(matProduct, --i);

    result = bfMatMul(factor, prev);
    HANDLE_ERROR();

    bfMatDelete(&prev);

    prev = result;
  }

  END_ERROR_HANDLING() {
    bfMatDelete(&prev);
    bfMatDelete(&result);
  }

  return result;
}

BfMat *bfMatProductSolve(BfMatProduct const *matProduct, BfMat const *otherMat) {
  BEGIN_ERROR_HANDLING();

  BfSize numFactors = bfMatProductNumFactors(matProduct);
  if (numFactors == 0)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMat const *factor = bfMatProductGetFactorConst(matProduct, 0);

  BfMat *result = bfMatSolve(factor, otherMat);
  HANDLE_ERROR();

  for (BfSize i = 1; i < numFactors; ++i) {
    factor = bfMatProductGetFactorConst(matProduct, i);

    BfMat *_ = bfMatSolve(factor, result);
    HANDLE_ERROR();

    bfMatDelete(&result);

    result = _;
  }

  END_ERROR_HANDLING() {
    assert(false);
  }

  return result;
}

/** Upcasting: */

BfMat *bfMatProductToMat(BfMatProduct *matProduct) {
  return &matProduct->super;
}

/** Downcasting: */

BfMatProduct *bfMatToMatProduct(BfMat *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_PRODUCT)) {
    bfSetError(BF_ERROR_RUNTIME_ERROR);
    return NULL;
  } else {
    return (BfMatProduct *)mat;
  }
}

BfMatProduct const *bfMatConstToMatProductConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_PRODUCT)) {
    bfSetError(BF_ERROR_RUNTIME_ERROR);
    return NULL;
  } else {
    return (BfMatProduct const *)mat;
  }
}

/** Implementation: MatProduct */

BfMatProduct *bfMatProductNew() {
  BEGIN_ERROR_HANDLING();

  BfMatProduct *prod = bfMemAlloc(1, sizeof(BfMatProduct));
  if (prod == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return prod;
}

void bfMatProductInit(BfMatProduct *mat) {
  BEGIN_ERROR_HANDLING();

  /* We don't store the number of rows or columns in `mat->super`
   * since we always look up the number of rows and columns from the
   * leftmost and rightmost factors at runtime. */
  bfMatInit(&mat->super, &MAT_VTABLE, BF_SIZE_BAD_VALUE, BF_SIZE_BAD_VALUE);

  mat->factorArr = bfGetUninitializedPtrArray();

  bfInitPtrArray(&mat->factorArr, /* capacity: */ 4);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfPtrArrayDeinit(&mat->factorArr);
}

void bfMatProductDeinit(BfMatProduct *prod) {
  for (BfSize i = 0; i < bfPtrArraySize(&prod->factorArr); ++i) {
    BfMat *mat = bfPtrArrayGet(&prod->factorArr, i);
    bfMatDelete(&mat);
  }

  bfPtrArrayDeinit(&prod->factorArr);
}

void bfMatProductDealloc(BfMatProduct **prod) {
  free(*prod);
  *prod = NULL;
}

void bfMatProductDeinitAndDealloc(BfMatProduct **prod) {
  bfMatProductDeinit(*prod);
  bfMatProductDealloc(prod);
}

BfSize bfMatProductNumFactors(BfMatProduct const *prod) {
  return bfPtrArraySize(&prod->factorArr);
}

BfMat *bfMatProductGetFactor(BfMatProduct *prod, BfSize i) {
  return bfPtrArrayGet(&prod->factorArr, i);
}

BfMat const *bfMatProductGetFactorConst(BfMatProduct const *prod, BfSize i) {
  return bfPtrArrayGet(&prod->factorArr, i);
}

void bfMatProductPostMultiply(BfMatProduct *prod, BfMat *mat) {
  // TODO: verify that shapes are compatible
  bfPtrArrayAppend(&prod->factorArr, mat);
}

BfMat *bfMatProductPopLastFactor(BfMatProduct *matProduct) {
  BEGIN_ERROR_HANDLING();

  BfMat *lastFactor = NULL;
  BfSize numFactors = bfMatProductNumFactors(matProduct);

  if (numFactors == 0)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  lastFactor = bfPtrArrayPopLast(&matProduct->factorArr);

  END_ERROR_HANDLING()
    lastFactor = NULL;

  return lastFactor;
}
