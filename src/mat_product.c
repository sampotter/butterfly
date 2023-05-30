#include <bf/mat_product.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>

/** Interface: MatProduct */

static BfMatVtable MAT_VTABLE = {
  .GetView = (__typeof__(&bfMatGetView))bfMatProductGetView,
  .Copy = (__typeof__(&bfMatProductCopy))bfMatProductCopy,
  .Steal = (__typeof__(&bfMatSteal))bfMatProductSteal,
  .Delete = (__typeof__(&bfMatProductDelete))bfMatProductDelete,
  .GetType = (__typeof__(&bfMatProductGetType))bfMatProductGetType,
  .NumBytes = (__typeof__(&bfMatNumBytes))bfMatProductNumBytes,
  .Dump = (__typeof__(&bfMatDump))bfMatProductDump,
  .GetNumRows = (__typeof__(&bfMatProductGetNumRows))bfMatProductGetNumRows,
  .GetNumCols = (__typeof__(&bfMatProductGetNumCols))bfMatProductGetNumCols,
  .ScaleCols = (__typeof__(&bfMatProductScaleCols))bfMatProductScaleCols,
  .Mul = (__typeof__(&bfMatProductMul))bfMatProductMul,
  .MulVec = (__typeof__(&bfMatMulVec))bfMatProductMulVec,
  .RmulVec = (__typeof__(&bfMatRmulVec))bfMatProductRmulVec,
  .Solve = (__typeof__(&bfMatSolve))bfMatProductSolve,
};

BfMat *bfMatProductGetView(BfMatProduct *matProduct) {
  BF_ERROR_BEGIN();

  BfMatProduct *matProductView = bfMatProductNew();
  HANDLE_ERROR();

  *matProductView = *matProduct;

  BfMat *matView = bfMatProductToMat(matProductView);

  matView->props |= BF_MAT_PROPS_VIEW;

  BF_ERROR_END()
    matView = NULL;

  return matView;
}

BfMat *bfMatProductCopy(BfMat const *mat) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END() {
    bfMatDelete(&factorCopy);
    bfMatProductDeinitAndDealloc(&matProductCopy);
  }

  return bfMatProductToMat(matProductCopy);
}

BfMat *bfMatProductSteal(BfMatProduct *matProduct) {
  BF_ERROR_BEGIN();

  BfMat *mat = bfMatProductToMat(matProduct);

  if (bfMatIsView(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatProduct *matProductNew = bfMatProductNew();
  HANDLE_ERROR();

  *matProductNew = *matProduct;

  mat->props |= BF_MAT_PROPS_VIEW;

  BF_ERROR_END() {
    BF_DIE();
  }

  return bfMatProductToMat(matProductNew);
}

void bfMatProductDelete(BfMat **mat) {
  bfMatProductDeinitAndDealloc((BfMatProduct **)mat);
}

BfType bfMatProductGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_PRODUCT;
}

BfSize bfMatProductNumBytes(BfMatProduct const *matProduct) {
  BfSize numBytes = 0;
  for (BfSize i = 0; i < bfMatProductNumFactors(matProduct); ++i) {
    BfMat const *factor = bfMatProductGetFactorConst(matProduct, i);
    numBytes += bfMatNumBytes(factor);
  }
  return numBytes;
}

void bfMatProductDump(BfMatProduct const *matProduct, FILE *fp) {
  BF_ERROR_BEGIN();

  /* Serialize the number of factors: */
  BfSize numFactors = bfMatProductNumFactors(matProduct);
  fwrite(&numFactors, sizeof(BfSize), 1, fp);

  /* Recursively write each of the factors: */
  for (BfSize k = 0; k < numFactors; ++k) {
    BfMat const *factor = bfMatProductGetFactorConst(matProduct, k);
    bfMatDump(factor, fp);
    HANDLE_ERROR();
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

bool bfMatProductInstanceOf(BfMat const *mat, BfType type) {
  return bfTypeDerivedFrom(bfMatGetType(mat), type);
}

BfSize bfMatProductGetNumRows(BfMat const *mat) {
  BF_ERROR_BEGIN();

  BfMat const *factor = NULL;
  BfSize numRows = BF_SIZE_BAD_VALUE;

  if (bfMatIsTransposed(mat))
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  factor = bfMatProductGetFactor((BfMatProduct *)mat, 0);
  HANDLE_ERROR();

  numRows = bfMatGetNumRows(factor);

  BF_ERROR_END() {}

  return numRows;
}

BfSize bfMatProductGetNumCols(BfMat const *mat) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END() {}

  return numCols;
}

void bfMatProductScaleCols(BfMat *mat, BfVec const *vec) {
  BF_ERROR_BEGIN();

  BfMatProduct *matProduct = bfMatToMatProduct(mat);
  HANDLE_ERROR();

  BfSize numFactors = bfMatProductNumFactors(matProduct);
  if (numFactors == 0)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMat *factor = bfMatProductGetFactor(matProduct, numFactors - 1);

  bfMatScaleCols(factor, vec);

  BF_ERROR_END() {}
}

BfMat *bfMatProductMul(BfMat const *mat, BfMat const *otherMat) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END() {
    bfMatDelete(&prev);
    bfMatDelete(&result);
  }

  return result;
}

BfVec *bfMatProductMulVec(BfMatProduct const *matProduct, BfVec const *vec) {
  BF_ERROR_BEGIN();

  BfMat const *factor = NULL;
  BfVec *prev = NULL;
  BfVec *result = NULL;

  BfSize numFactors = bfMatProductNumFactors(matProduct);

  BfSize i = numFactors - 1;

  factor = bfMatProductGetFactorConst(matProduct, i);

  prev = bfMatMulVec(factor, vec);
  HANDLE_ERROR();

  while (i > 0) {
    factor = bfMatProductGetFactorConst(matProduct, --i);

    result = bfMatMulVec(factor, prev);
    HANDLE_ERROR();

    bfVecDelete(&prev);

    prev = result;
  }

  BF_ERROR_END() {
    bfVecDelete(&prev);
    bfVecDelete(&result);
  }

  return result;
}

BfVec *bfMatProductRmulVec(BfMatProduct const *matProduct, BfVec const *vec) {
  BF_ERROR_BEGIN();

  BfMat const *factor = NULL;
  BfVec *prev = NULL;
  BfVec *result = NULL;

  BfSize numFactors = bfMatProductNumFactors(matProduct);

  BfSize i = 0;

  factor = bfMatProductGetFactorConst(matProduct, i);

  prev = bfMatRmulVec(factor, vec);
  HANDLE_ERROR();

  while (++i < numFactors) {
    factor = bfMatProductGetFactorConst(matProduct, i);

    result = bfMatRmulVec(factor, prev);
    HANDLE_ERROR();

    bfVecDelete(&prev);

    prev = result;
  }

  BF_ERROR_END() {
    bfVecDelete(&prev);
    bfVecDelete(&result);
  }

  return result;
}

BfMat *bfMatProductSolve(BfMatProduct const *matProduct, BfMat const *otherMat) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END() {
    BF_DIE();
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
  BF_ERROR_BEGIN();

  BfMatProduct *prod = bfMemAlloc(1, sizeof(BfMatProduct));
  if (prod == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BF_ERROR_END() {}

  return prod;
}

void bfMatProductInit(BfMatProduct *mat) {
  BF_ERROR_BEGIN();

  /* We don't store the number of rows or columns in `mat->super`
   * since we always look up the number of rows and columns from the
   * leftmost and rightmost factors at runtime. */
  bfMatInit(&mat->super, &MAT_VTABLE, BF_SIZE_BAD_VALUE, BF_SIZE_BAD_VALUE);

  mat->factorArr = bfGetUninitializedPtrArray();

  bfInitPtrArray(&mat->factorArr, /* capacity: */ 4);
  HANDLE_ERROR();

  BF_ERROR_END()
    bfPtrArrayDeinit(&mat->factorArr);
}

void bfMatProductDeinit(BfMatProduct *prod) {
  if (!bfMatIsView(&prod->super)) {
    for (BfSize i = 0; i < bfPtrArraySize(&prod->factorArr); ++i) {
      BfMat *mat = bfPtrArrayGet(&prod->factorArr, i);
      bfMatDelete(&mat);
    }
    bfPtrArrayDeinit(&prod->factorArr);
  }
}

void bfMatProductDealloc(BfMatProduct **prod) {
  bfMemFree(*prod);
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
  BF_ERROR_BEGIN();

  BfMat *lastFactor = NULL;
  BfSize numFactors = bfMatProductNumFactors(matProduct);

  if (numFactors == 0)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  lastFactor = bfPtrArrayPopLast(&matProduct->factorArr);

  BF_ERROR_END()
    lastFactor = NULL;

  return lastFactor;
}
