#pragma once

#include "mat.h"
#include "ptr_array.h"

/** Interface: Mat */

BfMat *bfMatProductGetView(BfMatProduct *matProduct);
BfMat *bfMatProductCopy(BfMat const *mat);
BfMat *bfMatProductSteal(BfMatProduct *matProduct);
void bfMatProductDelete(BfMat **mat);
BfType bfMatProductGetType(BfMat const *mat);
BfSize bfMatProductNumBytes(BfMatProduct const *matProduct);
void bfMatProductDump(BfMatProduct const *matProduct, FILE *fp);
BfSize bfMatProductGetNumRows(BfMatProduct const *matProduct);
BfSize bfMatProductGetNumCols(BfMat const *mat);
void bfMatProductScaleCols(BfMat *mat, BfVec const *vec);
BfMat *bfMatProductMul(BfMat const *mat, BfMat const *otherMat);
BfVec *bfMatProductMulVec(BfMatProduct const *matProduct, BfVec const *vec);
BfMat *bfMatProductRmul(BfMatProduct const *matProduct, BfMat const *mat);
BfVec *bfMatProductRmulVec(BfMatProduct const *matProduct, BfVec const *vec);
BfMat *bfMatProductSolve(BfMatProduct const *matProduct, BfMat const *otherMat);
void bfMatProductTranspose(BfMatProduct *matProduct);

/** Implementation: MatProduct */

typedef struct BfMatProduct {
  BfMat super;
  BfPtrArray factorArr;
} BfMatProduct;

BfMat *bfMatProductToMat(BfMatProduct *matProduct);
BfMat const *bfMatProductConstToMatConst(BfMatProduct const *matProduct);

BfMatProduct *bfMatToMatProduct(BfMat *mat);
BfMatProduct const *bfMatConstToMatProductConst(BfMat const *mat);

BfMatProduct *bfMatProductNew(void);
void bfMatProductInit(BfMatProduct *prod);
void bfMatProductDeinit(BfMatProduct *prod);
void bfMatProductDealloc(BfMatProduct **prod);
void bfMatProductDeinitAndDealloc(BfMatProduct **prod);
BfSize bfMatProductNumFactors(BfMatProduct const *prod);
BfMat *bfMatProductGetFactor(BfMatProduct *prod, BfSize i);
BfMat const *bfMatProductGetFactorConst(BfMatProduct const *prod, BfSize i);
void bfMatProductPostMultiply(BfMatProduct *prod, BfMat *mat);
BfMat *bfMatProductPopLastFactor(BfMatProduct *prod);
