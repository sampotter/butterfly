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
BfSize bfMatProductGetNumRows(BfMat const *mat);
BfSize bfMatProductGetNumCols(BfMat const *mat);
void bfMatProductScaleCols(BfMat *mat, BfVec const *vec);
BfMat *bfMatProductMul(BfMat const *mat, BfMat const *otherMat);
BfVec *bfMatProductMulVec(BfMatProduct const *matProduct, BfVec const *vec);
BfVec *bfMatProductRmulVec(BfMatProduct const *matProduct, BfVec const *vec);
BfMat *bfMatProductSolve(BfMatProduct const *matProduct, BfMat const *otherMat);

/** Implementation: MatProduct */

typedef struct BfMatProduct {
  BfMat super;
  BfPtrArray factorArr;
} BfMatProduct;

BfMat *bfMatProductToMat(BfMatProduct *matProduct);

BfMatProduct *bfMatToMatProduct(BfMat *mat);
BfMatProduct const *bfMatConstToMatProductConst(BfMat const *mat);

BfMatProduct *bfMatProductNew();
void bfMatProductInit(BfMatProduct *prod);
void bfMatProductDeinit(BfMatProduct *prod);
void bfMatProductDealloc(BfMatProduct **prod);
void bfMatProductDeinitAndDealloc(BfMatProduct **prod);
BfSize bfMatProductNumFactors(BfMatProduct const *prod);
BfMat *bfMatProductGetFactor(BfMatProduct *prod, BfSize i);
BfMat const *bfMatProductGetFactorConst(BfMatProduct const *prod, BfSize i);
void bfMatProductPostMultiply(BfMatProduct *prod, BfMat *mat);
BfMat *bfMatProductPopLastFactor(BfMatProduct *prod);
