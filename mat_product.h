#pragma once

#include "mat.h"
#include "ptr_array.h"

typedef struct BfMatProduct {
  BfMat super;
  BfPtrArray factorArr;
} BfMatProduct;

BfMatProduct *bfMatProductNew();
void bfMatProductInit(BfMatProduct *prod);
BfMat *bfMatProductGetMatPtr();
BfSize bfMatProductNumFactors(BfMatProduct *prod);
BfMat *bfMatProductGetFactor(BfMatProduct *prod, BfSize i);
void bfMatProductPostMultiply(BfMatProduct *prod, BfMat *mat);

/* BfMat interface: */
void bfMatProductDeinit(BfMatProduct *mat);
void bfMatProductDelete(BfMatProduct **mat);
void bfMatProductDeinitAndDelete(BfMatProduct **mat);
BfMatType bfMatProductGetType(BfMatProduct const *mat);
BfSize bfMatProductNumBytes(BfMatProduct const *mat);
void bfMatProductSave(BfMatProduct const *mat, char const *path);
BfMat *bfMatProductMul(BfMatProduct const *op1, BfMat const *op2);
BfMat *bfMatProductLstSq(BfMatProduct const *op1, BfMat const *op2);
