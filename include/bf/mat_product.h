#pragma once

#include "mat.h"
#include "ptr_array.h"

typedef struct BfMatProduct {
  BfMat super;
  BfPtrArray factorArr;
} BfMatProduct;

#define INTERFACE BF_INTERFACE_Mat
BF_DECLARE_INTERFACE(MatProduct)
#undef INTERFACE

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
