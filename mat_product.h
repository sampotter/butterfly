#pragma once

#include "mat.h"
#include "ptr_array.h"

typedef struct BfMatProduct {
  BfPtrArray factorArr;
} BfMatProduct;

BfMatProduct *bfMatProductNew();
void bfMatProductInit(BfMatProduct *prod);
void bfMatProductDeinit(BfMatProduct *prod);
void bfMatProductDelete(BfMatProduct **prod);
void bfMatProductDeinitAndDelete(BfMatProduct **prod);
BfSize bfMatProductNumFactors(BfMatProduct *prod);
BfMat *bfMatProductGetFactor(BfMatProduct *prod, BfSize i);
void bfMatProductPostMultiply(BfMatProduct *prod, BfMat *mat);
