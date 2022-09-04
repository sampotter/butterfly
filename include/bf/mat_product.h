#pragma once

#include "mat.h"
#include "ptr_array.h"

typedef struct BfMatProduct {
  BfMat super;
  BfPtrArray factorArr;
} BfMatProduct;

#define INTERFACE BF_INTERFACE_Mat
BF_DECLARE_INTERFACE(MatProduct);
#undef INTERFACE

BfMatProduct *bfMatProductNew();
void bfMatProductInit(BfMatProduct *prod);
void bfMatProductDeinit(BfMatProduct *prod);
void bfMatProductDealloc(BfMatProduct **prod);
void bfMatProductDeinitAndDealloc(BfMatProduct **prod);
BfMat *bfMatProductGetMatPtr(BfMatProduct *);
BfMat const *bfMatProductGetMatConstPtr(BfMatProduct const *);
BfSize bfMatProductNumFactors(BfMatProduct *prod);
BfMat *bfMatProductGetFactor(BfMatProduct *prod, BfSize i);
void bfMatProductPostMultiply(BfMatProduct *prod, BfMat *mat);
