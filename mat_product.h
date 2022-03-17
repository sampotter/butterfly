#pragma once

#include "mat.h"
#include "ptr_array.h"

typedef struct BfMatProduct {
  BfMat super;
  BfPtrArray factorArr;
} BfMatProduct;

BF_DECLARE_INTERFACE_MAT(MatProduct);

BfMatProduct *bfMatProductNew();
void bfMatProductInit(BfMatProduct *prod);
BfMat *bfMatProductGetMatPtr(BfMatProduct *);
BfMat const *bfMatProductGetMatConstPtr(BfMatProduct const *);
BfSize bfMatProductNumFactors(BfMatProduct *prod);
BfMat *bfMatProductGetFactor(BfMatProduct *prod, BfSize i);
void bfMatProductPostMultiply(BfMatProduct *prod, BfMat *mat);
