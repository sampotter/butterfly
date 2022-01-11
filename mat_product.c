#include "mat_product.h"

#include <stdlib.h>

#include "error.h"
#include "error_macros.h"

BfMatProduct *bfMatProductNew() {
  BEGIN_ERROR_HANDLING();

  BfMatProduct *prod = malloc(sizeof(BfMatProduct));
  if (prod == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return prod;
}

void bfMatProductInit(BfMatProduct *prod, BfSize numFactors) {
  bfInitPtrArray(prod->factorArr, /* capacity: */ 4);
}

void bfMatProductDeinit(BfMatProduct *prod) {
  for (BfSize i = 0; i < bfPtrArraySize(prod->factorArr); ++i) {
    BfMat *mat = bfPtrArrayGet(prod->factorArr, i);
    BfMatDeinitAndDelete(&mat);
  }

  bfPtrArrayDeinitAndDelete(&prod->factorArr);
}

void bfMatProductDelete(BfMatProduct **prod) {
  free(*prod);
  *prod = NULL;
}

void bfMatProductDeinitAndDelete(BfMatProduct **prod) {
  bfMatProductDeinit(*prod);
  bfMatProductDelete(prod);
}

void bfMatProductPostMultiply(BfMatProduct *prod, BfMat *mat) {
  bfPtrArrayAppend(prod, mat);
}
