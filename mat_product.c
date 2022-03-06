#include "mat_product.h"

#include <assert.h>
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

void bfMatProductInit(BfMatProduct *prod) {
  BEGIN_ERROR_HANDLING();

  prod->factorArr = bfGetUninitializedPtrArray();

  bfInitPtrArray(&prod->factorArr, /* capacity: */ 4);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfPtrArrayDeinit(&prod->factorArr);
}

void bfMatProductDeinit(BfMatProduct *prod) {
  for (BfSize i = 0; i < bfPtrArraySize(&prod->factorArr); ++i) {
    BfMat *mat = bfPtrArrayGet(&prod->factorArr, i);
    bfMatDeinitAndDelete(&mat);
  }

  bfPtrArrayDeinit(&prod->factorArr);
}

void bfMatProductDelete(BfMatProduct **prod) {
  free(*prod);
  *prod = NULL;
}

void bfMatProductDeinitAndDelete(BfMatProduct **prod) {
  bfMatProductDeinit(*prod);
  bfMatProductDelete(prod);
}

BfSize bfMatProductNumFactors(BfMatProduct *prod) {
  return bfPtrArraySize(&prod->factorArr);
}

BfMat *bfMatProductGetFactor(BfMatProduct *prod, BfSize i) {
  return bfPtrArrayGet(&prod->factorArr, i);
}

void bfMatProductPostMultiply(BfMatProduct *prod, BfMat *mat) {
  assert(false); // TODO: verify that shapes are compatible

  bfPtrArrayAppend(&prod->factorArr, mat);
}
