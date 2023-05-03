#pragma once

#include "def.h"
#include "types.h"

struct BfIndexedMat {
  /* The leading row index for `mat`. */
  BfSize i0;

  /* The leading column index for `mat`. */
  BfSize j0;

  /* The indexed matrix. */
  BfMat *mat;
};

BfIndexedMat *bfIndexedMatNewFromMat(BfSize i0, BfSize j0, BfMat *mat, BfPolicy policy);
void bfIndexedMatInitFromMat(BfIndexedMat *indexedMat, BfSize i0, BfSize j0, BfMat *mat, BfPolicy policy);
void bfIndexedMatDeinit(BfIndexedMat *indexedMat);
void bfIndexedMatDealloc(BfIndexedMat **indexedMat);
void bfIndexedMatDelete(BfIndexedMat **indexedMat);
