#pragma once

#include "geom.h"
#include "mat_dense_real.h"

struct BfVectors2 {
  BfVector2 *data;
  BfSize size;
};

BfVectors2 const *bfVectors2ConstViewFromMat(BfMat const *mat);
BfVectors2 const *bfVectors2ConstViewFromMatDenseReal(BfMatDenseReal const *matDenseReal);
BfVectors2 bfGetUninitializedVectors2();
void bfInitEmptyVectors2(BfVectors2 *vectors, BfSize numVectors);
void bfFreeVectors2(BfVectors2 *vectors);
void bfGetVectorsByIndex(BfVectors2 const *vectors,
                        BfSize numInds, BfSize const *inds,
                        BfVectors2 *indexedVectors);
