#pragma once

#include "geom.h"
#include "mat_dense_real.h"

void bfVector3Copy(BfVector3 u, BfVector3 const v);
void bfVector3Scale(BfVector3 u, BfReal alpha);
BfReal bfVector3Norm(BfVector3 const u);
BfReal bfVector3Dot(BfVector3 const u, BfVector3 const v);
void bfVector3Cross(BfVector3 const u, BfVector3 const v, BfVector3 w);

struct BfVectors2 {
  BfVector2 *data;
  BfSize size;
};

BfVectors2 const *bfVectors2ConstViewFromMat(BfMat const *mat);
BfVectors2 const *bfVectors2ConstViewFromMatDenseReal(BfMatDenseReal const *matDenseReal);
BfVectors2 bfGetUninitializedVectors2();
void bfInitEmptyVectors2(BfVectors2 *vectors, BfSize numVectors);
void bfReadVectors2FromFile(char const *path, BfVectors2 *vectors);
void bfFreeVectors2(BfVectors2 *vectors);
void bfGetVectorsByIndex(BfVectors2 const *vectors,
                        BfSize numInds, BfSize const *inds,
                        BfVectors2 *indexedVectors);
void bfSaveVectors2(BfVectors2 const *vectors, char const *path);

struct BfVectors3 {
  BfVector3 *data;
  BfSize size;
};

void bfVectors3InitEmpty(BfVectors3 *vectors, BfSize numVectors);
void bfVectors3Deinit(BfVectors3 *vectors);
void bfVectors3GetByIndex(BfVectors3 const *vectors, BfSize numInds, BfSize const *inds, BfVectors3 *indexedVectors);
