#pragma once

#include "geom.h"
#include "mat.h"

BfReal bfPoint2Dist(BfPoint2 const p, BfPoint2 const q);

struct BfPoints2 {
  BfPoint2 *data;
  BfSize size;
};

BfPoints2 const *bfPoints2ConstViewFromMat(BfMat const *mat);
BfPoints2 const *bfPoints2ConstViewFromMatDenseReal(BfMatDenseReal const *matDenseReal);
BfPoints2 bfGetUninitializedPoints2();
void bfInitEmptyPoints2(BfPoints2 *points, BfSize numPoints);
void bfReadPoints2FromFile(char const *path, BfPoints2 *points);
void bfFreePoints2(BfPoints2 *points);
bool bfPoints2Initialized(BfPoints2 const *points);
BfBbox2 bfGetPoints2BoundingBox(BfPoints2 const *points);
void bfGetPointsByIndex(BfPoints2 const *points,
                                BfSize numInds, BfSize const *inds,
                                BfPoints2 *indexedPoints);
void bfPrintPoints2(BfPoints2 const *points);
void bfSavePoints2(BfPoints2 const *points, char const *path);
BfReal *bfPoints2PairwiseDists(BfPoints2 const *X, BfPoints2 const *Y);
