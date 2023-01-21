#pragma once

#include "geom.h"
#include "mat.h"

BfReal bfPoint2Dist(BfPoint2 const p, BfPoint2 const q);

BfReal bfPoint3Dist(BfPoint3 const p, BfPoint3 const q);
void bfPoint3Sub(BfPoint3 const v, BfPoint3 const u, BfVector3 uv);
void bfPoint3GetPointOnRay(BfPoint3 const r0, BfVector3 const dr, BfReal t, BfPoint3 rt);
void bfPoint3Copy(BfPoint3 x, BfPoint3 const y);

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
BfBbox2 bfPoints2GetBoundingBox(BfPoints2 const *points);
void bfGetPointsByIndex(BfPoints2 const *points,
                                BfSize numInds, BfSize const *inds,
                                BfPoints2 *indexedPoints);
void bfPrintPoints2(BfPoints2 const *points);
void bfSavePoints2(BfPoints2 const *points, char const *path);
BfReal *bfPoints2PairwiseDists(BfPoints2 const *X, BfPoints2 const *Y);

struct BfPoints3 {
  BfPoint3 *data;
  BfSize size;
};

void bfPoints3InitEmpty(BfPoints3 *points, BfSize numPoints);
void bfPoints3InitFromBinaryFile(BfPoints3 *points, char const *path);
void bfPoints3Deinit(BfPoints3 *points);
BfBoundingBox3 bfPoints3GetBoundingBox(BfPoints3 const *points);
void bfPoints3GetByIndex(BfPoints3 const *points, BfSize numInds, BfSize const *inds, BfPoints3 *indexedPoints);
