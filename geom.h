#pragma once

#include "def.h"
#include "mat.h"

typedef struct BfPoints2 {
  BfPoint2 *data;
  BfSize size;
} BfPoints2;

typedef struct BfBbox2 {
  BfPoint2 min, max;
} BfBbox2;

typedef struct BfCircle2 {
  BfReal r;
  BfPoint2 center;
} BfCircle2;

BfReal *bfPoints2PairwiseDists(BfPoints2 const *X, BfPoints2 const *Y);

bool bfBbox2ContainsPoint(BfBbox2 const *bbox, BfPoint2 const point);
bool bfBbox2ContainsPoints(BfBbox2 const *bbox, BfPoints2 const *points);
BfReal bfPoint2Dist(BfPoint2 const p, BfPoint2 const q);
BfPoints2 bfSamplePointsOnCircle2(BfCircle2 const *circ, BfSize numPoints);
bool bfCircle2ContainsPoint(BfCircle2 const *circ, BfPoint2 const point);
bool bfCircle2ContainsPoints(BfCircle2 const *circ, BfPoints2 const *points);

BfBbox2 bfGetEmptyBbox2();
bool bfBbox2IsEmpty(BfBbox2 const *bbox);

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
