#pragma once

#include "def.h"
#include "mat.h"

typedef struct BfBbox2 {
  BfPoint2 min, max;
} BfBbox2;

typedef struct BfCircle2 {
  BfReal r;
  BfPoint2 center;
} BfCircle2;

bool bfBbox2ContainsPoint(BfBbox2 const *bbox, BfPoint2 const point);
BfReal bfPoint2Dist(BfPoint2 const p, BfPoint2 const q);
enum BfError bfSamplePointsOnCircle2(BfCircle2 const *circ, BfMat *X);
