#pragma once

#include "geom.h"

typedef struct BfCircle2 {
  BfReal r;
  BfPoint2 center;
} BfCircle2;

BfPoints2 bfSamplePointsOnCircle2(BfCircle2 const *circ, BfSize numPoints);
bool bfCircle2ContainsPoint(BfCircle2 const *circ, BfPoint2 const point);
bool bfCircle2ContainsPoints(BfCircle2 const *circ, BfPoints2 const *points);
