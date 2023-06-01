#pragma once

#include "geom.h"

struct BfCircle {
  BfReal r;
  BfPoint2 center;
};

BfPoints2 bfCircle2SamplePoints(BfCircle const *circ, BfSize numPoints);
BfVectors2 *bfCircle2SampleUnitNormals(BfCircle const *circ, BfSize numPoints);
bool bfCircle2ContainsPoint(BfCircle const *circ, BfPoint2 const point);
bool bfCircle2ContainsPoints(BfCircle const *circ, BfPoints2 const *points);
