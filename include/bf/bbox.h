#pragma once

#include "geom.h"

typedef struct BfBbox2 {
  BfPoint2 min, max;
} BfBbox2;

BfBbox2 bfGetEmptyBbox2();
bool bfBbox2IsEmpty(BfBbox2 const *bbox);
bool bfBbox2ContainsPoint(BfBbox2 const *bbox, BfPoint2 const point);
bool bfBbox2ContainsPoints(BfBbox2 const *bbox, BfPoints2 const *points);
