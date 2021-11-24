#pragma once

#include "def.h"

typedef struct BfBbox2 {
  BfPoint2 min, max;
} BfBbox2;

bool bfBbox2ContainsPoint(BfBbox2 const *bbox, BfPoint2 const point);
BfReal bfPoint2Dist(BfPoint2 const p, BfPoint2 const q);
