#pragma once

#include "geom.h"
#include "sphere.h"

typedef struct BfBbox2 {
  BfPoint2 min, max;
} BfBbox2;

BfBbox2 bfGetEmptyBbox2();
bool bfBbox2IsEmpty(BfBbox2 const *bbox);
void bfBbox2RescaleToSquare(BfBbox2 *bbox);
bool bfBbox2ContainsPoint(BfBbox2 const *bbox, BfPoint2 const point);
bool bfBbox2ContainsPoints(BfBbox2 const *bbox, BfPoints2 const *points);
void bfBbox2GetCenter(BfBbox2 const *bbox, BfPoint2 center);

typedef struct BfBoundingBox3 {
  BfPoint3 min, max;
} BfBoundingBox3;

void bfBoundingBox3InitEmpty(BfBoundingBox3 *boundingBox);
bool bfBoundingBox3IsEmpty(BfBoundingBox3 const *boundingBox);
void bfBoundingBox3RescaleToCube(BfBoundingBox3 *boundingBox);
bool bfBoundingBox3ContainsPoint(BfBoundingBox3 const *boundingBox, BfPoint3 const point);
bool bfBoundingBox3ContainsPoints(BfBoundingBox3 const *boundingBox, BfPoints3 const *points);
void bfBoundingBox3GetCenter(BfBoundingBox3 const *boundingBox, BfPoint3 center);
BfSphere bfBoundingBox3GetBoundingSphere(BfBoundingBox3 const *boundingBox);
BfReal bfBoundingBox3GetVolume(BfBoundingBox3 const *boundingBox);
void bfBoundingBox3GetSemiLengths(BfBoundingBox3 const *boundingBox, BfVector3 semiLengths);
BfReal bfBoundingBox3GetDistanceToPoint(BfBoundingBox3 const *boundingBox, BfPoint3 const point);
