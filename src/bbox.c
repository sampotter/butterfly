#include <bf/bbox.h>

#include <math.h>

#include <bf/points.h>

BfBbox2 bfGetEmptyBbox2() {
  return (BfBbox2) {
    .min = { INFINITY,  INFINITY},
    .max = {-INFINITY, -INFINITY}
  };
}

bool bfBbox2IsEmpty(BfBbox2 const *bbox) {
  return bbox->min[0] >= bbox->max[0] && bbox->min[1] >= bbox->max[1];
}

void bfBbox2RescaleToSquare(BfBbox2 *bbox) {
  BfReal w = bbox->max[0] - bbox->min[0];
  BfReal h = bbox->max[1] - bbox->min[1];
  if (w > h) {
    BfReal c = (bbox->min[1] + bbox->max[1])/2;
    bbox->min[1] = w*(bbox->min[1] - c)/h + c;
    bbox->max[1] = w*(bbox->max[1] - c)/h + c;
  } else {
    BfReal c = (bbox->min[0] + bbox->max[0])/2;
    bbox->min[0] = h*(bbox->min[0] - c)/w + c;
    bbox->max[0] = h*(bbox->max[0] - c)/w + c;
  }
}

bool bfBbox2ContainsPoint(BfBbox2 const *bbox, BfPoint2 const point) {
  return bbox->min[0] <= point[0] && point[0] <= bbox->max[0]
      && bbox->min[1] <= point[1] && point[1] <= bbox->max[1];
}

bool bfBbox2ContainsPoints(BfBbox2 const *bbox, BfPoints2 const *points) {
  for (BfSize i = 0; i < points->size; ++i)
    if (!bfBbox2ContainsPoint(bbox, points->data[i]))
      return false;
  return true;
}

void bfBbox2GetCenter(BfBbox2 const *bbox, BfPoint2 center) {
  center[0] = (bbox->min[0] + bbox->max[0])/2;
  center[1] = (bbox->min[1] + bbox->max[1])/2;
}

void bfBoundingBox3InitEmpty(BfBoundingBox3 *boundingBox) {
  boundingBox->min[0] = boundingBox->min[1] = boundingBox->min[2] = INFINITY;
  boundingBox->max[0] = boundingBox->max[1] = boundingBox->max[2] = -INFINITY;
}

bool bfBoundingBox3IsEmpty(BfBoundingBox3 const *boundingBox) {
  return boundingBox->min[0] >= boundingBox->max[0] ||
         boundingBox->min[1] >= boundingBox->max[1] ||
         boundingBox->min[2] >= boundingBox->max[2];
}

void bfBoundingBox3RescaleToCube(BfBoundingBox3 *boundingBox) {
  BfPoint3 c = {
    (boundingBox->min[0] + boundingBox->max[0])/2,
    (boundingBox->min[0] + boundingBox->max[0])/2,
    (boundingBox->min[0] + boundingBox->max[0])/2
  };

  BfReal dx = boundingBox->max[0] - boundingBox->min[0];
  BfReal dy = boundingBox->max[1] - boundingBox->min[1];
  BfReal dz = boundingBox->max[2] - boundingBox->min[2];

  BfReal dmax = fmax(dx, fmax(dy, dz));

  boundingBox->min[0] = c[0] - dmax/2;
  boundingBox->min[1] = c[1] - dmax/2;
  boundingBox->min[2] = c[2] - dmax/2;

  boundingBox->max[0] = c[0] + dmax/2;
  boundingBox->max[1] = c[1] + dmax/2;
  boundingBox->max[2] = c[2] + dmax/2;
}

void bfBoundingBox2GetCenter(BfBoundingBox3 const *boundingBox, BfPoint2 center) {
  center[0] = (boundingBox->min[0] + boundingBox->max[0])/2;
  center[1] = (boundingBox->min[1] + boundingBox->max[1])/2;
  center[2] = (boundingBox->min[2] + boundingBox->max[2])/2;
}
