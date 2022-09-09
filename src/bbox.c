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
