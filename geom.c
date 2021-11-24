#include "geom.h"

#include <math.h>

bool bfBbox2ContainsPoint(BfBbox2 const *bbox, BfPoint2 const point) {
  return bbox->min[0] <= point[0] && point[0] <= bbox->max[0]
      && bbox->min[1] <= point[1] && point[1] <= bbox->max[1];
}

BfReal bfPoint2Dist(BfPoint2 const p, BfPoint2 const q) {
  return hypot(q[0] - p[0], q[1] - p[1]);
}
