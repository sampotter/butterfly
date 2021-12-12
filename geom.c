#include "geom.h"

#include <math.h>

#include "const.h"

bool bfBbox2ContainsPoint(BfBbox2 const *bbox, BfPoint2 const point) {
  return bbox->min[0] <= point[0] && point[0] <= bbox->max[0]
      && bbox->min[1] <= point[1] && point[1] <= bbox->max[1];
}

BfReal bfPoint2Dist(BfPoint2 const p, BfPoint2 const q) {
  return hypot(q[0] - p[0], q[1] - p[1]);
}

enum BfError
bfSamplePointsOnCircle2(BfCircle2 const *circ, BfMat *X)
{
  if (X->numCols != 2)
    return BF_ERROR_INVALID_ARGUMENTS;

  BfSize num_points = X->numRows;

  BfReal scale = BF_TWO_PI/((BfReal)num_points), *x;
  for (BfSize i = 0; i < num_points; ++i) {
    BfReal theta = scale*i;
    x = ((BfReal *)X->data) + 2*i;
    x[0] = circ->r*cos(theta) + circ->center[0];
    x[1] = circ->r*sin(theta) + circ->center[1];
  }

  return BF_ERROR_NO_ERROR;
}
