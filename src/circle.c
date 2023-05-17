#include <bf/circle.h>

#include <math.h>

#include <bf/const.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/points.h>
#include <bf/vectors.h>

BfPoints2 bfCircle2SamplePoints(BfCircle const *circ, BfSize numPoints) {
  BfPoints2 points;
  bfInitEmptyPoints2(&points, numPoints);

  BfReal const scale = BF_TWO_PI/(BfReal)numPoints;

  BfPoint2 *point = points.data;
  for (BfSize i = 0; i < numPoints; ++i) {
    BfReal theta = scale*i;
    point[i][0] = circ->r*cos(theta) + circ->center[0];
    point[i][1] = circ->r*sin(theta) + circ->center[1];
  }

  return points;
}

BfVectors2 bfCircle2SampleUnitNormals(BfCircle const *circ, BfSize n) {
  BF_ERROR_BEGIN();

  (void)circ; /* Don't actually use this... just for consistency in
               * the interface */

  BfVectors2 unitNormals;

  bfInitEmptyVectors2(&unitNormals, n);
  HANDLE_ERROR();

  BfReal const scale = BF_TWO_PI/n;

  BfVector2 *vector = unitNormals.data;
  for (BfSize i = 0; i < n; ++i) {
    BfReal theta = scale*i;
    vector[i][0] = cos(theta);
    vector[i][1] = sin(theta);
  }

  BF_ERROR_END()
    bfFreeVectors2(&unitNormals);

  return unitNormals;
}

bool bfCircle2ContainsPoint(BfCircle const *circ, BfPoint2 const point) {
  return bfPoint2Dist(circ->center, point) <= circ->r;
}

bool bfCircle2ContainsPoints(BfCircle const *circ, BfPoints2 const *points) {
  for (BfSize i = 0; i < points->size; ++i)
    if (!bfCircle2ContainsPoint(circ, points->data[i]))
      return false;
  return true;
}
