#include <bf/circle.h>

#include <math.h>

#include <bf/assert.h>
#include <bf/const.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/points.h>
#include <bf/vectors.h>

BfPoints2 *bfCircle2SamplePoints(BfCircle const *circ, BfSize numPoints) {
  BF_ERROR_BEGIN();

  BfPoints2 *points = bfPoints2NewWithCapacity(numPoints);
  HANDLE_ERROR();

  BfReal const scale = BF_TWO_PI/(BfReal)numPoints;

  BfPoint2 point;
  for (BfSize i = 0; i < numPoints; ++i) {
    BfReal theta = scale*i;
    point[0] = circ->r*cos(theta) + circ->center[0];
    point[1] = circ->r*sin(theta) + circ->center[1];
    bfPoints2Append(points, point);
    HANDLE_ERROR();
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return points;
}

BfVectors2 *bfCircle2SampleUnitNormals(BfCircle const *circ, BfSize n) {
  BF_ERROR_BEGIN();

  (void)circ; /* Don't actually use this... just for consistency in
               * the interface */

  BfVectors2 *unitNormals = bfVectors2NewWithCapacity(n);
  HANDLE_ERROR();

  BfReal const scale = BF_TWO_PI/n;

  for (BfSize i = 0; i < n; ++i) {
    BfReal theta = scale*i;
    BfVector2 unitNormal = {cos(theta), sin(theta)};
    bfVectors2Append(unitNormals, unitNormal);
  }

  BF_ERROR_END() {
    BF_DIE();
  }

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
