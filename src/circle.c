#include <bf/circle.h>

#include <math.h>

#include <bf/const.h>
#include <bf/points.h>

BfPoints2 bfSamplePointsOnCircle2(BfCircle2 const *circ, BfSize numPoints) {
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

bool bfCircle2ContainsPoint(BfCircle2 const *circ, BfPoint2 const point) {
  return bfPoint2Dist(circ->center, point) <= circ->r;
}

bool bfCircle2ContainsPoints(BfCircle2 const *circ, BfPoints2 const *points) {
  for (BfSize i = 0; i < points->size; ++i)
    if (!bfCircle2ContainsPoint(circ, points->data[i]))
      return false;
  return true;
}
