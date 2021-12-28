#include "geom.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "const.h"

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

BfReal bfPoint2Dist(BfPoint2 const p, BfPoint2 const q) {
  return hypot(q[0] - p[0], q[1] - p[1]);
}

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

BfBbox2 bfGetEmptyBbox2() {
  return (BfBbox2) {
    .min = { INFINITY,  INFINITY},
    .max = {-INFINITY, -INFINITY}
  };
}

bool bfBbox2IsEmpty(BfBbox2 const *bbox) {
  return bbox->min[0] >= bbox->max[0] && bbox->min[1] >= bbox->max[1];
}

BfPoints2 bfGetUninitializedPoints2() {
  return (BfPoints2) {.data = NULL, .size = 0};
}

enum BfError bfInitEmptyPoints2(BfPoints2 *points, BfSize numPoints) {
  if (numPoints == 0)
    return BF_ERROR_INVALID_ARGUMENTS;

  points->size = numPoints;

  points->data = malloc(numPoints*sizeof(BfPoint2));
  if (points->data == NULL)
    return BF_ERROR_MEMORY_ERROR;

  return BF_ERROR_NO_ERROR;
}

enum BfError bfReadPoints2FromFile(char const *path, BfPoints2 *points) {
  enum BfError error = BF_ERROR_NO_ERROR;

  /* open the file for reading */
  FILE *fp = fopen(path, "r");
  if (fp == NULL) {
    error = BF_ERROR_RUNTIME_ERROR;
    goto cleanup;
  }

  /* get the size of the file */
  fseek(fp, 0, SEEK_END);
  BfSize size = ftell(fp);
  fseek(fp, 0, SEEK_SET);

  if (size % sizeof(BfPoint2) != 0) {
    error = BF_ERROR_RUNTIME_ERROR;
    goto cleanup;
  }

  /* get the number of points */
  points->size = size/sizeof(BfPoint2);

  /* allocate space for the points and read them in */
  points->data = malloc(size);
  fread(points->data, sizeof(BfPoint2), points->size, fp);

cleanup:
  fclose(fp);
  if (error)
    free(points->data);

  return error;
}

void bfFreePoints2(BfPoints2 *points) {
  free(points->data);
}

bool bfPoints2Initialized(BfPoints2 const *points) {
  return points->data != NULL && points->size != 0;
}

BfBbox2 bfGetPoints2BoundingBox(BfPoints2 const *points) {
  BfBbox2 bbox = bfGetEmptyBbox2();

  BfSize numPoints = points->size;

  for (size_t i = 0; i < numPoints; ++i) {
    BfReal const *point = points->data[i];

    BfReal x = point[0], y = point[1];

    bbox.min[0] = fmin(bbox.min[0], x);
    bbox.max[0] = fmax(bbox.max[0], x);
    bbox.min[1] = fmin(bbox.min[1], y);
    bbox.max[1] = fmax(bbox.max[1], y);
  }

  return bbox;
}

enum BfError
bfGetPointsByIndex(BfPoints2 const *points,
                   BfSize numInds, BfSize const *inds,
                   BfPoints2 *indexedPoints)
{
  enum BfError error = bfInitEmptyPoints2(indexedPoints, numInds);
  if (error)
    return error;

  BfPoint2 const *point = (BfPoint2 const *)points->data;
  BfPoint2 *indexedPoint = indexedPoints->data;

  for (BfSize i = 0, j; i < numInds; ++i) {
    j = inds[i];
    indexedPoint[i][0] = point[j][0];
    indexedPoint[i][1] = point[j][1];
  }

  return BF_ERROR_NO_ERROR;
}

void bfPrintPoints2(BfPoints2 const *points) {
  for (BfSize i = 0; i < points->size; ++i)
    printf("points[%lu] = (%g, %g)\n",i,points->data[i][0],points->data[i][1]);
}

enum BfError bfSavePoints2(BfPoints2 const *points, char const *path) {
  enum BfError error = BF_ERROR_NO_ERROR;

  FILE *fp = fopen(path, "w");
  if (fp == NULL) {
    error = BF_ERROR_FILE_ERROR;
    goto cleanup;
  }

  fwrite(points->data, points->size, sizeof(BfPoint2), fp);

cleanup:
  fclose(fp);

  return error;
}
