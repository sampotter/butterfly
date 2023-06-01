#include <bf/points.h>

#include <math.h>

#include <bf/assert.h>
#include <bf/bbox.h>
#include <bf/const.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_dense_real.h>
#include <bf/mem.h>
#include <bf/rand.h>
#include <bf/vec_real.h>

BfReal bfPoint2Dist(BfPoint2 const p, BfPoint2 const q) {
  return hypot(q[0] - p[0], q[1] - p[1]);
}

BfReal bfPoint2Magnitude(BfPoint2 const p) {
  return hypot(p[0], p[1]);
}

void bfPoint2SampleUniformlyFromBoundingBox(BfBbox2 const *bbox, BfPoint2 p) {
  p[0] = (bbox->max[0] - bbox->min[0])*bfRealUniform1() + bbox->min[0];
  p[1] = (bbox->max[1] - bbox->min[1])*bfRealUniform1() + bbox->min[1];
}

void bfPoint2RotateAboutOrigin(BfPoint2 p, BfReal theta) {
  BfReal c = cos(theta);
  BfReal s = sin(theta);
  BfReal tmp = c*p[0] - s*p[1];
  p[1] = s*p[0] + c*p[1];
  p[0] = tmp;
}

void bfPoint2Translate(BfPoint2 p, BfVector2 const u) {
  p[0] += u[0];
  p[1] += u[1];
}

BfReal bfPoint3Dist(BfPoint3 const p, BfPoint3 const q) {
  BfVector3 u = {q[0] - p[0], q[1] - p[1], q[2] - p[2]};
  return sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
}

void bfPoint3Sub(BfPoint3 const u, BfPoint3 const v, BfVector3 w) {
  w[0] = u[0] - v[0];
  w[1] = u[1] - v[1];
  w[2] = u[2] - v[2];
}

void bfPoint3GetPointOnRay(BfPoint3 const r0, BfVector3 const dr, BfReal t, BfPoint3 rt) {
  rt[0] = r0[0] + t*dr[0];
  rt[1] = r0[1] + t*dr[1];
  rt[2] = r0[2] + t*dr[2];
}

void bfPoint3Copy(BfPoint3 x, BfPoint3 const y) {
  x[0] = y[0];
  x[1] = y[1];
  x[2] = y[2];
}

BfPoints2 *bfPoints2NewEmpty() {
  BF_ERROR_BEGIN();

  BfPoints2 *points = bfMemAlloc(1, sizeof(BfPoints2));
  if (points == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  points->size = 0;
  points->capacity = 16;
  points->isView = BF_ARRAY_DEFAULT_CAPACITY;

  points->data = bfMemAlloc(points->capacity, sizeof(BfPoint2));

  BF_ERROR_END() {
    BF_DIE();
  }

  return points;
}

BfPoints2 *bfPoints2NewGrid(BfBbox2 const *bbox, BfSize nx, BfSize ny) {
  BF_ERROR_BEGIN();

  BfPoints2 *points = bfPoints2NewEmpty();
  HANDLE_ERROR();

  BfReal hx = (bbox->max[0] - bbox->min[0])/(nx - 1);
  BfReal hy = (bbox->max[1] - bbox->min[1])/(ny - 1);

  for (BfSize i = 0; i < nx; ++i) {
    for (BfSize j = 0; j < ny; ++j) {
      BfPoint2 p = {hx*i + bbox->min[0], hy*j + bbox->min[1]};
      bfPoints2Append(points, p);
      HANDLE_ERROR();
    }
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return points;
}

BfPoints2 const *bfPoints2ConstViewFromMat(BfMat const *mat) {
  return bfPoints2ConstViewFromMatDenseReal(bfMatConstToMatDenseRealConst(mat));
}

BfPoints2 const *bfPoints2ConstViewFromMatDenseReal(BfMatDenseReal const *matDenseReal) {
  BF_ERROR_BEGIN();

  BfPoints2 *points = NULL;

  BfMat const *mat = bfMatDenseRealConstToMatConst(matDenseReal);
  BfMatDense const *matDense = bfMatDenseRealConstToMatDenseConst(matDenseReal);

  if (bfMatDenseRealGetNumCols(mat) != 2)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (bfMatDenseGetColStride(matDense) != 1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  points = bfMemAlloc(1, sizeof(BfPoints2));
  points->data = (BfPoint2 *)matDenseReal->data;
  points->size = bfMatDenseRealGetNumRows(mat);
  points->capacity = BF_SIZE_BAD_VALUE;
  points->isView = true;

  BF_ERROR_END() {
    bfMemFree(points);
    points = NULL;
  }

  return points;
}

/** Implementation: Points1 */

BfPoints1 *bfPoints1New() {
  BF_ERROR_BEGIN();

  BfPoints1 *points = bfMemAlloc(1, sizeof(BfPoints1));
  if (points == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  points->data = NULL;
  points->size = BF_SIZE_BAD_VALUE;
  points->capacity = BF_SIZE_BAD_VALUE;

  BF_ERROR_END() {
    points = NULL;
  }

  return points;
}

BfPoints1 *bfPoints1Copy(BfPoints1 const *points) {
  BF_ERROR_BEGIN();

  BfPoints1 *pointsCopy = bfPoints1New();
  HANDLE_ERROR();

  bfPoints1InitEmpty(pointsCopy, points->size);
  HANDLE_ERROR();

  bfMemCopy(points->data, points->size, sizeof(BfReal), pointsCopy->data);

  pointsCopy->size = points->size;

  BF_ERROR_END() {
    BF_DIE();
  }

  return pointsCopy;
}

void bfPoints1InitEmpty(BfPoints1 *points, BfSize capacity) {
  BF_ERROR_BEGIN();

  points->data = bfMemAlloc(capacity, sizeof(BfPoint1));
  if (points->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  points->size = 0;
  points->capacity = capacity;
  points->isView = false;

  BF_ERROR_END() {
    bfPoints1Deinit(points);
  }
}

void bfPoints1InitViewFromVecReal(BfPoints1 *points, BfVecReal const *vecReal) {
  BfVec const *vec = bfVecRealConstToVecConst(vecReal);

  points->data = vecReal->data;
  points->size = vec->size;
  points->capacity = BF_SIZE_BAD_VALUE;
  points->isView = true;
}

void bfPoints1Deinit(BfPoints1 *points) {
  if (!points->isView)
    bfMemFree(points->data);
  points->data = NULL;
  points->size = BF_SIZE_BAD_VALUE;
  points->capacity = BF_SIZE_BAD_VALUE;
}

void bfPoints1Dealloc(BfPoints1 **points) {
  bfMemFree(*points);
  *points = NULL;
}

void bfPoints1Delete(BfPoints1 **points) {
  bfPoints1Deinit(*points);
  bfPoints1Dealloc(points);
}

bool bfPoints1IsSorted(BfPoints1 const *points) {
  for (BfSize i = 1; i < points->size; ++i)
    if (points->data[i - 1] > points->data[i])
      return false;
  return true;
}

void bfPoints1Append(BfPoints1 *points, BfPoint1 point) {
  BF_ERROR_BEGIN();

  /* Grow the array if we're at capacity */
  if (points->size == points->capacity) {
    points->capacity *= 2;

    BfReal *newData = bfMemRealloc(points->data, points->capacity, sizeof(BfReal));
    if (newData == NULL)
      RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

    points->data = newData;
  }

  /* Append new point */
  points->data[points->size++] = point;

  BF_ERROR_END() {}
}

void bfPoints1InsertPointsSorted(BfPoints1 *points, BfPoints1 const *newPoints) {
  BF_ERROR_BEGIN();

  if (!bfPoints1IsSorted(points))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (!bfPoints1IsSorted(newPoints))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize newSize = points->size + newPoints->size;

  BfReal *newData = bfMemAlloc(newSize, sizeof(BfReal));
  if (newData == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfSize i = 0, j = 0, k = 0;
  while ((i < points->size || j < newPoints->size) && k < newSize) {
    if (i < points->size && j < newPoints->size) {
      if (points->data[i] < newPoints->data[j]) {
        newData[k++] = points->data[i++];
      } else {
        newData[k++] = newPoints->data[j++];
      }
    } else if (i < points->size) {
      newData[k++] = points->data[i++];
    } else if (j < newPoints->size) {
      newData[k++] = newPoints->data[j++];
    } else {
      BF_DIE();
    }
  }
  BF_ASSERT(i == points->size && j == newPoints->size && k == newSize);

  if (!points->isView)
    bfMemFree(points->data);
  points->data = newData;
  points->size = newSize;
  points->capacity = newSize;

  BF_ERROR_END() {
    bfMemFree(newData);
  }
}

void bfPoints1Map(BfPoints1 *points, BfReal (*func)(BfReal)) {
  BF_ERROR_BEGIN();

  for (BfSize i = 0; i < points->size; ++i) {
    points->data[i] = func(points->data[i]);
    HANDLE_ERROR();
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfPoints1Save(BfPoints1 const *points, char const *path) {
  BF_ERROR_BEGIN();

  FILE *fp = fopen(path, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  fwrite(points->data, points->size, sizeof(BfPoint1), fp);
  // TODO: error-handling

  BF_ERROR_END() {}

  fclose(fp);
}

/** Implementation: Points2 */

BfPoints2 bfGetUninitializedPoints2() {
  return (BfPoints2) {.data = NULL, .size = 0};
}

void bfInitEmptyPoints2(BfPoints2 *points, BfSize numPoints) {
  if (numPoints == 0) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  points->size = numPoints;

  points->data = bfMemAlloc(numPoints, sizeof(BfPoint2));
  if (points->data == NULL)
    bfSetError(BF_ERROR_MEMORY_ERROR);

#if BF_DEBUG
  for (BfSize i = 0; i < numPoints; ++i)
    points->data[i][0] = points->data[i][1] = BF_NAN;
#endif
}

void bfReadPoints2FromFile(char const *path, BfPoints2 *points) {
  BF_ERROR_BEGIN();

  /* open the file for reading */
  FILE *fp = fopen(path, "r");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* get the size of the file */
  fseek(fp, 0, SEEK_END);
  BfSize size = ftell(fp);
  fseek(fp, 0, SEEK_SET);

  /* make sure the binary file is the right size */
  if (size % sizeof(BfPoint2) != 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* get the number of points */
  points->size = size/sizeof(BfPoint2);

  /* allocate space for the points */
  points->data = bfMemAlloc(size, sizeof(char));
  if (points->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* read them in */
  fread(points->data, sizeof(BfPoint2), points->size, fp);
  // TODO: error-handling

  BF_ERROR_END() {
    bfMemFree(points->data);
  }

  fclose(fp);
}

void bfFreePoints2(BfPoints2 *points) {
  if (!points->isView)
    bfMemFree(points->data);
}

void bfPoints2Dealloc(BfPoints2 **points) {
  bfMemFree(*points);
  *points = NULL;
}

void bfPoints2DeinitAndDealloc(BfPoints2 **points) {
  bfFreePoints2(*points);
  bfPoints2Dealloc(points);
}

bool bfPoints2Initialized(BfPoints2 const *points) {
  return points->data != NULL && points->size != 0;
}

BfBbox2 bfPoints2GetBoundingBox(BfPoints2 const *points) {
  BfBbox2 bbox = bfGetEmptyBbox2();

  BfSize numPoints = points->size;

  for (BfSize i = 0; i < numPoints; ++i) {
    BfReal const *point = points->data[i];

    BfReal x = point[0], y = point[1];

    bbox.min[0] = fmin(bbox.min[0], x);
    bbox.max[0] = fmax(bbox.max[0], x);
    bbox.min[1] = fmin(bbox.min[1], y);
    bbox.max[1] = fmax(bbox.max[1], y);
  }

  return bbox;
}

void bfGetPointsByIndex(BfPoints2 const *points,
                        BfSize numInds, BfSize const *inds,
                        BfPoints2 *indexedPoints)
{
  BF_ERROR_BEGIN();

  bfInitEmptyPoints2(indexedPoints, numInds);
  HANDLE_ERROR();

  BfPoint2 const *point = (BfPoint2 const *)points->data;
  BfPoint2 *indexedPoint = indexedPoints->data;

  for (BfSize i = 0, j; i < numInds; ++i) {
    j = inds[i];
    indexedPoint[i][0] = point[j][0];
    indexedPoint[i][1] = point[j][1];
  }

  BF_ERROR_END() {
    bfFreePoints2(indexedPoints);
  }
}

void bfPrintPoints2(BfPoints2 const *points) {
  for (BfSize i = 0; i < points->size; ++i)
    printf("points[%lu] = (%g, %g)\n",i,points->data[i][0],points->data[i][1]);
}

void bfSavePoints2(BfPoints2 const *points, char const *path) {
  BF_ERROR_BEGIN();

  FILE *fp = fopen(path, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  fwrite(points->data, points->size, sizeof(BfPoint2), fp);
  // TODO: error-handling

  BF_ERROR_END() {}

  fclose(fp);
}

/* Compute all pairwise distances between the point sets `X` and
 * `Y`. If `|X| == m` and `|Y| == n`, then the result is a length
 * `m*n` array of `BfReal`s where the `m*i + j` element is the
 * distance between `X[i]` and `Y[j]`.
 *
 * It is the responsibility of the caller to free the returned
 * array. */
BfReal *bfPoints2PairwiseDists(BfPoints2 const *X, BfPoints2 const *Y) {
  BF_ERROR_BEGIN();

  BfSize m = X->size, n = Y->size;

  BfReal *r = bfMemAlloc(m*n, sizeof(BfReal));
  if (r == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfSize k = 0;
  for (BfSize i = 0; i < m; ++i) {
    BfReal const *x = &X->data[i][0];
    for (BfSize j = 0; j < n; ++j) {
      BfReal const *y = &Y->data[j][0];
      r[k++] = hypot(y[0] - x[0], y[1] - x[1]);
    }
  }

  BF_ERROR_END() {}

  return r;
}

void bfPoints2Append(BfPoints2 *points, BfPoint2 const p) {
  BF_ERROR_BEGIN();

  /* Grow the array if we're at capacity */
  if (points->size == points->capacity) {
    points->capacity *= 2;

    BfPoint2 *newData = bfMemRealloc(points->data, points->capacity, sizeof(BfPoint2));
    if (newData == NULL)
      RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

    points->data = newData;
  }

  /* Append new point */
  bfMemCopy(p, 1, sizeof(BfPoint2), &points->data[points->size++]);

  BF_ERROR_END() {}
}

void bfPoints2Extend(BfPoints2 *points, BfPoints2 const *newPoints) {
  /* TODO: bad implementation to start... fix */
  for (BfSize i = 0; i < newPoints->size; ++i) {
    BfPoint2 p;
    bfPoints2Get(newPoints, i, p);
    bfPoints2Append(points, p);
  }
}

void bfPoints2Get(BfPoints2 const *points, BfSize i, BfPoint2 p) {
  BF_ERROR_BEGIN();

  if (i >= points->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  bfMemCopy(&points->data[i], 1, sizeof(BfPoint2), p);

  BF_ERROR_END() {
    BF_DIE();
  }
}

BfSize bfPoints2GetSize(BfPoints2 const *points) {
  return points->size;
}

BfPoints2 *bfPoints2GetRangeView(BfPoints2 *points, BfSize i0, BfSize i1) {
  BF_ERROR_BEGIN();

  if ((i0 >= points->size && i1 > points->size) || i0 > points->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (i0 > i1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfPoints2 *pointsView = bfMemAlloc(1, sizeof(BfPoints2));
  if (pointsView == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  pointsView->size = i1 - i0;
  pointsView->capacity = BF_SIZE_BAD_VALUE;
  pointsView->isView = true;
  pointsView->data = &points->data[i0];

  BF_ERROR_END() {
    BF_DIE();
  }

  return pointsView;
}

/** Implementation: Points3 */

void bfPoints3InitEmpty(BfPoints3 *points, BfSize numPoints) {
  if (numPoints == 0) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  points->size = numPoints;

  points->data = bfMemAlloc(numPoints, sizeof(BfPoint3));
  if (points->data == NULL)
    bfSetError(BF_ERROR_MEMORY_ERROR);
}

void bfPoints3InitFromBinaryFile(BfPoints3 *points, char const *path) {
  BF_ERROR_BEGIN();

  /* open the file for reading */
  FILE *fp = fopen(path, "r");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* get the size of the file */
  fseek(fp, 0, SEEK_END);
  BfSize size = ftell(fp);
  fseek(fp, 0, SEEK_SET);

  /* make sure the binary file is the right size */
  if (size % sizeof(BfPoint3) != 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* get the number of points */
  points->size = size/sizeof(BfPoint3);

  /* allocate space for the points */
  points->data = bfMemAlloc(size, sizeof(char));
  if (points->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* read them in */
  fread(points->data, sizeof(BfPoint3), points->size, fp);
  // TODO: error-handling

  BF_ERROR_END() {
    bfMemFree(points->data);
  }

  fclose(fp);
}

void bfPoints3Deinit(BfPoints3 *points) {
  bfMemFree(points->data);
  points->data = NULL;

  points->size = BF_SIZE_BAD_VALUE;
}

BfBoundingBox3 bfPoints3GetBoundingBox(BfPoints3 const *points) {
  BfBoundingBox3 boundingBox;
  bfBoundingBox3InitEmpty(&boundingBox);

  BfSize numPoints = points->size;

  for (BfSize i = 0; i < numPoints; ++i) {
    BfReal const *point = points->data[i];
    for (BfSize j = 0; j < 3; ++j) {
      boundingBox.min[j] = fmin(boundingBox.min[j], point[j]);
      boundingBox.max[j] = fmax(boundingBox.max[j], point[j]);
    }
  }

  return boundingBox;
}

void bfPoints3GetByIndex(BfPoints3 const *points, BfSize numInds, BfSize const *inds, BfPoints3 *indexedPoints) {
  BF_ERROR_BEGIN();

  bfPoints3InitEmpty(indexedPoints, numInds);
  HANDLE_ERROR();

  for (BfSize i = 0; i < numInds; ++i)
    bfPoint3Copy(indexedPoints->data[i], points->data[inds[i]]);

  BF_ERROR_END() {
    bfPoints3Deinit(indexedPoints);
  }
}
