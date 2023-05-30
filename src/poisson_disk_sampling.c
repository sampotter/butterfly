#include <bf/poisson_disk_sampling.h>

#include <math.h>

#include <bf/assert.h>
#include <bf/bbox.h>
#include <bf/const.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/rand.h>
#include <bf/size_array.h>

#define OFFSET2_SIZE 21

static int offset2[OFFSET2_SIZE][2] = {
  /* Neighbor offsets: */
            {-1, -2}, { 0, -2}, { 1, -2},
  {-2, -1}, {-1, -1}, { 0, -1}, { 1, -1}, { 2, -1},
  {-2,  0}, {-1,  0},           { 1,  0}, { 2,  0},
  {-2,  1}, {-1,  1}, { 0,  1}, { 1,  1}, { 2,  1},
            {-1,  2}, { 0,  2}, { 1,  2},

  /* No offset: */ { 0,  0}
};

typedef struct {
  BfBbox2 const *bbox;
  BfReal h;
  BfSize nx;
  BfSize ny;
  BfSize n;
  BfReal minDist;
  BfSize k;
  BfSize *cellIndex;
  BfPoints2 *samples;
} Workspace;

static void getCellCoords(BfPoint2 const p, Workspace const *w, BfSize *i, BfSize *j) {
  BfReal iFrac = (p[0] - w->bbox->min[0])/w->h;
  BF_ASSERT(0 <= iFrac && iFrac < w->nx);

  BfReal jFrac = (p[1] - w->bbox->min[1])/w->h;
  BF_ASSERT(0 <= jFrac && jFrac < w->ny);

  *i = floor(iFrac);
  *j = floor(jFrac);

  BF_ASSERT(*i < w->nx);
  BF_ASSERT(*j < w->ny);
}

static BfSize getCellIndex(BfPoint2 const p, Workspace const *w) {
  BfSize i, j;
  getCellCoords(p, w, &i, &j);
  return w->ny*i + j;
}

static bool pointValid(BfPoint2 const p, Workspace const *w) {
  BfSize i0, j0;
  getCellCoords(p, w, &i0, &j0);
  for (BfSize k = 0; k < OFFSET2_SIZE; ++k) {
    int64_t i = i0 + offset2[k][0];
    int64_t j = j0 + offset2[k][1];
    if (0 <= i && i < w->nx && 0 <= j && j < w->ny) {
      BfSize l = w->cellIndex[w->ny*i + j];
      if (l != BF_SIZE_BAD_VALUE) {
        BfPoint2 q;
        bfPoints2Get(w->samples, l, q);
        if (bfPoint2Dist(p, q) < w->minDist)
          return false;
      }
    }
  }
  return true;
}

static bool sampleClosePointUniformly(BfPoint2 x, Workspace const *w, BfPoint2 y) {
  BfReal r = w->minDist;
  BfBbox2 bbox = {.min = {-2*r, -2*r}, .max = {2*r, 2*r}};
  for (BfSize i = 0; i < w->k; ++i) {
    BfPoint2 dy;

  sample:
    /* Sample y uniformly in the box [-2r, 2r] x [-2r, 2r]. */
    bfPoint2SampleUniformlyFromBoundingBox(&bbox, dy);

    /* Resample if dy is outside the [r, 2r] annulus. */
    BfReal R = bfPoint2Magnitude(dy);
    if (!(r <= R && R <= 2*r))
      goto sample;

    /* Resample if y doesn't lie in the domain. */
    y[0] = x[0] + dy[0];
    y[1] = x[1] + dy[1];
    if (!bfBbox2ContainsPoint(w->bbox, y))
      goto sample;

    /* If the point is far enough away from the other current points,
     * accept it. */
    if (pointValid(y, w))
      return true;
  }

  y[0] = y[1] = BF_NAN;
  return false;
}

BfPoints2 *bfPoints2SamplePoissonDisk(BfBbox2 const *bbox, BfReal minDist, BfSize k) {
  BF_ERROR_BEGIN();

  BfReal width = bbox->max[0] - bbox->min[0];
  BfReal height = bbox->max[1] - bbox->min[1];

  Workspace w;
  w.bbox = bbox;
  w.h = minDist/sqrt(2.0);
  w.nx = width/w.h + 1;
  w.ny = height/w.h + 1;
  w.n = w.nx*w.ny;
  w.minDist = minDist;
  w.k = k;

  w.samples = bfPoints2NewEmpty();
  HANDLE_ERROR();

  w.cellIndex = bfMemAlloc(w.n, sizeof(BfSize));
  for (BfSize i = 0; i < w.n; ++i)
    w.cellIndex[i] = BF_SIZE_BAD_VALUE;

  BfPoint2 sample;
  bfPoint2SampleUniformlyFromBoundingBox(bbox, sample);
  bfPoints2Append(w.samples, sample);

  BfSizeArray *activeList = bfSizeArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  w.cellIndex[getCellIndex(sample, &w)] = 0;
  bfSizeArrayAppend(activeList, 0);

  BfSize i, j;
  BfPoint2 x;
  while (!bfSizeArrayIsEmpty(activeList)) {
    i = bfSizeArrayGetRand(activeList);
    bfPoints2Get(w.samples, i, x);
    if (sampleClosePointUniformly(x, &w, sample)) {
      j = bfPoints2GetSize(w.samples);
      bfPoints2Append(w.samples, sample);
      bfSizeArrayAppend(activeList, j);
      w.cellIndex[getCellIndex(sample, &w)] = j;
    } else {
      bfSizeArrayDeleteFirst(activeList, i);
    }
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  bfSizeArrayDeinitAndDealloc(&activeList);

  bfMemFree(w.cellIndex);

  return w.samples;
}
