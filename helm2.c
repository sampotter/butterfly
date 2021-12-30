#include "helm2.h"

#include <assert.h>
#include <math.h>

#include "const.h"
#include "error_macros.h"
#include "geom.h"

BfSize bfHelm2RankEstForTwoCircles(BfCircle2 const *circ1,
                                   BfCircle2 const *circ2,
                                   BfReal k, BfReal C, BfReal eps)
{
  assert(k > 0);
  assert(C > 0);
  assert(eps > 0);

  BfReal r1 = circ1->r;
  BfReal r2 = circ2->r;
  BfReal R = bfPoint2Dist(circ1->center, circ2->center);
  assert(r1 > 0 && r2 > 0 && R > 0);

  /* Reading Michielssen & Boag, seems like there should be a factor
   * of two pi in front of this, but that leads to a rank estimate
   * which is way higher than necessary. */
  BfReal p = k*r1*r2/(R - r1 - r2) - C*log10(eps);
  assert(p > 0);

  BfSize rank = (BfSize)ceil(p);
  assert(rank > 0);
  return rank;
}

BfComplex
bfHelm2GetKernelValue(BfPoint2 const srcPt, BfPoint2 const tgtPt, BfReal K)
{
  BfReal r = bfPoint2Dist(srcPt, tgtPt);
  BfReal arg = K*r;
  return (I*j0(arg) - y0(arg))/4;
}

void
bfGetHelm2KernelMatrix(BfMat *kernelMat,
                       BfPoints2 const *srcPts, BfPoints2 const *tgtPts,
                       BfReal K)
{
  enum BfError error;
  bool erred = false;

  if (bfMatIsInitialized(kernelMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize m = tgtPts->size; /* number of rows */
  BfSize n = srcPts->size; /* number of columns */

  bfInitEmptyMat(kernelMat, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE, m, n);
  HANDLE_ERROR();

  BfPoint2 *srcPt = srcPts->data;
  BfPoint2 *tgtPt = tgtPts->data;

  BfVec row;
  for (BfSize i = 0; i < m; ++i) {
    row = bfGetMatRow(kernelMat, i);
    for (BfSize j = 0; j < n; ++j)
      *(BfComplex *)bfVecGetEltPtr(&row, j) =
        bfHelm2GetKernelValue(srcPt[j], tgtPt[i], K);
  }

cleanup:
  if (erred)
    bfFreeMat(kernelMat);
}
