#include "helm2.h"

#include <math.h>

#include "const.h"
#include "geom.h"

enum BfError
bfHelm2RankEstForTwoCircles(BfCircle2 const circ1, BfCircle2 const circ2,
                            BfReal k, BfReal C, BfReal eps, BfSize *rank)
{
  if (k <= 0 || C <= 0 || eps <= 0)
    return BF_ERROR_INVALID_ARGUMENTS;

  BfReal R = bfPoint2Dist(circ1.center, circ2.center);
  if (R <= 0)
    return BF_ERROR_INVALID_ARGUMENTS;

  /* Reading Michielssen & Boag, seems like there should be a factor
   * of two pi in front of this, but that leads to a rank estimate
   * which is way higher than necessary. */
  BfReal p = k*circ1.r*circ2.r/(R - circ1.r - circ2.r) - C*log10(eps);
  if (p <= 0)
    return BF_ERROR_INVALID_ARGUMENTS;

  *rank = ceil(p);

  return BF_ERROR_NO_ERROR;
}

BfComplex
bfHelm2GetKernelValue(BfPoint2 const srcPt, BfPoint2 const tgtPt, BfReal K)
{
  BfReal r = bfPoint2Dist(srcPt, tgtPt);
  BfReal arg = K*r;
  return (I*j0(arg) - y0(arg))/4;
}

enum BfError
bfGetHelm2KernelMatrix(BfMat *kernelMat,
                       BfPoints2 const *srcPts, BfPoints2 const *tgtPts,
                       BfReal K)
{
  if (bfMatIsInitialized(kernelMat))
    return BF_ERROR_INVALID_ARGUMENTS;

  BfSize m = tgtPts->size; /* number of rows */
  BfSize n = srcPts->size; /* number of columns */

  bfInitEmptyMat(kernelMat, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE, m, n);

  BfPoint2 *srcPt = srcPts->data;
  BfPoint2 *tgtPt = tgtPts->data;

  BfVec row;
  for (BfSize i = 0; i < m; ++i) {
    row = bfGetMatRow(kernelMat, i);
    for (BfSize j = 0; j < n; ++j)
      *(BfComplex *)bfVecGetEltPtr(&row, j) =
        bfHelm2GetKernelValue(srcPt[j], tgtPt[i], K);
  }

  return BF_ERROR_NO_ERROR;
}
