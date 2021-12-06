#include "helm2.h"

#include <math.h>

#include "const.h"
#include "geom.h"

enum BfError
bfHelm2RankEstForTwoCircles(BfCircle2 const circ1, BfCircle2 const circ2,
                            BfReal k, BfReal C, BfReal eps, BfReal *p)
{
  if (k <= 0 || C <= 0 || eps <= 0)
    return BF_ERROR_INVALID_ARGUMENTS;

  BfReal R = bfPoint2Dist(circ1.center, circ2.center);
  if (R <= 0)
    return BF_ERROR_INVALID_ARGUMENTS;

  /* Reading Michielssen & Boag, seems like there should be a factor
   * of two pi in front of this, but that leads to a rank estimate
   * which is way higher than necessary. */
  *p = k*circ1.r*circ2.r/(R - circ1.r - circ2.r) - C*log10(eps);
  if (*p <= 0)
    return BF_ERROR_INVALID_ARGUMENTS;

  return BF_ERROR_NO_ERROR;
}

BfComplex
bfHelm2GetKernelValue(BfVec const *x, BfVec const *y, BfReal k)
{
  BfReal r = bfVecDist(x, y);
  BfReal arg = k*r;
  return (I*j0(arg) - y0(arg))/4;
}

enum BfError
bfHelm2KernelMatrixFromPoints(BfMat *K, BfMat const *X, BfMat const *Y, BfReal k)
{
  if (K->dtype != BF_DTYPE_COMPLEX)
    return BF_ERROR_BAD_DTYPE;

  if (X->dtype != BF_DTYPE_REAL)
    return BF_ERROR_BAD_DTYPE;

  if (Y->dtype != BF_DTYPE_REAL)
    return BF_ERROR_BAD_DTYPE;

  BfSize m = bfMatNumRows(K), n = bfMatNumCols(K);

  if (bfMatNumRows(X) != n)
    return BF_ERROR_INVALID_ARGUMENTS;

  if (bfMatNumCols(X) != 2)
    return BF_ERROR_INVALID_ARGUMENTS;

  if (bfMatNumRows(Y) != m)
    return BF_ERROR_INVALID_ARGUMENTS;

  if (bfMatNumCols(Y) != 2)
    return BF_ERROR_INVALID_ARGUMENTS;

  BfVec x, y, row;
  for (BfSize i = 0; i < m; ++i) {
    row = bfGetMatRow(K, i);
    y = bfGetMatRow(Y, i);
    for (BfSize j = 0; j < n; ++j) {
      x = bfGetMatRow(X, j);
      *(BfComplex *)bfVecGetEltPtr(&row, j) = bfHelm2GetKernelValue(&x, &y, k);
    }
  }

  return BF_ERROR_NO_ERROR;
}
