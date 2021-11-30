#include "helm2.h"

#include <math.h>

#include "const.h"
#include "geom.h"

BfComplex
bfHelm2GetKernelValue(BfPoint2 const p, BfPoint2 const q, BfReal k)
{
  BfReal arg = k*bfPoint2Dist(p, q);
  return (I*j0(arg) - y0(arg))/4;
}

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

enum BfError
bfHelm2KernelMatrixFromPoints(BfMat *K, BfMat const *X, BfMat const *Y, BfReal k)
{
  BfSize m = K->shape[0], n = K->shape[1];

  if (X->shape[0] != n)
    return BF_ERROR_INVALID_ARGUMENTS;

  if (X->shape[1] != 2)
    return BF_ERROR_INVALID_ARGUMENTS;

  if (Y->shape[0] != m)
    return BF_ERROR_INVALID_ARGUMENTS;

  if (Y->shape[1] != 2)
    return BF_ERROR_INVALID_ARGUMENTS;

  BfReal *x, *y;
  BfComplex *row;

  enum BfError error = BF_ERROR_NO_ERROR;

  for (BfSize i = 0; i < m; ++i)
  {
    error = bfGetMatRow(K, i, (void **)&row);
    if (error)
      return error;

    error = bfGetMatRow(Y, i, (void **)&y);
    if (error)
      return error;

    for (BfSize j = 0; j < n; ++j)
    {
      error = bfGetMatRow(X, j, (void **)&x);
      if (error)
        return error;

      row[j] = bfHelm2GetKernelValue(x, y, k);
    }
  }

  return BF_ERROR_NO_ERROR;
}
