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

BfReal
bfHelm2RankEstForTwoCircles(BfCircle2 const circ1, BfCircle2 const circ2,
                            BfReal k, BfReal C, BfReal eps)
{
  BfReal R = bfPoint2Dist(circ1.center, circ2.center);

  return BF_TWO_PI*k*circ1.r*circ2.r/(R - circ1.r - circ2.r) - C*log10(eps);
}

enum BfError
bfHelm2KernelMatrixFromPoints(BfMat *K, BfMat const *X, BfMat const *Y, BfReal k)
{
  BfSize m = K->shape[0], n = K->shape[1];

  BfReal *x, *y;
  BfComplex *row;

  for (BfSize i = 0; i < m; ++i) {
    x = ((BfReal *)X->data) + 2*i;
    row = ((BfComplex *)K->data) + n*i;
    for (BfSize j = 0; j < n; ++j) {
      y = ((BfReal *)Y->data) + 2*j;
      row[j] = bfHelm2GetKernelValue(x, y, k);
    }
  }

  return BF_ERROR_NO_ERROR;
}
