#include "helm2.h"

#include <math.h>

#include "const.h"
#include "geom.h"

BfComplex bfHelm2GetKernelValue(BfPoint2 const p, BfPoint2 const q)
{
  BfReal arg = BF_TWO_PI*bfPoint2Dist(p, q);
  return (I*j0(arg) - y0(arg))/4;
}
