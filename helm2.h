#pragma once

#include "def.h"
#include "error.h"
#include "geom.h"

BfComplex bfHelm2GetKernelValue(BfPoint2 const srcPt, BfPoint2 const tgtPt, BfReal K);

BfSize bfHelm2RankEstForTwoCircles(BfCircle2 const *circ1,
                                   BfCircle2 const *circ2,
                                   BfReal k, BfReal C, BfReal eps);

void
bfGetHelm2KernelMatrix(BfMat *kernelMat,
                       BfPoints2 const *srcPts, BfPoints2 const *tgtPts,
                       BfReal K);
