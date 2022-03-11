#pragma once

#include "def.h"
#include "error.h"
#include "geom.h"
#include "mat_dense_complex.h"

BfComplex bfHelm2GetKernelValue(BfReal r, BfReal K);

BfSize bfHelm2RankEstForTwoCircles(BfCircle2 const *circ1,
                                   BfCircle2 const *circ2,
                                   BfReal k, BfReal C, BfReal eps);

BfMatDenseComplex *
bfGetHelm2KernelMatrix(BfPoints2 const *srcPts, BfPoints2 const *tgtPts,
                       BfReal K);

BfMatDenseComplex *
bfHelm2GetShiftMatrix(BfPoints2 const *srcPtsOrig, BfPoints2 const *srcPtsEquiv,
                      BfPoints2 const *tgtPts, BfReal K);
