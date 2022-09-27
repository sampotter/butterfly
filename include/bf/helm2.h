#pragma once

#include "def.h"
#include "error.h"
#include "geom.h"
#include "mat_dense_complex.h"
#include "quadrature.h"
#include "quadtree.h"

BfComplex bfHelm2GetKernelValue(BfReal r, BfReal K);

BfSize bfHelm2RankEstForTwoCircles(BfCircle2 const *circ1,
                                   BfCircle2 const *circ2,
                                   BfReal k, BfReal C, BfReal eps);

BfMatDenseComplex *bfGetHelm2KernelMatrix(BfPoints2 const *srcPts,
                                          BfPoints2 const *tgtPts, BfReal K);

BfMat *bf_hh2_get_dGdN(BfPoints2 const *Xsrc, BfPoints2 const *Xtgt, BfReal K,
                       BfVectors2 const *Ntgt);

BfComplex bf_hh2_get_dGdN_1(BfPoint2 const xsrc, BfPoint2 const xtgt, BfReal K,
                            BfVector2 const ntgt);

BfMatDenseComplex *
bfHelm2GetReexpansionMatrix(BfPoints2 const *srcPtsOrig, BfPoints2 const *srcPtsEquiv,
                      BfPoints2 const *tgtPts, BfReal K);
