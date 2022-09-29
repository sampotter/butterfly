#pragma once

#include "def.h"
#include "error.h"
#include "geom.h"
#include "mat_dense_complex.h"
#include "quadrature.h"
#include "quadtree.h"

BfComplex bfHelm2GetKernelValue(BfPoint2 const xsrc, BfPoint2 const xtgt, BfVector2 const ntgt, BfReal K, BfLayerPotential layerPot);
BfSize bfHelm2RankEstForTwoCircles(BfCircle2 const *circ1, BfCircle2 const *circ2, BfReal K, BfReal C, BfReal eps);
BfMat *bfGetHelm2KernelMatrix(BfPoints2 const *Xsrc, BfPoints2 const *Xtgt, BfVectors2 const *Ntgt, BfReal K, BfLayerPotential layerPot);
BfMat *bfHelm2GetReexpansionMatrix(BfPoints2 const *srcPtsOrig, BfPoints2 const *srcPtsEquiv, BfPoints2 const *tgtPts, BfVectors2 const *tgtNormals, BfReal K, BfLayerPotential layerPot);
