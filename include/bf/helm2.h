#pragma once

#include "def.h"
#include "error.h"
#include "geom.h"
#include "mat_dense_complex.h"
#include "quadrature.h"
#include "quadtree.h"

BfComplex bfHelm2GetKernelValue(BfPoint2 const xsrc, BfPoint2 const xtgt, BfVector2 const nsrc, BfVector2 const ntgt, BfReal K, BfLayerPotential layerPot);
BfSize bfHelm2RankEstForTwoCircles(BfCircle const *circ1, BfCircle const *circ2, BfReal K, BfReal C, BfReal eps);
BfMat *bfHelm2GetKernelMatrix(BfPoints2 const *Xsrc, BfPoints2 const *Xtgt, BfVectors2 const *Nsrc, BfVectors2 const *Ntgt, BfReal K, BfLayerPotential layerPot, BfComplex const *alpha, BfComplex const *beta);
BfMat *bfHelm2GetReexpansionMatrix(BfPoints2 const *srcPtsOrig, BfPoints2 const *srcPtsEquiv, BfVectors2 const *srcNormalsOrig, BfVectors2 const *srcNormalsEquiv, BfPoints2 const *tgtPts, BfReal K, BfLayerPotential layerPot, BfComplex const *alpha, BfComplex const *beta);
