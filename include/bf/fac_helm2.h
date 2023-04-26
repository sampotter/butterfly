#pragma once

#include "layer_pot.h"
#include "mat_product.h"
#include "quadtree.h"
#include "tree_level_iter.h"

#if BF_DEBUG
#include "points.h"
#include "vectors.h"
typedef struct {
  BfPoints2 srcPts[2];
  BfPoints2 tgtPts;
} BfFacAux;
#endif

BfSize bfFacHelm2Prepare(BfQuadtreeNode const *srcNode, BfQuadtreeNode const *tgtNode, BfReal K, BfTreeLevelIter *srcLevelIter, BfTreeLevelIter *tgtLevelIter);
BfMatProduct *bfFacHelm2Make(BfQuadtree const *srcTree, BfQuadtree const *tgtTree, BfReal K, BfLayerPotential layerPot, BfComplex const *alpha, BfComplex const *beta, BfTreeLevelIter *srcLevelIter, BfTreeLevelIter *tgtLevelIter, BfSize numFactors);
BfMat *bfFacHelm2MakeMultilevel(BfQuadtree const *srcTree, BfQuadtree const *tgtTree, BfReal K, BfLayerPotential layerPot, BfComplex const *alpha, BfComplex const *beta);
