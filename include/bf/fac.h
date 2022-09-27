#pragma once

#include "layer_pot.h"
#include "mat_product.h"
#include "quadtree.h"

BfSize
bfFacHelm2Prepare(BfQuadtreeNode const *srcNode,
                  BfQuadtreeNode const *tgtNode,
                  BfReal K,
                  BfQuadtreeLevelIter *srcLevelIter,
                  BfQuadtreeLevelIter *tgtLevelIter);

BfMatProduct *bfFacHelm2Make(BfQuadtree const *tree, BfReal K,
                             BfQuadtreeLevelIter *srcLevelIter,
                             BfQuadtreeLevelIter *tgtLevelIter,
                             BfSize numFactors);

BfMat *bfFacHelm2MakeMultilevel(BfQuadtree const *tree, BfReal K, BfLayerPotential layerPot);
