#pragma once

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

BfMatBlockDense *bfFacHelm2MakeMultilevel(BfQuadtree const *tree, BfReal K);
