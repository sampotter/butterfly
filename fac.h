#pragma once

#include "mat_product.h"
#include "quadtree.h"

BfSize
bfFacHelm2Prepare(BfQuadtree const *tree, BfQuadtreeNode const *srcNode,
                  BfQuadtreeNode const *tgtNode, BfReal K,
                  BfQuadtreeLevelIter *srcLevelIter,
                  BfQuadtreeLevelIter *tgtLevelIter);

BfMatProduct *bfFacHelm2Make(BfQuadtree const *tree,
                             BfQuadtreeNode const *srcNode,
                             BfQuadtreeNode const *tgtNode,
                             BfReal K,
                             BfQuadtreeLevelIter *srcLevelIter,
                             BfQuadtreeLevelIter *tgtLevelIter,
                             BfSize numFactors);
