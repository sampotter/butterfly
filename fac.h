#pragma once

#include "mat_block_coo.h"
#include "quadtree.h"

typedef struct BfFactor {
  BfMatBlockCoo *mat;
#if BF_DEBUG
  BfPoints2 *srcPtsOrig, *srcPtsEquiv, *tgtPts;
#endif
} BfFactor;

void bfFreeFactor(BfFactor *factor);

void bfMakeFac (BfQuadtree const *tree,
                BfQuadtreeNode const *srcNode, BfQuadtreeNode const *tgtNode,
                BfReal K, BfSize *numFactors, BfFactor **factors);

void bfFreeFac(BfSize numFactors, BfFactor **factorPtr);

void bfMulFac(BfFactor const *factor, BfMat const *X, BfMat *Y);
