#pragma once

#include "helm2.h"
#include "mat_product.h"
#include "quadtree.h"
#include "tree_level_iter.h"

#if BF_DEBUG
#include "points.h"
#include "vectors.h"
typedef struct {
  BfPoints2 *srcPts[2];
  BfPoints2 *tgtPts;
} BfFacAux;
#endif

BfMat *bfFacHelm2MakeSingleLevel(BfHelm2 const *helm, BfQuadtreeNode const *srcNode, BfQuadtreeNode const *tgtNode);
BfMat *bfFacHelm2MakeMultilevel(BfHelm2 const *helm, BfQuadtree const *srcTree, BfQuadtree const *tgtTree);
