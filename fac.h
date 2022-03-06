#pragma once

#include "mat_product.h"
#include "quadtree.h"

BfMatProduct *bfFacMakeSingleLevelHelm2(BfQuadtree const *tree,
                                        BfQuadtreeNode const *srcNode,
                                        BfQuadtreeNode const *tgtNode,
                                        BfReal K);
