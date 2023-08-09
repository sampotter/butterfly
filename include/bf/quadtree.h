#pragma once

#include "bbox.h"
#include "error.h"
#include "mat.h"
#include "ptr_array.h"
#include "tree.h"

typedef void (*BfQuadtreeMapFunc)(BfQuadtree *, BfQuadtreeNode *, void *);
typedef void (*BfQuadtreeMapConstFunc)(BfQuadtree const *, BfQuadtreeNode const *, void *);

/** Quadtree: */

struct BfQuadtree {
  BfTree super;

  /* A pointer to the point set upon which the quadtree has been
   * built. The points in this array in the original order, and have
   * not been permuted into the quadtree order. */
  BfPoints2 const *points;

  /* Pointer to an array of unit normals, each associated with a point
   * in `points`. */
  BfVectors2 const *unitNormals;
};

/** Interface(Tree, Quadtree) */

void bfQuadtreeCopyInto(BfQuadtree *quadtree, BfQuadtree *dstQuadtree);
BfType bfQuadtreeGetType(BfQuadtree const *);

/** Upcasting: Quadtree -> Tree */

BfTree *bfQuadtreeToTree(BfQuadtree *quadtree);
BfTree const *bfQuadtreeConstToTreeConst(BfQuadtree const *quadtree);

/** Downcasting: Tree -> Quadtree */

BfQuadtree *bfTreeToQuadtree(BfTree *tree);
BfQuadtree const *bfTreeConstToQuadtreeConst(BfTree const *tree);

/** Implementation: Quadtree */

BfQuadtree *bfQuadtreeNew(void);
void bfQuadtreeInit(BfQuadtree *tree, BfPoints2 const *points, BfVectors2 const *unitNormals);
void bfQuadtreeDeinit(BfQuadtree *tree);
void bfQuadtreeDealloc(BfQuadtree **quadtree);
void bfQuadtreeDeinitAndDealloc(BfQuadtree **quadtree);
void bfQuadtreeSaveBoxesToTextFile(BfQuadtree const *tree, char const *path);
