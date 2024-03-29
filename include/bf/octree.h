#pragma once

#include "bbox.h"
#include "error.h"
#include "mat.h"
#include "ptr_array.h"
#include "tree.h"

/** Octree: */

struct BfOctree {
  BfTree super;

  /* A pointer to the point set upon which the octree has been
   * built. The points in this array in the original order, and have
   * not been permuted into the octree order. */
  BfPoints3 const *points;

  /* Pointer to an array of unit normals, each associated with a point
   * in `points`. */
  BfVectors3 const *unitNormals;
};

/** Interface(Tree, Octree) */

BfType bfOctreeGetType(BfOctree const *);

/** Upcasting: Octree -> Tree */

BfTree *bfOctreeToTree(BfOctree *octree);
BfTree const *bfOctreeConstToTreeConst(BfOctree const *octree);

/** Downcasting: Tree -> Octree */

BfOctree *bfTreeToOctree(BfTree *tree);

BfOctree *bfOctreeNew(void);
BfOctree *bfOctreeNewFromPoints(BfPoints3 const *points, BfSize maxLeafSize);
void bfOctreeInit(BfOctree *tree, BfPoints3 const *points, BfVectors3 const *unitNormals, BfSize maxLeafSize);
void bfOctreeDeinit(BfOctree *tree);
void bfOctreeDealloc(BfOctree **octree);
void bfOctreeDelete(BfOctree **octree);
void bfOctreeSaveBoxesToTextFile(BfOctree const *tree, char const *path);
BfSizeArray *bfOctreeGetNearestNeighbors(BfOctree const *tree, BfPoint3 const point, BfSize numNbs);
