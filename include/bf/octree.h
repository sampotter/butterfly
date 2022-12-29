#pragma once

#include <stdlib.h>

#include "bbox.h"
#include "error.h"
#include "mat.h"
#include "ptr_array.h"
#include "tree.h"

/** Interface: Tree -> Octree */

BfType bfOctreeGetType(BfTree const *tree);

/** Upcasting: Octree -> Tree */

BfTree *bfOctreeToTree(BfOctree *octree);

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

void bfOctreeInitFromPoints(BfOctree *tree, BfPoints3 const *points, BfVectors3 const *unitNormals);
void bfOctreeDeinit(BfOctree *tree);
void bfOctreeSaveBoxesToTextFile(BfOctree const *tree, char const *path);
