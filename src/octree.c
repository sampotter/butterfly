#include <bf/octree.h>

#include <bf/bbox.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/octree_node.h>
#include <bf/points.h>

#include "macros.h"

/** Interface: */

static BfTreeVtable TreeVtable = {
  .GetType = (__typeof__(&bfTreeGetType))bfOctreeGetType,
};

BfType bfOctreeGetType(BfTree const *tree) {
  (void)tree;
  return BF_TYPE_OCTREE;
}

/** Upcasting: */

BfTree *bfOctreeToTree(BfOctree *octree) {
  return &octree->super;
}

/** Implementation: */

void bfOctreeInitFromPoints(BfOctree *tree, BfPoints3 const *points, BfVectors3 const *unitNormals) {
  BEGIN_ERROR_HANDLING();

  BfOctreeNode *root = bfOctreeNodeNew();
  HANDLE_ERROR();

  bfTreeInit(&tree->super, &TreeVtable, bfOctreeNodeToTreeNode(root), points->size);
  HANDLE_ERROR();

  tree->points = points;
  tree->unitNormals = unitNormals;

  bfOctreeNodeInitRoot(root, tree);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    bfOctreeNodeDeinit(root);
    bfOctreeDeinit(tree);
  }
}

void bfOctreeDeinit(BfOctree *octree) {
  (void)octree;
  assert(false);
}
