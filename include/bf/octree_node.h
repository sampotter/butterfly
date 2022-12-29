#include "octree.h"
#include "tree_node.h"

/** Interface: */

BfType bfOctreeNodeGetType(BfTreeNode const *treeNode);

/** Upcasting: OctreeNode -> TreeNode */

BfTreeNode *bfOctreeNodeToTreeNode(BfOctreeNode *node);
BfTreeNode const *bfOctreeNodeConstToTreeNodeConst(BfOctreeNode const *node);

/** Downcasting: TreeNode -> OctreeNode */

BfOctreeNode *bfTreeNodeToOctreeNode(BfTreeNode *treeNode);
BfOctreeNode const *bfTreeNodeConstToOctreeNodeConst(BfTreeNode const *treeNode);

/** Implementation: OctreeNode */

struct BfOctreeNode {
  BfTreeNode super;

  /* The bounding box for this octree node. All of the points indexed
   * by this node must lie in the rectangle determined by the product
   * of `[boundingBox.min[0], boundingBox.max[0]]` and
   * `[boundingBox.min[1], boundingBox.max[1]]`. */
  BfBoundingBox3 boundingBox;

  /* The split for this octree node. If a `point` is in `child[0]`,
   * it satisfies `point[0] <= split[0]` and `point[1] <=
   * split[1]`. If it's in `child[1]`, then `point[0] <= split[0]` and
   * `point[1] > split[1]`. Likewise, if it's in `child[2]`, then
   * `point[0] > split[0]` and `point[1] <= split[1]`; and if it's in
   * `child[3]`, then `point[0] > split[0]` and `point[1]` >
   * `split[1]`. */
  BfPoint3 split;
};

BfOctreeNode *bfOctreeNodeNew();
void bfOctreeNodeInitRoot(BfOctreeNode *node, BfOctree const *tree);
void bfOctreeNodeDeinit(BfOctreeNode *node);
BfSphere bfOctreeNodeGetBoundingSphere(BfOctreeNode const *node);
void bfOctreeNodeGetPoints(BfOctree const *tree, BfOctreeNode const *node, BfPoints3 *points);
void bfOctreeNodeGetUnitNormals(BfOctree const *tree, BfOctreeNode const *node, BfVectors3 *vectors);
bool bfOctreeNodeIsSeparated(BfOctreeNode const *node1, BfOctreeNode const *node2);
