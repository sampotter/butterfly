#include "quadtree.h"
#include "tree_node.h"

/** Interface: TreeNode */

BfQuadtreeNode *bfQuadtreeNodeCopy(BfQuadtreeNode const *quadtreeNode);
BfType bfQuadtreeNodeGetType(BfTreeNode const *treeNode);
void bfQuadtreeNodeDelete(BfQuadtreeNode **quadtreeNode);

/** Upcasting: QuadtreeNode -> TreeNode */

BfTreeNode *bfQuadtreeNodeToTreeNode(BfQuadtreeNode *node);
BfTreeNode const *bfQuadtreeNodeConstToTreeNodeConst(BfQuadtreeNode const *node);

/** Downcasting: TreeNode -> QuadtreeNode */

BfQuadtreeNode *bfTreeNodeToQuadtreeNode(BfTreeNode *treeNode);
BfQuadtreeNode const *bfTreeNodeConstToQuadtreeNodeConst(BfTreeNode const *treeNode);

/** Implementation: QuadtreeNode */

struct BfQuadtreeNode {
  BfTreeNode super;

  /* The bounding box for this quadtree node. All of the points
   * indexed by this node must lie in the rectangle determined by the
   * product of `[bbox.min[0], bbox.max[0]]` and `[bbox.min[1],
   * bbox.max[1]]`. */
  BfBbox2 bbox;

  /* The split for this quadtree node. If a `point` is in `child[0]`,
   * it satisfies `point[0] <= split[0]` and `point[1] <=
   * split[1]`. If it's in `child[1]`, then `point[0] <= split[0]` and
   * `point[1] > split[1]`. Likewise, if it's in `child[2]`, then
   * `point[0] > split[0]` and `point[1] <= split[1]`; and if it's in
   * `child[3]`, then `point[0] > split[0]` and `point[1]` >
   * `split[1]`. */
  BfPoint2 split;
};

BfTreeNodeVtable *bfQuadtreeNodeGetDefaultTreeNodeVtable(void);
BfQuadtreeNode *bfQuadtreeNodeNew(void);
void bfQuadtreeNodeInitRoot(BfQuadtreeNode *node, BfQuadtree const *tree);
void bfQuadtreeNodeDeinit(BfQuadtreeNode *node);
void bfQuadtreeNodeDealloc(BfQuadtreeNode **node);
void bfQuadtreeNodeDeinitAndDealloc(BfQuadtreeNode **node);
BfCircle bfQuadtreeNodeGetBoundingCircle(BfQuadtreeNode const *node);
BfPoints2 *bfQuadtreeNodeGetPoints(BfQuadtreeNode const *node, BfQuadtree const *tree);
BfVectors2 *bfQuadtreeNodeGetUnitNormals(BfQuadtreeNode const *node, BfQuadtree const *tree);
bool bfQuadtreeNodesAreSeparated(BfQuadtreeNode const *node1, BfQuadtreeNode const *node2);
BfQuadtree *bfQuadtreeNodeGetQuadtree(BfQuadtreeNode *node);
BfQuadtree const *bfQuadtreeNodeGetQuadtreeConst(BfQuadtreeNode const *node);
BfBbox2 bfQuadtreeNodeGetBbox(BfQuadtreeNode const *node);
void bfQuadtreeNodeGetSplit(BfQuadtreeNode const *node, BfPoint2 split);
