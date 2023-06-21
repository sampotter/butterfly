#include "octree.h"
#include "tree_node.h"

/** Interface: TreeNode */

BfType bfOctreeNodeGetType(BfTreeNode const *treeNode);
void bfOctreeNodeDelete(BfOctreeNode **octreeNode);

/** OctreeNode: */

/** Upcasting: OctreeNode -> TreeNode */

BfTreeNode *bfOctreeNodeToTreeNode(BfOctreeNode *node);
BfTreeNode const *bfOctreeNodeConstToTreeNodeConst(BfOctreeNode const *node);

/** Downcasting: TreeNode -> OctreeNode */

BfOctreeNode *bfTreeNodeToOctreeNode(BfTreeNode *treeNode);
BfOctreeNode const *bfTreeNodeConstToOctreeNodeConst(BfTreeNode const *treeNode);

/** Implementation: OctreeNode */

struct BfOctreeNode {
  BfTreeNode super;

  BfBoundingBox3 boundingBox;
  BfPoint3 split;
};

BfTreeNodeVtable *bfOctreeNodeGetDefaultTreeNodeVtable();
BfOctreeNode *bfOctreeNodeNew();
void bfOctreeNodeInitRoot(BfOctreeNode *node, BfOctree const *tree, BfSize maxLeafSize);
void bfOctreeNodeDeinit(BfOctreeNode *node);
void bfOctreeNodeDealloc(BfOctreeNode **node);
void bfOctreeNodeDeinitAndDealloc(BfOctreeNode **node);
BfCircle bfOctreeNodeGetBoundingCircle(BfOctreeNode const *node);
BfPoints3 bfOctreeNodeGetPoints(BfOctreeNode const *node, BfOctree const *tree);
BfVectors3 bfOctreeNodeGetUnitNormals(BfOctreeNode const *node, BfOctree const *tree);
bool bfOctreeNodesAreSeparated(BfOctreeNode const *node1, BfOctreeNode const *node2);
