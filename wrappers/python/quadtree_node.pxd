from bbox cimport BfBbox2
from geom cimport BfPoint2
from points cimport BfPoints2
from quadtree cimport BfQuadtree
from tree_node cimport BfTreeNode

cdef extern from "bf/quadtree_node.h":
    cdef struct BfQuadtreeNode:
        pass

    BfTreeNode *bfQuadtreeNodeToTreeNode(BfQuadtreeNode *node)
    BfQuadtreeNode *bfTreeNodeToQuadtreeNode(BfTreeNode *treeNode)
    BfPoints2 *bfQuadtreeNodeGetPoints(const BfQuadtreeNode *node, const BfQuadtree *tree)
    BfQuadtree *bfQuadtreeNodeGetQuadtree(BfQuadtreeNode *node)
    BfBbox2 bfQuadtreeNodeGetBbox(const BfQuadtreeNode *node)
    void bfQuadtreeNodeGetSplit(const BfQuadtreeNode *node, BfPoint2 split)
