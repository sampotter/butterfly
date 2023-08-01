from defs cimport BfSize

cdef extern from "bf/tree_node.h":
    cdef struct BfTreeNode:
        pass

    BfTreeNode *bfTreeNodeGetParent(BfTreeNode *treeNode)
    BfSize bfTreeNodeGetNumPoints(const BfTreeNode *treeNode)
    BfSize bfTreeNodeGetFirstIndex(const BfTreeNode *treeNode)
    BfSize bfTreeNodeGetLastIndex(const BfTreeNode *treeNode)
