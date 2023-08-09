from defs cimport BfSize
from types cimport BfTree, BfTreeNode, BfType

cdef extern from "bf/tree_node.h":
    BfType bfTreeNodeGetType(const BfTreeNode *treeNode)
    BfTreeNode *bfTreeNodeGetParent(BfTreeNode *treeNode)
    BfSize bfTreeNodeGetDepth(const BfTreeNode *treeNode)
    BfSize bfTreeNodeGetNumPoints(const BfTreeNode *treeNode)
    BfSize bfTreeNodeGetFirstIndex(const BfTreeNode *treeNode)
    BfSize bfTreeNodeGetLastIndex(const BfTreeNode *treeNode)
    BfTree *bfTreeNodeGetTree(BfTreeNode *treeNode)
    BfTreeNode *bfTreeNodeGetChild(BfTreeNode *treeNode, BfSize i)
    BfSize bfTreeNodeGetMaxNumChildren(const BfTreeNode *treeNode)
    bint bfTreeNodeHasChild(const BfTreeNode *treeNode, BfSize i)
