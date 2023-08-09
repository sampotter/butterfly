from defs cimport BfSize
from perm cimport BfPerm
from ptr_array cimport BfPtrArray
from types cimport BfNodeSpan, BfTree, BfTreeNode, BfType

cdef extern from "bf/tree.h":
    BfType bfTreeGetType(const BfTree *tree)

    BfTree *bfTreeNewFromNodeSpan(BfNodeSpan *nodeSpan, BfPerm **perm)
    BfTree *bfTreeNewFromNode(BfTreeNode *node, BfPerm **perm)
    BfTree *bfTreeNewForMiddleFac(const BfTree *templateTree, BfSize p)
    BfTreeNode *bfTreeGetRootNode(BfTree *)
    BfPerm *bfTreeGetPerm(BfTree *tree)
    BfSize bfTreeGetMaxDepth(const BfTree *tree)
    BfPtrArray *bfTreeGetLevelPtrArray(BfTree *tree, BfSize level)
