from tree_traversals cimport BfTreeTraversal
from types cimport BfPtrArray, BfTree, BfTreeLevelIter

cdef extern from "bf/tree_level_iter.h":
    BfTreeLevelIter *bfTreeLevelIterNewFromTree(BfTreeTraversal traversal, BfTree *tree)
    bint bfTreeLevelIterIsDone(const BfTreeLevelIter *iter)
    void bfTreeLevelIterNext(BfTreeLevelIter *iter)
    BfPtrArray *bfTreeLevelIterGetLevelNodes(BfTreeLevelIter *iter)
