from perm cimport BfPerm
from tree cimport BfTree

cdef extern from "bf/tree.h":
    struct BfTree:
        pass

    BfPerm *bfTreeGetPerm(BfTree *tree)
