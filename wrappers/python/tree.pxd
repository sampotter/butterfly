from defs cimport BfSize
from perm cimport BfPerm
from ptr_array cimport BfPtrArray
from types cimport BfTree

cdef extern from "bf/tree.h":
    BfPerm *bfTreeGetPerm(BfTree *tree)
    BfPtrArray bfTreeGetLevelPtrArray(BfTree *tree, BfSize level)
