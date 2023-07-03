from perm cimport BfPerm
from types cimport BfTree

cdef extern from "bf/tree.h":
    BfPerm *bfTreeGetPerm(BfTree *tree)
