from defs cimport BfSize, BfReal
from types cimport BfTree

cdef extern from "bf/fac.h":
    ctypedef struct BfFacSpec:
        BfTree *rowTree
        BfTree *colTree
        BfSize rowTreeInitDepth
        BfSize colTreeInitDepth
        BfReal tol
        BfSize minNumRows
        BfSize minNumCols
        bint compareRelativeErrors
