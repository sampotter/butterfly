from defs cimport BfSize, BfReal, BfPolicy
from types cimport BfMatProduct, BfTree, BfFac

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

    void bfFacDeinitAndDealloc(BfFac **fac)
    BfMatProduct *bfFacGetMatProduct(BfFac *fac, BfPolicy policy)
