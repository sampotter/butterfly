from defs cimport BfPolicy
from types cimport BfMat, BfMatDiff

cdef extern from "bf/mat_diff.h":
    BfMat *bfMatDiffToMat(BfMatDiff *matDiff)
    BfMatDiff *bfMatToMatDiff(BfMat *mat)

    BfMatDiff *bfMatDiffNew(BfMat *first, BfMat *second, BfPolicy policy)
    BfMat *bfMatDiffGetFirst(BfMatDiff *matDiff)
    BfMat *bfMatDiffGetSecond(BfMatDiff *matDiff)
