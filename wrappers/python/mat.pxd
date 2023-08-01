from cython cimport typeof

from defs cimport BfSize
from types cimport BfMat, BfType

cdef extern from "bf/mat.h":
    BfType bfMatGetType(const BfMat *mat)
    BfSize bfMatGetNumRows(const BfMat *mat)
    BfSize bfMatGetNumCols(const BfMat *mat)
    void bfMatAddInplace(BfMat *mat, const BfMat *otherMat)
    BfMat *bfMatToType(const BfMat *mat, BfType type)
