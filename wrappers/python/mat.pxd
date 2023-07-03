from types cimport BfMat, BfType

cdef extern from "bf/mat.h":
    BfType bfMatGetType(const BfMat *mat)
