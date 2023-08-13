from defs cimport BfSize
from types cimport BfMat

cdef extern from "bf/indexed_mat.h":
    cdef struct BfIndexedMat:
        BfSize i0
        BfSize j0
        BfMat *mat

    BfIndexedMat *bfIndexedMatAlloc()
