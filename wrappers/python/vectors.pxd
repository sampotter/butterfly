from defs cimport BfSize

cdef extern from "bf/vectors.h":
    struct BfVectors2:
        pass

    BfVectors2 *bfVectors2NewEmpty()
    void bfVectors2Extend(BfVectors2 *vectors, const BfVectors2 *newVectors)
