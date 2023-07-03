from types cimport BfMat, BfMatBlockDense

cdef extern from "bf/mat_block_dense.h":
    BfMat *bfMatBlockDenseToMat(BfMatBlockDense *matBlockDense)
    BfMatBlockDense *bfMatToMatBlockDense(BfMat *mat)
