from types cimport BfMat, BfMatDenseComplex

cdef extern from "bf/mat_dense_complex.h":
    BfMat *bfMatDenseComplexToMat(BfMatDenseComplex *matDenseComplex)
    BfMatDenseComplex *bfMatToMatDenseComplex(BfMat *mat)
