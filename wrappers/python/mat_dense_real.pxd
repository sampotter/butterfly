from defs cimport BfReal, BfSize
from types cimport BfMat, BfMatDenseReal

cdef extern from "bf/mat_dense_real.h":
    BfMat *bfMatDenseRealToMat(BfMatDenseReal *matDenseReal)
    BfMatDenseReal *bfMatDenseRealNewViewFromPtr(BfSize numRows, BfSize numCols, BfReal *data, BfSize rowStride, BfSize colStride)
