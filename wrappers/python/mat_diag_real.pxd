from defs cimport BfSize, BfReal
from types cimport BfMat, BfMatDiagReal

cdef extern from "bf/mat_diag_real.h":
    BfMat *bfMatDiagRealToMat(BfMatDiagReal *mat)
    BfMatDiagReal *bfMatDiagRealNewConstant(BfSize numRows, BfSize numCols, BfReal diagValue)
