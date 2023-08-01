from defs cimport BfSize
from types cimport BfMat, BfMatFunc

cdef extern from "bf/mat_func.h":
    ctypedef BfMat *(*MatMulFunc)(const BfMat *, void *);

    BfMatFunc *bfMatFuncNew();
    void bfMatFuncInit(BfMatFunc *matFunc, BfSize numRows, BfSize numCols, MatMulFunc matMul, void *);
