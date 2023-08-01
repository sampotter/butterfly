from defs cimport BfReal, BfSize
from types cimport BfMat

cdef extern from "bf/linalg.h":
    BfMat *bfSolveGMRES(const BfMat *A, const BfMat *B, BfMat *X0, BfReal tol, BfSize maxNumIter, BfSize *numIter, const BfMat *M)
