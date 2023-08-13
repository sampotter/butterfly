from defs cimport BfPolicy
from types cimport BfPtrArray, BfMat, BfMatBlockDiag

cdef extern from "bf/mat_block_diag.h":
    BfMat *bfMatBlockDiagToMat(BfMatBlockDiag *matBlockDiag)

    BfMatBlockDiag *bfMatBlockDiagNewFromBlocks(BfPtrArray *blocks, BfPolicy policy)
