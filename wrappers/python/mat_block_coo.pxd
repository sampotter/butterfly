from defs cimport BfSize, BfPolicy
from types cimport BfPtrArray, BfMat, BfMatBlockCoo

cdef extern from "bf/mat_block_coo.h":
    BfMat *bfMatBlockCooToMat(BfMatBlockCoo *matBlockCoo)

    BfMatBlockCoo *bfMatBlockCooNewFromIndexedBlocks(BfSize numRows, BfSize numCols, BfPtrArray *indexedBlocks, BfPolicy policy)
