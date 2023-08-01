from defs cimport BfSize, BfPolicy
from ptr_array cimport BfPtrArray
from types cimport BfMat, BfMatBlockDense, BfVec

cdef extern from "bf/mat_block_dense.h":
    BfMatBlockDense *bfMatBlockDenseNewFromBlocks(BfSize numRowBlocks, BfSize numColBlocks, BfPtrArray *blocks, BfPolicy policy)
    void bfMatBlockDenseScaleCols(BfMatBlockDense *matBlockDense, const BfVec *vec)
    BfMat *bfMatBlockDenseGetBlock(BfMatBlockDense *matBlockDense, BfSize i, BfSize j)
    BfMat *bfMatBlockDenseToMat(BfMatBlockDense *matBlockDense)
    BfMatBlockDense *bfMatToMatBlockDense(BfMat *mat)
