from defs cimport BfPolicy
from types cimport BfMat, BfMatProduct

cdef extern from "bf/mat_product.h":
    BfMat *bfMatProductToMat(BfMatProduct *matProduct)
    BfMatProduct *bfMatToMatProduct(BfMat *mat)
    BfMatProduct *bfMatProductNew()
    void bfMatProductInit(BfMatProduct *prod)
    void bfMatProductPostMultiply(BfMatProduct *prod, BfMat *mat)
    BfMat *bfMatProductGetTransposed(BfMatProduct *matProduct, BfPolicy policy)
