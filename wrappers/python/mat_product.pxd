from defs cimport BfSize, BfPolicy
from types cimport BfMat, BfMatProduct

cdef extern from "bf/mat_product.h":
    BfMat *bfMatProductToMat(BfMatProduct *matProduct)
    BfMatProduct *bfMatToMatProduct(BfMat *mat)
    BfMatProduct *bfMatProductNew()
    void bfMatProductInit(BfMatProduct *prod)
    BfSize bfMatProductNumFactors(const BfMatProduct *prod)
    BfMat *bfMatProductGetFactor(BfMatProduct *prod, BfSize i)
    void bfMatProductPostMultiply(BfMatProduct *prod, BfMat *mat)
