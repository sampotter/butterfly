from types cimport BfMat, BfMatProduct

cdef extern from "bf/mat_product.h":
    BfMat *bfMatProductToMat(BfMatProduct *matProduct)
    BfMatProduct *bfMatProductNew()
    void bfMatProductInit(BfMatProduct *prod)
    void bfMatProductPostMultiply(BfMatProduct *prod, BfMat *mat)
