from defs cimport BfSize

cdef extern from "bf/size_array.h":
    cdef struct BfSizeArray:
        pass

    BfSizeArray *bfSizeArrayNewWithCapacity(BfSize capacity)
    BfSizeArray *bfSizeArrayNewView(BfSize n, BfSize *inds)
    void bfSizeArrayDeinitAndDealloc(BfSizeArray **sizeArray)
    void bfSizeArrayAppend(BfSizeArray *sizeArray, BfSize elt)
    BfSize bfSizeArrayGetSize(const BfSizeArray *sizeArray)
    BfSize *bfSizeArrayGetDataPtr(const BfSizeArray *sizeArray)
