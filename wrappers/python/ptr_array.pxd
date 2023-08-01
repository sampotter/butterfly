from defs cimport BfPtr, BfSize

cdef extern from "bf/ptr_array.h":
    cdef struct BfPtrArray:
        pass

    void bfInitPtrArray(BfPtrArray *arr, BfSize capacity)
    void bfPtrArrayDeinit(BfPtrArray *arr)
    BfSize bfPtrArraySize(const BfPtrArray *arr)
    void bfPtrArrayAppend(BfPtrArray *arr, BfPtr ptr)
    BfPtr bfPtrArrayGet(const BfPtrArray *arr, BfSize pos)
