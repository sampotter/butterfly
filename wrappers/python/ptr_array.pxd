from defs cimport BfPtr, BfSize
from types cimport BfPtrArray

cdef extern from "bf/ptr_array.h":
    void bfInitPtrArray(BfPtrArray *arr, BfSize capacity)
    void bfPtrArrayDeinit(BfPtrArray *arr)
    void bfPtrArrayDeinitAndDealloc(BfPtrArray **arr)
    BfSize bfPtrArraySize(const BfPtrArray *arr)
    void bfPtrArrayAppend(BfPtrArray *arr, BfPtr ptr)
    BfPtr bfPtrArrayGet(const BfPtrArray *arr, BfSize pos)
