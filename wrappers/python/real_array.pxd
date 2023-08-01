from defs cimport BfReal, BfSize
from perm cimport BfPerm
from types cimport BfVec

cdef extern from "bf/real_array.h":
    struct BfRealArray:
        pass

    BfRealArray *bfRealArrayNew()
    void bfRealArrayInitWithDefaultCapacity(BfRealArray *realArray)
    void bfRealArrayInitCopy(BfRealArray *realArray, const BfRealArray *otherRealArray)
    void bfRealArrayExtend(BfRealArray *realArray, const BfRealArray *otherRealArray)
    BfVec *bfRealArrayGetVecView(BfRealArray *realArray)
    BfSize bfRealArrayGetSize(const BfRealArray *realArray)
    void bfRealArrayPermute(BfRealArray *realArray, const BfPerm *perm)
    BfReal *bfRealArrayGetDataPtr(BfRealArray *realArray)
