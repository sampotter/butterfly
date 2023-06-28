cdef extern from "bf/real_array.h":
    struct BfRealArray:
        pass

    BfRealArray *bfRealArrayNew()
    void bfRealArrayInitWithDefaultCapacity(BfRealArray *realArray)
    void bfRealArrayExtend(BfRealArray *realArray, const BfRealArray *otherRealArray)
