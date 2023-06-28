cdef extern from "bf/perm.h":
    struct BfPerm:
        pass

    BfPerm *bfPermGetView(BfPerm *perm)
    BfPerm *bfPermGetReversePerm(const BfPerm *perm)
