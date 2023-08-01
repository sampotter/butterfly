ctypedef double BfReal
ctypedef size_t BfSize
ctypedef double complex BfComplex
ctypedef void *BfPtr

cdef extern from "bf/def.h":
    cdef enum BfPolicy:
        BF_POLICY_VIEW
        BF_POLICY_COPY
        BF_POLICY_STEAL
