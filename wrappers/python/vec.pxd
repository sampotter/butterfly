from types cimport BfVec

cdef extern from "bf/vec.h":
    void bfVecDelete(BfVec **)
