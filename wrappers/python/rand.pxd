from defs cimport BfSize

cdef extern from "bf/rand.h":
    void bfSeed(BfSize seed)
