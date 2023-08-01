from defs cimport BfReal, BfSize

cdef extern from "bf/rand.h":
    void bfSeed(BfSize seed)
    BfReal bfRealUniform1()
