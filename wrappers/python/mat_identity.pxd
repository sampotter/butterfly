from defs cimport BfSize
from types cimport BfMat, BfMatIdentity

cdef extern from "bf/mat_identity.h":
    BfMat *bfMatIdentityToMat(BfMatIdentity *matIdentity);
    BfMatIdentity *bfMatIdentityNew()
    void bfMatIdentityInit(BfMatIdentity *matIdentity, BfSize n)
