from defs cimport BfReal

cdef extern from "bf/geom.h":
    ctypedef BfReal BfPoint2[2]
