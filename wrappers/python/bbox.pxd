from geom cimport BfPoint2

cdef extern from "bf/bbox.h":
    struct BfBbox2:
        BfPoint2 min
        BfPoint2 max
