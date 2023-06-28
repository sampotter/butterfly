from defs cimport BfReal, BfSize
from geom cimport BfPoint2
from points cimport BfPoints2
from real_array cimport BfRealArray
from vectors cimport BfVectors2

cdef extern from "bf/ellipse.h":
    struct BfEllipse:
        BfReal semiMajorAxis
        BfReal semiMinorAxis
        BfPoint2 center
        BfReal theta

    BfReal bfEllipseGetPerimeter(const BfEllipse *ellipse)
    void bfEllipseSampleLinspaced(const BfEllipse *ellipse, BfSize numPoints, BfPoints2 *points, BfVectors2 *unitTangents, BfVectors2 *unitNormals, BfRealArray *weights)
