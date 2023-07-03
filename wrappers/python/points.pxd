from bbox cimport BfBbox2
from defs cimport BfReal, BfSize
from geom cimport BfPoint2

cdef extern from "bf/points.h":
    struct BfPoints2:
        pass

cdef extern from "bf/poisson_disk_sampling.h":
    BfPoints2 *bfPoints2NewEmpty()
    void bfPoints2DeinitAndDealloc(BfPoints2 **points)
    BfPoints2 *bfPoints2SamplePoissonDisk(const BfBbox2 *bbox, BfReal minDist, BfSize k)
    void bfPoints2Append(BfPoints2 *points, const BfPoint2 p)
    void bfPoints2Extend(BfPoints2 *points, const BfPoints2 *newPoints)
    void bfPoints2Get(const BfPoints2 *points, BfSize i, BfPoint2 p)
    BfSize bfPoints2GetSize(const BfPoints2 *points)
