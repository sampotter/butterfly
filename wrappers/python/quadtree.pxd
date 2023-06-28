from points cimport BfPoints2
from tree cimport BfTree
from vectors cimport BfVectors2

cdef extern from "bf/quadtree.h":
    struct BfQuadtree:
        pass

    BfQuadtree *bfQuadtreeNew()
    void bfQuadtreeInit(BfQuadtree *tree, const BfPoints2 *points, const BfVectors2 *unitNormals)
    void bfQuadtreeDeinitAndDealloc(BfQuadtree **quadtree)
