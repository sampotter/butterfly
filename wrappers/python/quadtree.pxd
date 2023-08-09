from points cimport BfPoints2
from types cimport BfTree, BfQuadtree
from vectors cimport BfVectors2

cdef extern from "bf/quadtree.h":
    BfTree *bfQuadtreeToTree(BfQuadtree *quadtree)
    BfQuadtree *bfTreeToQuadtree(BfTree *tree)

    BfQuadtree *bfQuadtreeNew()
    void bfQuadtreeInit(BfQuadtree *tree, const BfPoints2 *points, const BfVectors2 *unitNormals)
    void bfQuadtreeDeinitAndDealloc(BfQuadtree **quadtree)
