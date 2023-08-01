from defs cimport BfReal, BfComplex
from helm2 cimport BfHelm2
from layer_pot cimport BfLayerPotential
from types cimport BfMat, BfQuadtree

cdef extern from "bf/fac_helm2.h":
    BfMat *bfFacHelm2MakeMultilevel(const BfHelm2 *helm, const BfQuadtree *srcTree, const BfQuadtree *tgtTree)
