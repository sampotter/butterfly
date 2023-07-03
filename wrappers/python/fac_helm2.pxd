from defs cimport BfReal, BfComplex
from layer_pot cimport BfLayerPotential
from types cimport BfMat, BfQuadtree

cdef extern from "bf/fac_helm2.h":
    BfMat *bfFacHelm2MakeMultilevel(const BfQuadtree *srcTree, const BfQuadtree *tgtTree, BfReal K, BfLayerPotential layerPot, const BfComplex *alpha, const BfComplex *beta)
