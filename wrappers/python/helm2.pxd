from defs cimport BfReal, BfComplex
from layer_pot cimport BfLayerPotential
from points cimport BfPoints2
from types cimport BfMat
from vectors cimport BfVectors2

cdef extern from "bf/helm2.h":
    BfMat *bfHelm2GetKernelMatrix(const BfPoints2 *Xsrc, const BfPoints2 *Xtgt, const BfVectors2 *Nsrc, const BfVectors2 *Ntgt, BfReal K, BfLayerPotential layerPot, const BfComplex *alpha, const BfComplex *beta)
