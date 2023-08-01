from defs cimport BfReal, BfComplex, BfSize
from layer_pot cimport BfLayerPotential
from points cimport BfPoints2
from size_array cimport BfSizeArray
from types cimport BfMat, BfTree
from vectors cimport BfVectors2

cdef extern from "bf/helm2.h":
    cdef struct BfHelm2:
        BfReal k
        BfLayerPotential layerPot
        BfComplex alpha
        BfComplex beta

    BfMat *bfHelm2GetKernelMatrix(const BfHelm2 *helm, const BfPoints2 *Xsrc, const BfPoints2 *Xtgt, const BfVectors2 *Nsrc, const BfVectors2 *Ntgt)
    void bfHelm2ApplyKrCorrection(const BfHelm2 *helm2, BfSize krOrder, const BfPoints2 *pts, const BfVectors2 *normals, BfMat *mat)
    void bfHelm2ApplyKrCorrectionTree(const BfHelm2 *helm2, BfSize krOrder, const BfPoints2 *pts, const BfVectors2 *normals, const BfTree *tree, BfMat *mat)
    void bfHelm2ApplyBlockCorrection(const BfHelm2 *helm2, const BfSizeArray *offsets, BfSize krOrder, const BfPoints2 *pts, const BfVectors2 *normals, BfMat *mat)
    void bfHelm2ApplyBlockCorrectionTree(const BfHelm2 *helm2, const BfSizeArray *offsets, BfSize krOrder, const BfPoints2 *pts, const BfVectors2 *normals, const BfTree *tree, BfMat *mat)
