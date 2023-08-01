#pragma once

#include "def.h"
#include "error.h"
#include "geom.h"
#include "mat_dense_complex.h"
#include "quadrature.h"
#include "quadtree.h"

typedef struct BfHelm2 {
  BfReal k;
  BfLayerPotential layerPot;
  BfComplex alpha;
  BfComplex beta;
} BfHelm2;

BfComplex bfHelm2GetKernelValue(BfHelm2 const *helm, BfPoint2 const xsrc, BfPoint2 const xtgt, BfVector2 const nsrc, BfVector2 const ntgt);
BfSize bfHelm2RankEstForTwoCircles(BfHelm2 const *helm, BfCircle const *circ1, BfCircle const *circ2, BfReal C, BfReal eps);
BfMat *bfHelm2GetKernelMatrix(BfHelm2 const *helm, BfPoints2 const *Xsrc, BfPoints2 const *Xtgt, BfVectors2 const *Nsrc, BfVectors2 const *Ntgt);
BfMat *bfHelm2GetReexpansionMatrix(BfHelm2 const *helm, BfPoints2 const *srcPtsOrig, BfPoints2 const *srcPtsEquiv, BfVectors2 const *srcNormalsOrig, BfVectors2 const *srcNormalsEquiv, BfPoints2 const *tgtPts);
void bfHelm2ApplyKrCorrection(BfHelm2 const *helm, BfSize krOrder, BfPoints2 const *points, BfVectors2 const *normals, BfMat *mat);
void bfHelm2ApplyKrCorrectionTree(BfHelm2 const *helm, BfSize krOrder, BfPoints2 const *points, BfVectors2 const *normals, BfTree const *tree, BfMat *mat);
void bfHelm2ApplyBlockCorrection(BfHelm2 const *helm, BfSizeArray const *offsets, BfSize krOrder, BfPoints2 const *points, BfVectors2 const *normals, BfMat *mat);
void bfHelm2ApplyBlockCorrectionTree(BfHelm2 const *helm, BfSizeArray const *offsets, BfSize krOrder, BfPoints2 const *points, BfVectors2 const *normals, BfTree const *tree, BfMat *mat);
