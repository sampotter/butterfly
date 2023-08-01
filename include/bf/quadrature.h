#pragma once

#include <bf/layer_pot.h>
#include <bf/mat.h>
#include <bf/points.h>
#include <bf/quadtree.h>
#include <bf/vectors.h>

typedef enum BfQuadratures {
  BF_QUADRATURE_KAPUR_ROKHLIN,
} BfQuadrature;

typedef BfComplex (*BfKernelComplex)(BfSize, BfSize, void *);

void bfQuadKrAccumCorrection(BfSize order, BfKernelComplex K, BfPtr *aux, BfSize n, BfComplex const *x, BfReal const *h, BfComplex *y);
void bfQuadKrApplyCorrection(BfMat *mat, BfSize order, BfKernelComplex K, BfPtr *aux);
void bfQuadKrApplyCorrectionTree(BfMat *mat, BfSize order, BfTree const *tree, BfKernelComplex K, BfPtr *aux);
void bfQuadKrApplyBlockCorrection(BfMat *mat, BfSizeArray const *offsets, BfSize order, BfKernelComplex K, BfPtr *aux);
void bfQuadKrApplyBlockCorrectionTree(BfMat *mat, BfSizeArray const *offsets, BfSize order, BfTree const *tree, BfKernelComplex K, BfPtr *aux);
