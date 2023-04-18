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

void bf_accum_with_KR_correction(BfSize order, BfKernelComplex K, BfPtr *aux, BfSize n, BfComplex const *x, BfReal const *h, BfComplex *y);
void bf_apply_KR_correction(BfMat *mat, BfSize order, BfKernelComplex K, BfPtr *aux);
void bf_apply_KR_correction_quadtree(BfMat *mat, BfSize order, BfTree const *tree, BfKernelComplex K, BfPtr *aux);
void bf_apply_block_KR_correction(BfMat *mat, BfSizeArray const *offsets, BfSize order, BfKernelComplex K, BfPtr *aux);
void bf_apply_block_KR_correction_quadtree(BfMat *mat, BfSizeArray const *offsets, BfSize order, BfTree const *Tree, BfKernelComplex K, BfPtr *aux);
