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

void bf_apply_KR_correction(BfMat *mat, BfSize order, BfKernelComplex K,
                            BfPtr *aux);

void bf_apply_KR_correction_quadtree(BfMat *mat, BfSize order,
                                     BfTree const *tree, BfKernelComplex K,
                                     BfPtr *aux);
