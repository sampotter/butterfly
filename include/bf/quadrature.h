#pragma once

#include <bf/mat.h>

typedef enum BfQuadratures {
  BF_QUADRATURE_KAPUR_ROKHLIN,
} BfQuadrature;

void bf_apply_KR_correction(BfMat *mat, BfSize order);
