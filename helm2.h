#pragma once

#include "def.h"
#include "error.h"
#include "geom.h"

BfComplex bfHelm2GetKernelValue(BfPoint2 const p, BfPoint2 const q, BfReal k);

BfReal bfHelm2RankEstForTwoCircles(BfCircle2 const circ1, BfCircle2 const circ2,
                                   BfReal k, BfReal C, BfReal eps);

enum BfError
bfHelm2KernelMatrixFromPoints(BfMat *K, BfMat const *X, BfMat const *Y, BfReal k);
