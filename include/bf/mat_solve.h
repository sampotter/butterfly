#pragma once

#include <bf/mat.h>

BfMat *bfMatSolveGMRES(BfMat const *A, BfMat const *B, BfMat *X0, BfReal tol,
                       BfSize maxNumIter, BfSize *numIter);
