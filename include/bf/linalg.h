#pragma once

#include "mat.h"
#include "mat_diag_real.h"

typedef struct BfTruncSpec {
  bool usingTol;
  union {
    BfReal tol;
    BfSize k;
  };
} BfTruncSpec;

BfSize bfTruncSpecGetNumTerms(BfTruncSpec const *truncSpec, BfMatDiagReal const *S);

BfMat *bfSolveGMRES(BfMat const *A, BfMat const *B, BfMat *X0, BfReal tol, BfSize maxNumIter, BfSize *numIter, BfMat const *M);
BfReal bfGetMaxEigenvalue(BfMat const *L, BfMat const *M);
void bfGetShiftedEigs(BfMat const *A, BfMat const *M, BfReal sigma, BfSize k, BfMat **PhiPtr, BfVecReal **LambdaPtr);
void bfGetEigenband(BfMat const *A, BfMat const *M, BfReal lam0, BfReal lam1, BfReal sigma, BfMat **PhiPtr, BfVecReal **LambdaPtr);
bool bfGetTruncatedSvd(BfMat const *mat, BfMat **U, BfMatDiagReal **S, BfMat **V, BfTruncSpec const *truncSpec, BfBackend backend);
