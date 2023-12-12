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

/* Compute the first `k` eigenpairs of the symmetric generalized eigenvalue
 * problem:
 *
 *     A phi = lam_k M phi
 *
 * which are closest in magnitude to `sigma`.
 *
 * NOTE: we return the matrix Phi_k^T, where:
 *
 *     Phi_k = [phi_1 ... phi_k]
 *
 * since we store dense matrices with row major ordering. This makes
 * some things more convenient. */
void bfGetShiftedEigs(BfMat const *A, BfMat const *M, BfReal sigma, BfSize k, BfMat **PhiTransposePtr, BfVecReal **LambdaPtr);

typedef enum BfEigenbandMethod {
  BF_EIGENBAND_METHOD_DOUBLING,
  BF_EIGENBAND_METHOD_COVERING
} BfEigenbandMethod;

/* Compute all eigenpairs of the symmetric generalized eigenvalue problem:
 *
 *     A phi = lam M phi
 *
 * such that `lam` lies in `interval`. */
void bfGetEigenband(BfMat const *A, BfMat const *M, BfInterval const *interval, BfEigenbandMethod method, BfMat **PhiPtr, BfVecReal **LambdaPtr);

bool bfGetTruncatedSvd(BfMat const *mat, BfMat **U, BfMatDiagReal **S, BfMat **V, BfTruncSpec const *truncSpec, BfBackend backend);
