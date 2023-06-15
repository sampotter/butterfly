#include <bf/vec.h>

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/vec_complex.h>
#include <bf/vec_real.h>

/** Interface: Vec */

BfVec *bfVecCopy(BfVec const *vec) {
  return vec->vtbl->Copy(vec);
}

void bfVecDelete(BfVec **vec) {
  (*vec)->vtbl->Delete(vec);
}

BfType bfVecGetType(BfVec const *vec) {
  return vec->vtbl->GetType(vec);
}

BfPtr bfVecGetEltPtr(BfVec *vec, BfSize i) {
  return vec->vtbl->GetEltPtr(vec, i);
}

BfVec *bfVecGetSubvecCopy(BfVec const *vec, BfSize i0, BfSize i1) {
  return vec->vtbl->GetSubvecCopy(vec, i0, i1);
}

BfVec *bfVecGetSubvecView(BfVec *vec, BfSize i0, BfSize i1) {
  return vec->vtbl->GetSubvecView(vec, i0, i1);
}

BfVec const *bfVecGetSubvecViewConst(BfVec const *vec, BfSize i0, BfSize i1) {
  return vec->vtbl->GetSubvecViewConst(vec, i0, i1);
}

void bfVecPrint(BfVec const *vec, FILE *fp) {
  vec->vtbl->Print(vec, fp);
}

/* Set the elements of `vec` in the index range [i0, i1) to the values
 * given by `otherVec`. Sets `BF_ERROR_INVALID_ARGUMENTS` if
 * `otherVec` is the wrong size. */
void bfVecSetRange(BfVec *vec, BfSize i0, BfSize i1, BfVec const *otherVec) {
  vec->vtbl->SetRange(vec, i0, i1, otherVec);
}

void bfVecSetMask(BfVec *vec, bool const *mask, BfVec const *otherVec) {
  vec->vtbl->SetMask(vec, mask, otherVec);
}

/* Compute the Euclidean distance between `vec` and `otherVec`. */
BfReal bfVecDist(BfVec const *vec, BfVec const *otherVec) {
  return vec->vtbl->Dist(vec, otherVec);
}

/* Compute the $\ell_\infty$ distance between `vec` and `otherVec`. */
BfReal bfVecDistMax(BfVec const *vec, BfVec const *otherVec) {
  return vec->vtbl->DistMax(vec, otherVec);
}

/* Compute the Euclidean norm of `vec`. */
BfReal bfVecNormMax(BfVec const *vec) {
  return vec->vtbl->NormMax(vec);
}

/* Scale each element of `vec` by `factor`. */
void bfVecScaleByReal(BfVec *vec, BfReal factor) {
  vec->vtbl->ScaleByReal(vec, factor);
}

/* Accumulate `otherVec` componentwise into `vec`. */
void bfVecAddInplace(BfVec *vec, BfVec const *otherVec) {
  vec->vtbl->AddInplace(vec, otherVec);
}

/* Compute the left matrix-vector product `mat`*`vec` inplace,
 * overwriting `vec` with the results. */
void bfVecMulInplace(BfVec *vec, BfMat const *mat) {
  vec->vtbl->MulInplace(vec, mat);
}

/* Solve the system `mat`*lhs = `vec` inplace, storing the solution of
 * the linear system in `vec`. */
void bfVecSolveInplace(BfVec *vec, BfMat const *mat) {
  vec->vtbl->SolveInplace(vec, mat);
}

/* Compute the reciprocal of `vec` inplace. */
void bfVecRecipInplace(BfVec *vec) {
  vec->vtbl->RecipInplace(vec);
}

/* Find and return the Givens rotation acting on the indices `srcInd`
 * and `elimInd`, and eliminating `vec[elimInd]`. */
BfMat *bfVecGetGivensRotation(BfVec const *vec, BfSize srcInd, BfSize elimInd) {
  return vec->vtbl->GetGivensRotation(vec, srcInd, elimInd);
}

/* Permute the `vec` by `perm`. */
void bfVecPermute(BfVec *vec, BfPerm const *perm) {
  vec->vtbl->Permute(vec, perm);
}

/* Concatenate `vec` together with `otherVec`, returning a new vector. */
BfVec *bfVecConcat(BfVec const *vec, BfVec const *otherVec) {
  return vec->vtbl->Concat(vec, otherVec);
}

void bfVecSave(BfVec const *vec, char const *path) {
  vec->vtbl->Save(vec, path);
}

void bfVecDaxpy(BfVec *vec, BfReal scale, BfVec const *otherVec) {
  vec->vtbl->Daxpy(vec, scale, otherVec);
}

void bfVecDscal(BfVec *vec, BfReal scale) {
  vec->vtbl->Dscal(vec, scale);
}

/** Implementation: Vec */

void bfVecInit(BfVec *vec, BfVecVtable *vtbl, BfSize size) {
  vec->vtbl = vtbl;
  vec->size = size;
  vec->props = BF_VEC_PROPS_NONE;
}

void bfVecDeinit(BfVec *vec) {
  vec->vtbl = NULL;
  vec->size = BF_SIZE_BAD_VALUE;
  vec->props = BF_VEC_PROPS_NONE;
}

BfVec *bfVecFromFile(char const *path, BfSize size, BfDtype dtype) {
  switch (dtype) {
  case BF_DTYPE_COMPLEX:
    return (BfVec *)bfVecComplexFromFile(path, size);
  case BF_DTYPE_REAL:
    return (BfVec *)bfVecRealFromFile(path, size);
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

bool bfVecInstanceOf(BfVec const *vec, BfType type) {
  return bfTypeDerivedFrom(bfVecGetType(vec), type);
}
