#include <bf/vec_real.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat.h>
#include <bf/mat_givens.h>

/** Interface: Vec */

#define INTERFACE BF_INTERFACE_Vec
BF_DEFINE_VTABLE(Vec, VecReal)
#undef INTERFACE

BF_STUB(BfVec *, VecRealCopy, BfVec const *)

void bfVecRealDelete(BfVec **mat) {
  bfVecRealDeinitAndDealloc((BfVecReal **)mat);
}

BfType bfVecRealGetType(BfVec const *vec) {
  (void)vec;
  return BF_TYPE_VEC_REAL;
}

BF_STUB(bool, VecRealInstanceOf, BfVec const *, BfType)
BF_STUB(BfPtr, VecRealGetEltPtr, BfVec *, BfSize)

BfVec *bfVecRealGetSubvecCopy(BfVec const *vec, BfSize i0, BfSize i1) {
  BEGIN_ERROR_HANDLING();

  BfVecReal const *vecReal = NULL;
  BfVecReal *subvec = NULL;

  vecReal = bfVecConstToVecRealConst(vec);
  HANDLE_ERROR();

  subvec = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(subvec, i1 - i0);
  HANDLE_ERROR();

  memcpy(subvec->data, vecReal->data + i0, (i1 - i0)*sizeof(BfReal));

  END_ERROR_HANDLING() {
    bfVecRealDeinitAndDealloc(&subvec);
    subvec = NULL;
  }

  return bfVecRealToVec(subvec);
}

BF_STUB(BfVec *, VecRealGetSubvecView, BfVec *, BfSize, BfSize)

void bfVecRealPrint(BfVec const *vec, FILE *fp) {
  BEGIN_ERROR_HANDLING();

  BfVecReal const *vecReal = bfVecConstToVecRealConst(vec);
  HANDLE_ERROR();

  for (BfSize i = 0; i < vec->size; ++i)
    fprintf(fp, "%g\n", *(vecReal->data + i*vecReal->stride));

  END_ERROR_HANDLING() {}
}

BF_STUB(BfReal, VecRealDist, BfVec const *, BfVec const *)

BfReal bfVecRealNormMax(BfVec const *vec) {
  BEGIN_ERROR_HANDLING();

  BfVecReal const *vecReal = NULL;
  BfReal norm;

  vecReal = bfVecConstToVecRealConst(vec);
  HANDLE_ERROR();

  norm = -INFINITY;
  for (BfSize i = 0; i < vec->size; ++i) {
    BfReal y = *(vecReal->data + i*vecReal->stride);
    BfReal yabs = fabs(y);
    if (yabs > norm) norm = yabs;
  }

  END_ERROR_HANDLING()
    norm = NAN;

  return norm;
}

BF_STUB(void, VecRealScaleByReal, BfVec *, BfReal)
BF_STUB(void, VecRealAddInplace, BfVec *, BfVec const *)
BF_STUB(void, VecRealMulInplace, BfVec *, BfMat const *)
BF_STUB(void, VecRealSolveInplace, BfVec *, BfMat const *)

void bfVecRealRecipInplace(BfVec *vec) {
  BEGIN_ERROR_HANDLING();

  BfVecReal *vecReal = bfVecToVecReal(vec);
  HANDLE_ERROR();

  BfReal *ptr = vecReal->data;
  for (BfSize i = 0; i < vec->size; ++i) {
    BfReal value = *ptr;
    *ptr = 1/value;
  }

  END_ERROR_HANDLING() {}
}

BF_STUB(BfMat *, VecRealGetGivensRotation, BfVec const *, BfSize, BfSize)

/** Upcasting: */

BfVec *bfVecRealToVec(BfVecReal *vecReal) {
  return &vecReal->super;
}

/** Downcasting: */

BfVecReal *bfVecToVecReal(BfVec *vec) {
  if (!bfVecInstanceOf(vec, BF_TYPE_VEC_REAL)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfVecReal *)vec;
  }
}

BfVecReal const *bfVecConstToVecRealConst(BfVec const *vec) {
  if (!bfVecInstanceOf(vec, BF_TYPE_VEC_REAL)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfVecReal const *)vec;
  }
}

/** Implementation: VecReal */

BfVecReal *bfVecRealNew() {
  BEGIN_ERROR_HANDLING();

  BfVecReal *vecReal = malloc(sizeof(BfVecReal));
  if (vecReal == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return vecReal;
}

BfVecReal *bfVecRealFromFile(char const *path, BfSize size) {
  BEGIN_ERROR_HANDLING();

  FILE *fp = fopen(path, "rb");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  BfVecReal *vecReal = bfVecRealNew();
  if (vecReal == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  if (size == BF_SIZE_BAD_VALUE) {
    fseek(fp, 0, SEEK_END);
    BfSize numBytes = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    if (numBytes % sizeof(BfReal) != 0)
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);
    size = numBytes/sizeof(BfReal);
  }

  bfVecRealInit(vecReal, size);
  HANDLE_ERROR();

  fread(vecReal->data, sizeof(BfReal), size, fp);
  if (ferror(fp)) {
    clearerr(fp);
    RAISE_ERROR(BF_ERROR_FILE_ERROR);
  }

  END_ERROR_HANDLING()
    bfVecRealDeinitAndDealloc(&vecReal);

  fclose(fp);

  return vecReal;
}

void bfVecRealInit(BfVecReal *vec, BfSize size) {
  BEGIN_ERROR_HANDLING();

  bfVecInit(&vec->super, &VecVtbl, size);

  vec->stride = 1;

  vec->data = malloc(size*sizeof(BfReal));
  if (vec->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING()
    bfVecDeinit(&vec->super);
}

void bfVecRealInitView(BfVecReal *vecReal, BfSize size, BfSize stride,
                       BfReal *data) {
  BEGIN_ERROR_HANDLING();

  if (data == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  bfVecInit(&vecReal->super, &VecVtbl, size);

  vecReal->super.props |= BF_VEC_PROPS_VIEW;

  vecReal->stride = stride;
  vecReal->data = data;

  END_ERROR_HANDLING()
    bfVecDeinit(&vecReal->super);
}

void bfVecRealDeinit(BfVecReal *vecReal) {
  if (!(vecReal->super.props & BF_VEC_PROPS_VIEW))
    free(vecReal->data);

  vecReal->data = NULL;
}

void bfVecRealDealloc(BfVecReal **vecReal) {
  free(*vecReal);
  *vecReal = NULL;
}

void bfVecRealDeinitAndDealloc(BfVecReal **vecReal) {
  bfVecRealDeinit(*vecReal);
  bfVecRealDealloc(vecReal);
}