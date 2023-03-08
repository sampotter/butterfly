#include <bf/vec_real.h>

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat.h>
#include <bf/mat_givens.h>
#include <bf/rand.h>

/** Interface: Vec */

static BfVecVtable VEC_VTABLE = {
  .Copy = (__typeof__(&bfVecRealCopy))bfVecRealCopy,
  .Delete = (__typeof__(&bfVecRealDelete))bfVecRealDelete,
  .GetType = (__typeof__(&bfVecRealGetType))bfVecRealGetType,
  .GetEltPtr = (__typeof__(&bfVecRealGetEltPtr))bfVecRealGetEltPtr,
  .GetSubvecCopy = (__typeof__(&bfVecRealGetSubvecCopy))bfVecRealGetSubvecCopy,
  .GetSubvecView = (__typeof__(&bfVecGetSubvecView))bfVecRealGetSubvecView,
  .GetSubvecViewConst = (__typeof__(&bfVecGetSubvecViewConst))bfVecRealGetSubvecViewConst,
  .SetRange = (__typeof__(&bfVecSetRange))bfVecRealSetRange,
  .Print = (__typeof__(&bfVecRealPrint))bfVecRealPrint,
  .DistMax = (__typeof__(&bfVecDistMax))bfVecRealDistMax,
  .NormMax = (__typeof__(&bfVecRealNormMax))bfVecRealNormMax,
  .RecipInplace = (__typeof__(&bfVecRealRecipInplace))bfVecRealRecipInplace,
  .AddInplace = (__typeof__(&bfVecAddInplace))bfVecRealAddInplace,
  .Permute = (__typeof__(&bfVecRealPermute))bfVecRealPermute,
  .Concat = (__typeof__(&bfVecRealConcat))bfVecRealConcat,
};

BfVec *bfVecRealCopy(BfVec const *vec) {
  BEGIN_ERROR_HANDLING();

  BfVecReal const *vecReal = NULL;
  BfVecReal *copy = NULL;

  vecReal = bfVecConstToVecRealConst(vec);
  HANDLE_ERROR();

  copy = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(copy, vec->size);
  HANDLE_ERROR();

  BfReal const *inPtr = vecReal->data;
  BfReal *outPtr = copy->data;
  for (BfSize i = 0; i < vec->size; ++i) {
    *outPtr = *inPtr;
    inPtr += vecReal->stride;
    outPtr += copy->stride;
  }

  END_ERROR_HANDLING()
    bfVecRealDeinitAndDealloc(&copy);

  return bfVecRealToVec(copy);
}

void bfVecRealDelete(BfVec **mat) {
  bfVecRealDeinitAndDealloc((BfVecReal **)mat);
}

BfType bfVecRealGetType(BfVec const *vec) {
  (void)vec;
  return BF_TYPE_VEC_REAL;
}

BfPtr bfVecRealGetEltPtr(BfVec *vec, BfSize i) {
  BEGIN_ERROR_HANDLING();

  BfVecReal *vecReal = NULL;
  BfPtr ptr = NULL;

  vecReal = bfVecToVecReal(vec);
  HANDLE_ERROR();

  if (i >= vec->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  ptr = vecReal->data + i*vecReal->stride;

  END_ERROR_HANDLING() {}

  return ptr;
}

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

BfVecReal *bfVecRealGetSubvecView(BfVecReal *vecReal, BfSize i0, BfSize i1) {
  BEGIN_ERROR_HANDLING();

  BfVecReal *vecRealView = NULL;

  vecRealView = bfVecRealNew();
  HANDLE_ERROR();

  BfSize size = i1 - i0;
  BfSize stride = vecReal->stride;
  BfReal *data = vecReal->data + i0*vecReal->stride;

  bfVecRealInitView(vecRealView, size, stride, data);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfVecRealDeinitAndDealloc(&vecRealView);

  return vecRealView;
}

BfVecReal const *bfVecRealGetSubvecViewConst(BfVecReal const *vecReal, BfSize i0, BfSize i1) {
  BEGIN_ERROR_HANDLING();

  BfVecReal *vecRealView = NULL;

  vecRealView = bfVecRealNew();
  HANDLE_ERROR();

  BfSize size = i1 - i0;
  BfSize stride = vecReal->stride;
  BfReal *data = vecReal->data + i0*vecReal->stride;

  bfVecRealInitView(vecRealView, size, stride, data);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfVecRealDeinitAndDealloc(&vecRealView);

  return vecRealView;
}

void bfVecRealSetRange(BfVecReal *vecReal, BfSize i0, BfSize i1,
                       BfVec const *otherVec) {
  BEGIN_ERROR_HANDLING();

  if (i0 > i1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (bfVecGetType(otherVec) != BF_TYPE_VEC_REAL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfVec *vec = bfVecRealToVec(vecReal);

  BfSize n = vec->size;
  if (i1 > n || i0 > n)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfVecReal const *otherVecReal = bfVecConstToVecRealConst(otherVec);

  for (BfSize i = i0; i < i1; ++i)
    *(vecReal->data + i*vecReal->stride)
      = *(otherVecReal->data + (i - i0)*otherVecReal->stride);

  END_ERROR_HANDLING() {
    assert(false);
  }
}

void bfVecRealPrint(BfVec const *vec, FILE *fp) {
  BEGIN_ERROR_HANDLING();

  BfVecReal const *vecReal = bfVecConstToVecRealConst(vec);
  HANDLE_ERROR();

  for (BfSize i = 0; i < vec->size; ++i)
    fprintf(fp, "%g\n", *(vecReal->data + i*vecReal->stride));

  END_ERROR_HANDLING() {}
}

BfReal bfVecRealDistMax(BfVecReal const *vecReal, BfVec const *otherVec) {
  BEGIN_ERROR_HANDLING();

  if (bfVecGetType(otherVec) != BF_TYPE_VEC_REAL)
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfVec const *vec = bfVecRealConstToVecConst(vecReal);
  BfVecReal const *otherVecReal = bfVecConstToVecRealConst(otherVec);

  if (vec->size != otherVec->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfReal dist = -INFINITY;
  for (BfSize i = 0; i < vec->size; ++i) {
    BfReal x = *(vecReal->data + i*vecReal->stride);
    BfReal y = *(otherVecReal->data + i*otherVecReal->stride);
    BfReal abs_x_minus_y = fabs(x - y);
    if (abs_x_minus_y > dist) dist = abs_x_minus_y;
  }

  END_ERROR_HANDLING() {}

  return dist;
}

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

void bfVecRealAddInplace(BfVecReal *vecReal, BfVec const *otherVec) {
  BEGIN_ERROR_HANDLING();

  BfVec *vec = bfVecRealToVec(vecReal);

  BfSize n = vec->size;
  if (vec->size != otherVec->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (bfVecGetType(otherVec) != BF_TYPE_VEC_REAL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfVecReal const *otherVecReal = bfVecConstToVecRealConst(otherVec);

  for (BfSize i = 0; i < n; ++i)
    *(vecReal->data + i*vecReal->stride)
      += *(otherVecReal->data + i*otherVecReal->stride);

  END_ERROR_HANDLING() {
    assert(false);
  }
}

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

void bfVecRealPermute(BfVec *vec, BfPerm const *perm) {
  BEGIN_ERROR_HANDLING();

  BfVec *vecCopy = bfVecCopy(vec);
  HANDLE_ERROR();

  BfVecReal *vecReal = bfVecToVecReal(vec);
  HANDLE_ERROR();

  BfVecReal *vecRealCopy = bfVecToVecReal(vecCopy);
  HANDLE_ERROR();

  for (BfSize i = 0; i < vec->size; ++i) {
    BfReal const *inPtr = vecRealCopy->data + i*vecRealCopy->stride;
    BfReal *outPtr = vecReal->data + perm->index[i]*vecReal->stride;
    *outPtr = *inPtr;
  }

  END_ERROR_HANDLING() {}

  bfVecDelete(&vecCopy);
}

BfVec *bfVecRealConcat(BfVec const *vec, BfVec const *otherVec) {
  switch (bfVecGetType(otherVec)) {
    (void)vec;
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

/** Upcasting: */

BfVec *bfVecRealToVec(BfVecReal *vecReal) {
  return &vecReal->super;
}

BfVec const *bfVecRealConstToVecConst(BfVecReal const *vecReal) {
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

BfVecReal *bfVecRealNewWithValue(BfSize n, BfReal value) {
  BEGIN_ERROR_HANDLING();

  BfVecReal *vecReal = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(vecReal, n);
  HANDLE_ERROR();

  for (BfSize i = 0; i < n; ++i)
    *(vecReal->data + i*vecReal->stride) = value;

  END_ERROR_HANDLING() {
    assert(false);
  }

  return vecReal;
}

BfVecReal *bfVecRealNewRandn(BfSize n) {
  BEGIN_ERROR_HANDLING();

  BfVecReal *vecReal = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(vecReal, n);
  HANDLE_ERROR();

  bfRealRandn(n, vecReal->data);

  END_ERROR_HANDLING() {
    assert(false);
  }

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

  bfVecInit(&vec->super, &VEC_VTABLE, size);

  vec->stride = 1;

  vec->data = malloc(size*sizeof(BfReal));
  if (vec->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING()
    bfVecDeinit(&vec->super);
}

void bfVecRealInitFrom(BfVecReal *vec, BfSize size, BfSize stride, BfReal const *data) {
  BEGIN_ERROR_HANDLING();

  bfVecInit(&vec->super, &VEC_VTABLE, size);

  vec->stride = 1;

  vec->data = malloc(size*sizeof(BfReal));
  if (vec->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfReal *writePtr = vec->data;
  BfReal const *readPtr = data;
  for (BfSize i = 0; i < size; ++i) {
    *writePtr = *readPtr;
    writePtr += vec->stride;
    readPtr += stride;
  }

  END_ERROR_HANDLING()
    bfVecDeinit(&vec->super);
}

void bfVecRealInitView(BfVecReal *vecReal, BfSize size, BfSize stride, BfReal *data) {
  BEGIN_ERROR_HANDLING();

  if (data == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  bfVecInit(&vecReal->super, &VEC_VTABLE, size);

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

void bfVecRealDump(BfVecReal const *vecReal, char const *path) {
  BEGIN_ERROR_HANDLING();

  FILE *fp = fopen(path, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  fwrite(vecReal->data, sizeof(BfReal), vecReal->super.size, fp);

  END_ERROR_HANDLING() {}

  fclose(fp);
}
