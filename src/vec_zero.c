#include <bf/vec_zero.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/vec_complex.h>

/** Interface: Vec */

static BfVecVtable VEC_VTABLE = {
  .Delete = (__typeof__(&bfVecZeroDelete))bfVecZeroDelete,
  .GetType = (__typeof__(&bfVecZeroGetType))bfVecZeroGetType,
  .Concat = (__typeof__(&bfVecZeroConcat))bfVecZeroConcat,
};

void bfVecZeroDelete(BfVec **mat) {
  bfVecZeroDeinitAndDealloc((BfVecZero **)mat);
}

BfType bfVecZeroGetType(BfVec const *vec) {
  (void)vec;
  return BF_TYPE_VEC_ZERO;
}

static BfVec *concat_vecComplex(BfVec const *vec, BfVec const *otherVec) {
  BEGIN_ERROR_HANDLING();

  BfVecComplex *cat = NULL;

  assert(bfVecGetType(vec) == BF_TYPE_VEC_ZERO);

  BfVecComplex const *vecComplex = bfVecConstToVecComplexConst(otherVec);
  HANDLE_ERROR();

  cat = bfVecComplexNew();
  HANDLE_ERROR();

  bfVecComplexInit(cat, vec->size + otherVec->size);
  HANDLE_ERROR();

  BfComplex *writePtr = cat->data;

  /* First, pad with zeros to match the size of `vec` */
  for (BfSize i = 0; i < vec->size; ++i) {
    *writePtr = 0;
    writePtr += cat->stride;
  }

  /* Next, copy over data from vecComplex: */
  BfComplex *readPtr = vecComplex->data;
  for (BfSize i = 0; i < otherVec->size; ++i) {
    *writePtr = *readPtr;
    writePtr += cat->stride;
    readPtr += vecComplex->stride;
  }

  END_ERROR_HANDLING() {}

  return bfVecComplexToVec(cat);
}

BfVec *concat_vecZero(BfVec const *vec, BfVec const *otherVec) {
  BEGIN_ERROR_HANDLING();

  assert(bfVecGetType(vec) == BF_TYPE_VEC_ZERO);
  assert(bfVecGetType(otherVec) == BF_TYPE_VEC_ZERO);

  BfVecZero *cat = bfVecZeroNew();
  HANDLE_ERROR();

  bfVecZeroInit(cat, vec->size + otherVec->size);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {}

  return bfVecZeroToVec(cat);
}

BfVec *bfVecZeroConcat(BfVec const *vec, BfVec const *otherVec) {
  switch (bfVecGetType(otherVec)) {
  case BF_TYPE_VEC_COMPLEX:
    return concat_vecComplex(vec, otherVec);
  case BF_TYPE_VEC_ZERO:
    return concat_vecZero(vec, otherVec);
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

/** Upcasting: */

BfVec *bfVecZeroToVec(BfVecZero *vecZero) {
  return &vecZero->super;
}

/** Downcasting: */

BfVecZero const *bfVecConstToVecZeroConst(BfVec const *vec) {
  if (!bfVecInstanceOf(vec, BF_TYPE_VEC_ZERO)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfVecZero const *)vec;
  }
}

/** Implementation: VecZero */

BfVecZero *bfVecZeroNew() {
  BEGIN_ERROR_HANDLING();

  BfVecZero *vecZero = malloc(sizeof(BfVecZero));
  if (vecZero == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return vecZero;
}

void bfVecZeroInit(BfVecZero *vecZero, BfSize size) {
  BEGIN_ERROR_HANDLING();

  bfVecInit(&vecZero->super, &VEC_VTABLE, size);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfVecDeinit(&vecZero->super);
}

void bfVecZeroDeinit(BfVecZero *vecZero) {
  (void)vecZero;
}

void bfVecZeroDealloc(BfVecZero **vecZero) {
  free(*vecZero);
  *vecZero = NULL;
}

void bfVecZeroDeinitAndDealloc(BfVecZero **vecZero) {
  bfVecZeroDeinit(*vecZero);
  bfVecZeroDealloc(vecZero);
}
