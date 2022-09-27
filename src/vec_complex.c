#include <bf/vec_complex.h>

#include <math.h>
#include <stdlib.h>

#include <bf/blas.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_givens.h>

/** Interface: Vec */

#define INTERFACE BF_INTERFACE_Vec
BF_DEFINE_VTABLE(Vec, VecComplex)
#undef INTERFACE

BfVec *bfVecComplexCopy(BfVec const *vec) {
  BEGIN_ERROR_HANDLING();

  BfVecComplex const *vecComplex = NULL;
  BfVecComplex *copy = NULL;

  vecComplex = bfVecConstToVecComplexConst(vec);
  HANDLE_ERROR();

  copy = bfVecComplexNew();
  HANDLE_ERROR();

  bfVecComplexInit(copy, vec->size);
  HANDLE_ERROR();

  BfComplex const *inPtr = vecComplex->data;
  BfComplex *outPtr = copy->data;
  for (BfSize i = 0; i < vec->size; ++i) {
    *outPtr = *inPtr;
    inPtr += vecComplex->stride;
    outPtr += copy->stride;
  }

  END_ERROR_HANDLING()
    bfVecComplexDeinitAndDealloc(&copy);

  return bfVecComplexToVec(copy);
}

void bfVecComplexDelete(BfVec **mat) {
  bfVecComplexDeinitAndDealloc((BfVecComplex **)mat);
}

BfType bfVecComplexGetType(BfVec const *vec) {
  (void)vec;
  return BF_TYPE_VEC_COMPLEX;
}

BF_STUB(bool, VecComplexInstanceOf, BfVec const *, BfType)
BF_STUB(BfPtr, VecComplexGetEltPtr, BfVec *, BfSize)
BF_STUB(BfVec *, VecComplexGetSubvecCopy, BfVec const *, BfSize, BfSize)

BfVec *bfVecComplexGetSubvecView(BfVec *vec, BfSize i0, BfSize i1) {
  BEGIN_ERROR_HANDLING();

  BfVecComplex *vecComplex = NULL;
  BfVecComplex *view = NULL;

  vecComplex = bfVecToVecComplex(vec);
  HANDLE_ERROR();

  view = bfVecComplexNew();
  HANDLE_ERROR();

  BfSize size = i1 - i0;
  BfSize stride = vecComplex->stride;
  BfComplex *data = vecComplex->data + i0*vecComplex->stride;

  bfVecComplexInitView(view, size, stride, data);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfVecComplexDeinitAndDealloc(&view);

  return bfVecComplexToVec(view);
}

void bfVecComplexPrint(BfVec const *vec, FILE *fp) {
  BEGIN_ERROR_HANDLING();

  BfVecComplex const *vecComplex = bfVecConstToVecComplexConst(vec);
  HANDLE_ERROR();

  for (BfSize i = 0; i < vec->size; ++i) {
    BfComplex z = *(vecComplex->data + i*vecComplex->stride);
    fprintf(fp, "%g + i*%g\n", creal(z), cimag(z));
  }

  END_ERROR_HANDLING() {}
}

BF_STUB(BfReal, VecComplexDist, BfVec const *, BfVec const *)

BfReal bfVecComplexNormMax(BfVec const *vec) {
  BEGIN_ERROR_HANDLING();

  BfVecComplex const *vecComplex = NULL;
  BfReal norm;

  vecComplex = bfVecConstToVecComplexConst(vec);
  HANDLE_ERROR();

  norm = -INFINITY;
  for (BfSize i = 0; i < vec->size; ++i) {
    BfComplex z = *(vecComplex->data + i*vecComplex->stride);
    BfReal zabs = cabs(z);
    if (zabs > norm) norm = zabs;
  }

  END_ERROR_HANDLING()
    norm = NAN;

  return norm;
}

BF_STUB(void, VecComplexScaleByReal, BfVec *, BfReal)

void bfVecComplexAddInplace(BfVec *vec, BfVec const *otherVec) {
  BEGIN_ERROR_HANDLING();

  BfVecComplex *vecComplex = NULL;
  BfVecComplex const *otherVecComplex = NULL;

  if (vec->size != otherVec->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  vecComplex = bfVecToVecComplex(vec);
  HANDLE_ERROR();

  otherVecComplex = bfVecConstToVecComplexConst(otherVec);
  HANDLE_ERROR();

  BfComplex *outPtr = vecComplex->data;
  BfComplex const *inPtr = otherVecComplex->data;
  for (BfSize i = 0; i < vec->size; ++i) {
    *outPtr += *inPtr;
    outPtr += vecComplex->stride;
    inPtr += otherVecComplex->stride;
  }

  END_ERROR_HANDLING() {}
}

static void
mulInplace_givensComplex(BfVecComplex *vecComplex, BfMatGivensComplex const *givens) {
  BfComplex *z0_ptr = vecComplex->data + givens->srcInd*vecComplex->stride;
  BfComplex *z1_ptr = vecComplex->data + givens->elimInd*vecComplex->stride;

  BfComplex z0 = *z0_ptr;
  BfComplex z1 = *z1_ptr;

  BfReal c = givens->c;
  BfComplex s = givens->s;

  *z0_ptr = c*z0 - s*z1;
  *z1_ptr = s*z0 + c*z1;
}

void bfVecComplexMulInplace(BfVec *vec, BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfVecComplex *vecComplex = bfVecToVecComplex(vec);
  HANDLE_ERROR();

  if (vec->size != bfMatGetNumCols(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  switch (bfMatGetType(mat)) {
  case BF_TYPE_MAT_GIVENS_COMPLEX: {
    BfMatGivensComplex const *givens = bfMatConstToMatGivensComplexConst(mat);
    HANDLE_ERROR();
    mulInplace_givensComplex(vecComplex, givens);
    break;
  } default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING() {}
}

static void
solveInplace_givensComplex(BfVecComplex *vecComplex, BfMatGivensComplex const *givens) {
  BfComplex *z0_ptr = vecComplex->data + givens->srcInd*vecComplex->stride;
  BfComplex *z1_ptr = vecComplex->data + givens->elimInd*vecComplex->stride;

  BfComplex z0 = *z0_ptr;
  BfComplex z1 = *z1_ptr;

  BfReal c = givens->c;
  BfComplex s = givens->s;

  *z0_ptr =        c*z0 + conj(s)*z1;
  *z1_ptr = -conj(s)*z0 +       c*z1;
}

void bfVecComplexSolveInplace(BfVec *vec, BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfVecComplex *vecComplex = bfVecToVecComplex(vec);
  HANDLE_ERROR();

  if (vec->size != bfMatGetNumRows(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  switch (bfMatGetType(mat)) {
  case BF_TYPE_MAT_GIVENS_COMPLEX: {
    BfMatGivensComplex const *givens = bfMatConstToMatGivensComplexConst(mat);
    HANDLE_ERROR();
    solveInplace_givensComplex(vecComplex, givens);
    break;
  } default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING() {}
}

BF_STUB(void, VecComplexRecipInplace, BfVec *)

BfMat *bfVecComplexGetGivensRotation(BfVec const *vec, BfSize srcInd, BfSize elimInd) {
  BEGIN_ERROR_HANDLING();

  BfVecComplex const *vecComplex = NULL;
  BfMatGivensComplex *givens = NULL;

  if (srcInd == elimInd)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  vecComplex = bfVecConstToVecComplexConst(vec);
  HANDLE_ERROR();

  if (srcInd >= vec->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (elimInd >= vec->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  givens = bfMatGivensComplexNew();
  HANDLE_ERROR();


// def zrotg(a, b):
//     if a == 0:
//         c, s = 0, 1
//     else:
//         tmp = b/a
//         c = 1/np.sqrt(1 + abs(tmp)**2)
//         s = tmp*c
//     return c, np.conj(s)

  BfComplex a = *(vecComplex->data + srcInd*vecComplex->stride);
  BfComplex b = *(vecComplex->data + elimInd*vecComplex->stride);
  BfReal c;
  BfComplex s;

  /* This way of calculating the Givens is from Flatiron's BIE3D. */
  if (a == 0) {
    c = 0;
    s = 1;
  } else {
    BfComplex tmp = b/a;
    c = 1/sqrt(1 + pow(cabs(tmp), 2));
    s = conj(tmp*c); // conjugate here because of how we define a
                     // Givens rotation matrix
  }

  bfMatGivensComplexInit(givens, vec->size, srcInd, elimInd, c, s);

  END_ERROR_HANDLING()
    bfMatGivensComplexDeinitAndDealloc(&givens);

  return bfMatGivensComplexToMat(givens);
}

/** Upcasting: */

BfVec *bfVecComplexToVec(BfVecComplex *vecComplex) {
  return &vecComplex->super;
}

/** Downcasting: */

BfVecComplex *bfVecToVecComplex(BfVec *vec) {
  if (!bfVecInstanceOf(vec, BF_TYPE_VEC_COMPLEX)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfVecComplex *)vec;
  }
}

BfVecComplex const *bfVecConstToVecComplexConst(BfVec const *vec) {
  if (!bfVecInstanceOf(vec, BF_TYPE_VEC_COMPLEX)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfVecComplex const *)vec;
  }
}

/** Implementation: VecComplex */

BfVecComplex *bfVecComplexNew() {
  BEGIN_ERROR_HANDLING();

  BfVecComplex *vecComplex = malloc(sizeof(BfVecComplex));
  if (vecComplex == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return vecComplex;
}

BfVecComplex *bfVecComplexFromFile(char const *path, BfSize size) {
  BEGIN_ERROR_HANDLING();

  FILE *fp = fopen(path, "rb");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  BfVecComplex *vecComplex = bfVecComplexNew();
  if (vecComplex == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  if (size == BF_SIZE_BAD_VALUE) {
    fseek(fp, 0, SEEK_END);
    BfSize numBytes = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    if (numBytes % sizeof(BfComplex) != 0)
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);
    size = numBytes/sizeof(BfComplex);
  }

  bfVecComplexInit(vecComplex, size);
  HANDLE_ERROR();

  fread(vecComplex->data, sizeof(BfComplex), size, fp);
  if (ferror(fp)) {
    clearerr(fp);
    RAISE_ERROR(BF_ERROR_FILE_ERROR);
  }

  END_ERROR_HANDLING()
    bfVecComplexDeinitAndDealloc(&vecComplex);

  fclose(fp);

  return vecComplex;
}

void bfVecComplexInit(BfVecComplex *vecComplex, BfSize size) {
  BEGIN_ERROR_HANDLING();

  bfVecInit(&vecComplex->super, &VecVtbl, size);

  vecComplex->stride = 1;

  vecComplex->data = malloc(size*sizeof(BfComplex));
  if (vecComplex->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING()
    bfVecDeinit(&vecComplex->super);
}

void bfVecComplexInitView(BfVecComplex *vecComplex, BfSize size, BfSize stride, BfComplex *data) {
  BEGIN_ERROR_HANDLING();

  if (data == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  bfVecInit(&vecComplex->super, &VecVtbl, size);

  vecComplex->super.props |= BF_VEC_PROPS_VIEW;

  vecComplex->stride = stride;
  vecComplex->data = data;

  END_ERROR_HANDLING()
    bfVecDeinit(&vecComplex->super);
}

void bfVecComplexDeinit(BfVecComplex *vecComplex) {
  if (!(vecComplex->super.props & BF_VEC_PROPS_VIEW))
    free(vecComplex->data);

  vecComplex->data = NULL;
}

void bfVecComplexDealloc(BfVecComplex **vecComplex) {
  free(*vecComplex);
  *vecComplex = NULL;
}

void bfVecComplexDeinitAndDealloc(BfVecComplex **vecComplex) {
  bfVecComplexDeinit(*vecComplex);
  bfVecComplexDealloc(vecComplex);
}
