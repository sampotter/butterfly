#include <bf/vec_complex.h>

#include <assert.h>
#include <math.h>

#include <bf/blas.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_givens.h>
#include <bf/mem.h>
#include <bf/vec_real.h>
#include <bf/vec_zero.h>

/** Interface: Vec */

static BfVecVtable VEC_VTABLE = {
  .Copy = (__typeof__(&bfVecComplexCopy))bfVecComplexCopy,
  .Delete = (__typeof__(&bfVecComplexDelete))bfVecComplexDelete,
  .GetType = (__typeof__(&bfVecComplexGetType))bfVecComplexGetType,
  .GetSubvecView = (__typeof__(&bfVecComplexGetSubvecView))bfVecComplexGetSubvecView,
  .Print = (__typeof__(&bfVecComplexPrint))bfVecComplexPrint,
  .NormMax = (__typeof__(&bfVecComplexNormMax))bfVecComplexNormMax,
  .AddInplace = (__typeof__(&bfVecComplexAddInplace))bfVecComplexAddInplace,
  .MulInplace = (__typeof__(&bfVecComplexMulInplace))bfVecComplexMulInplace,
  .SolveInplace = (__typeof__(&bfVecComplexSolveInplace))bfVecComplexSolveInplace,
  .GetGivensRotation = (__typeof__(&bfVecComplexGetGivensRotation))bfVecComplexGetGivensRotation,
  .Concat = (__typeof__(&bfVecComplexConcat))bfVecComplexConcat,
  .Save = (__typeof__(&bfVecSave))bfVecComplexSave,
};

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

  BfComplex c = givens->c;
  BfComplex s = givens->s;

  *z0_ptr = conj(c)*z0 + -s*z1;
  *z1_ptr =       s*z0 +  c*z1;
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

  BfComplex c = givens->c;
  BfComplex s = givens->s;

  *z0_ptr =        c*z0 + conj(s)*z1;
  *z1_ptr = -conj(s)*z0 + conj(c)*z1;
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

  BfComplex a = *(vecComplex->data + srcInd*vecComplex->stride);
  BfComplex b = *(vecComplex->data + elimInd*vecComplex->stride);
  BfComplex c;
  BfComplex s;

  // TODO: we should replace this with LAPACK's zlartg, but it appears
  // to have been added pretty recently---like 2021 or so. It isn't
  // available in OpenBLAS. May be best to conditionally compile and
  // use zlartg if it's available and the implementation below
  // otherwise as a fallback.

  /* This method of computing a Givens rotation comes from SciPy's
   * implementation of GMRES:
   *
   *   https://github.com/scipy/scipy/blob/main/scipy/sparse/linalg/_isolve/iterative/GMRESREVCOM.f.src
   *
   * which is a fork of gmresrevcom from "templated.f":
   *
   *   https://people.sc.fsu.edu/~jburkardt/f77_src/templated/templated.f
   *
   * which is a FORTRAN77 implementation of the stuff described in the
   * "Templates for the Solution of Linear Systems" book (although the
   * original F77 code just calls out to drotg, which could be of poor
   * quality). */
  if (cabs(b) == 0) {
    c = 1;
    s = 0;
  } else if (cabs(b) > cabs(a)) {
    BfComplex tmp = -a/b;
    s = 1/sqrt(1 + pow(cabs(tmp), 2));
    c = tmp*s;
  } else {
    BfComplex tmp = -b/a;
    c = 1/sqrt(1 + pow(cabs(tmp), 2));
    s = tmp*c;
  }

  bfMatGivensComplexInit(givens, vec->size, srcInd, elimInd, c, s);

  END_ERROR_HANDLING()
    bfMatGivensComplexDeinitAndDealloc(&givens);

  return bfMatGivensComplexToMat(givens);
}

static BfVec *concat_vecComplex(BfVec const *vec, BfVec const *otherVec) {
  BEGIN_ERROR_HANDLING();

  BfVecComplex const *vecComplex = bfVecConstToVecComplexConst(vec);
  HANDLE_ERROR();

  BfVecComplex const *otherVecComplex = bfVecConstToVecComplexConst(otherVec);
  HANDLE_ERROR();

  BfVecComplex *cat = bfVecComplexNew();
  HANDLE_ERROR();

  bfVecComplexInit(cat, vec->size + otherVec->size);
  HANDLE_ERROR();

  BfComplex *writePtr = cat->data;
  BfComplex *readPtr = NULL;

  readPtr = vecComplex->data;
  for (BfSize i = 0; i < vec->size; ++i) {
    *writePtr = *readPtr;
    writePtr += cat->stride;
    readPtr += vecComplex->stride;
  }

  readPtr = otherVecComplex->data;
  for (BfSize i = 0; i < otherVec->size; ++i) {
    *writePtr = *readPtr;
    writePtr += cat->stride;
    readPtr += otherVecComplex->stride;
  }

  END_ERROR_HANDLING() {}

  return bfVecComplexToVec(cat);
}

static BfVec *concat_vecReal(BfVec const *vec, BfVec const *otherVec) {
  BEGIN_ERROR_HANDLING();

  BfVecComplex const *vecComplex = bfVecConstToVecComplexConst(vec);
  HANDLE_ERROR();

  BfVecReal const *vecReal = bfVecConstToVecRealConst(otherVec);
  HANDLE_ERROR();

  BfVecComplex *cat = bfVecComplexNew();
  HANDLE_ERROR();

  bfVecComplexInit(cat, vec->size + otherVec->size);
  HANDLE_ERROR();

  BfComplex *writePtr = cat->data;

  BfComplex *complexReadPtr = vecComplex->data;
  for (BfSize i = 0; i < vec->size; ++i) {
    *writePtr = *complexReadPtr;
    writePtr += cat->stride;
    complexReadPtr += vecComplex->stride;
  }

  BfReal *realReadPtr = vecReal->data;
  for (BfSize i = 0; i < otherVec->size; ++i) {
    *writePtr = *realReadPtr;
    writePtr += cat->stride;
    realReadPtr += vecReal->stride;
  }

  END_ERROR_HANDLING() {}

  return bfVecComplexToVec(cat);
}

static BfVec *concat_vecZero(BfVec const *vec, BfVec const *otherVec) {
  BEGIN_ERROR_HANDLING();

  BfVecComplex *cat = NULL;

  BfVecComplex const *vecComplex = bfVecConstToVecComplexConst(vec);
  HANDLE_ERROR();

  assert(bfVecGetType(otherVec) == BF_TYPE_VEC_ZERO);

  cat = bfVecComplexNew();
  HANDLE_ERROR();

  bfVecComplexInit(cat, vec->size + otherVec->size);
  HANDLE_ERROR();

  /* First, copy over data from vecComplex: */
  BfComplex *writePtr = cat->data;
  BfComplex *readPtr = vecComplex->data;
  for (BfSize i = 0; i < vec->size; ++i) {
    *writePtr = *readPtr;
    writePtr += cat->stride;
    readPtr += vecComplex->stride;
  }

  /* Pad what's left over with zeros */
  for (BfSize i = 0; i < otherVec->size; ++i) {
    *writePtr = 0;
    writePtr += cat->stride;
  }

  END_ERROR_HANDLING() {}

  return bfVecComplexToVec(cat);
}

BfVec *bfVecComplexConcat(BfVec const *vec, BfVec const *otherVec) {
  switch (bfVecGetType(otherVec)) {
  case BF_TYPE_VEC_COMPLEX:
    return concat_vecComplex(vec, otherVec);
  case BF_TYPE_VEC_REAL:
    return concat_vecReal(vec, otherVec);
  case BF_TYPE_VEC_ZERO:
    return concat_vecZero(vec, otherVec);
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

void bfVecComplexSave(BfVecComplex const *vecComplex, char const *path) {
  BEGIN_ERROR_HANDLING();

  FILE *fp = fopen(path, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  for (BfSize i = 0; i < vecComplex->super.size; ++i) {
    BfComplex const *ptr = vecComplex->data + i*vecComplex->stride;
    fwrite(ptr, sizeof(BfComplex), 1, fp);
  }

  END_ERROR_HANDLING() {
    assert(false);
  }

  fclose(fp);
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

  BfVecComplex *vecComplex = bfMemAlloc(1, sizeof(BfVecComplex));
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

  bfVecInit(&vecComplex->super, &VEC_VTABLE, size);

  vecComplex->stride = 1;

  vecComplex->data = bfMemAlloc(size, sizeof(BfComplex));
  if (vecComplex->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING()
    bfVecDeinit(&vecComplex->super);
}

void bfVecComplexInitView(BfVecComplex *vecComplex, BfSize size, BfSize stride, BfComplex *data) {
  BEGIN_ERROR_HANDLING();

  if (data == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  bfVecInit(&vecComplex->super, &VEC_VTABLE, size);

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
