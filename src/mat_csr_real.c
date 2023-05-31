#include <bf/mat_csr_real.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/vec_real.h>

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .GetView = (__typeof__(&bfMatGetView))bfMatCsrRealGetView,
  .Copy = (__typeof__(&bfMatCsrRealCopy))bfMatCsrRealCopy,
  .Delete = (__typeof__(&bfMatCsrRealDelete))bfMatCsrRealDelete,
  .GetType = (__typeof__(&bfMatCsrRealGetType))bfMatCsrRealGetType,
  .GetNumRows = (__typeof__(&bfMatCsrRealGetNumRows))bfMatCsrRealGetNumRows,
  .GetNumCols = (__typeof__(&bfMatCsrRealGetNumCols))bfMatCsrRealGetNumCols,
  .AddInplace = (__typeof__(&bfMatCsrRealAddInplace))bfMatCsrRealAddInplace,
  .Scale = (__typeof__(&bfMatCsrRealScale))bfMatCsrRealScale,
  .MulVec = (__typeof__(&bfMatCsrRealMulVec))bfMatCsrRealMulVec,
  .IsZero = (__typeof__(&bfMatCsrRealIsZero))bfMatCsrRealIsZero,
};

BfMat *bfMatCsrRealGetView(BfMatCsrReal *matCsrReal) {
  BF_ERROR_BEGIN();

  BfMatCsrReal *matCsrRealView = bfMatCsrRealNew();
  HANDLE_ERROR();

  *matCsrRealView = *matCsrReal;

  BfMat *matView = bfMatCsrRealToMat(matCsrRealView);

  matView->props |= BF_MAT_PROPS_VIEW;

  BF_ERROR_END()
    matView = NULL;

  return matView;
}

BfMat *bfMatCsrRealCopy(BfMat const *mat) {
  BF_ERROR_BEGIN();

  BfMatCsrReal const *matCsrReal = bfMatConstToMatCsrRealConst(mat);
  HANDLE_ERROR();

  BfMatCsrReal *copy = bfMatCsrRealNew();
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);

  bfMatCsrRealInit(copy, m, n, matCsrReal->rowptr, matCsrReal->colind, matCsrReal->data);
  HANDLE_ERROR();

  BF_ERROR_END() {
    bfMatCsrRealDeinitAndDealloc(&copy);
  }

  return bfMatCsrRealToMat(copy);
}

void bfMatCsrRealDelete(BfMat **mat) {
  bfMatCsrRealDeinitAndDealloc((BfMatCsrReal **)mat);
}

BfType bfMatCsrRealGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_CSR_REAL;
}

BfSize bfMatCsrRealGetNumRows(BfMat const *mat) {
  if (bfMatGetType(mat) != BF_TYPE_MAT_CSR_REAL) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return BF_SIZE_BAD_VALUE;
  } else {
    return mat->numRows;
  }
}

BfSize bfMatCsrRealGetNumCols(BfMat const *mat) {
  if (bfMatGetType(mat) != BF_TYPE_MAT_CSR_REAL) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return BF_SIZE_BAD_VALUE;
  } else {
    return mat->numCols;
  }
}

void bfMatCsrRealScale(BfMat *mat, BfComplex scalar) {
  BF_ERROR_BEGIN();

  if (cimag(scalar) != 0)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfReal scalar_ = creal(scalar);

  BfMatCsrReal *matCsrReal = bfMatToMatCsrReal(mat);
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(mat);
  BfSize nnz = matCsrReal->rowptr[m];
  for (BfSize i = 0; i < nnz; ++i)
    matCsrReal->data[i] *= scalar_;

  BF_ERROR_END() {}
}

void bfMatCsrRealAddInplace(BfMat *mat, BfMat const *otherMat) {
  BF_ERROR_BEGIN();

  BfMatCsrReal *matCsrReal = bfMatToMatCsrReal(mat);
  HANDLE_ERROR();

  BfMatCsrReal const *otherMatCsrReal = bfMatConstToMatCsrRealConst(otherMat);
  HANDLE_ERROR();

  if (!bfMatCsrRealHasSameSparsityPattern(matCsrReal, otherMatCsrReal))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize m = bfMatGetNumRows(mat);
  BfSize nnz = matCsrReal->rowptr[m];
  for (BfSize i = 0; i < nnz; ++i)
    matCsrReal->data[i] += otherMatCsrReal->data[i];

  BF_ERROR_END() {}
}

static BfVec *mulVec_vecReal(BfMat const *mat, BfVecReal const *vecReal) {
  BF_ERROR_BEGIN();

  BfMatCsrReal const *matCsrReal = bfMatConstToMatCsrRealConst(mat);
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(mat);

  BfSize n = bfMatGetNumCols(mat);
  if (vecReal->super.size != n)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfVecReal *result = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(result, m);
  HANDLE_ERROR();

  BfReal const *src = vecReal->data;
  BfReal *dst = result->data;
  for (BfSize i = 0; i < m; ++i) {
    *dst = 0;
    for (BfSize j = matCsrReal->rowptr[i]; j < matCsrReal->rowptr[i + 1]; ++j) {
      BfReal elt = *(src + vecReal->stride*matCsrReal->colind[j]);
      *dst += elt*matCsrReal->data[j];
    }
    dst += result->stride;
  }

  BF_ERROR_END() {
    bfVecRealDeinitAndDealloc(&result);
  }

  return bfVecRealToVec(result);
}

BfVec *bfMatCsrRealMulVec(BfMat const *mat, BfVec const *vec) {
  BF_ERROR_BEGIN();

  BfVec *result = NULL;

  switch (bfVecGetType(vec)) {
  case BF_TYPE_VEC_REAL:
    result = mulVec_vecReal(mat, bfVecConstToVecRealConst(vec));
    HANDLE_ERROR();
    break;
  default:
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);
  }

  BF_ERROR_END() {
    bfVecDelete(&result);
  }

  return result;
}

bool bfMatCsrRealIsZero(BfMat const *mat) {
  return bfMatConstToMatCsrRealConst(mat)->rowptr[mat->numRows] == 0;
}

/** Upcasting: */

BfMat *bfMatCsrRealToMat(BfMatCsrReal *matCsrReal) {
  return &matCsrReal->super;
}

BfMat const *bfMatCsrRealConstToMatConst(BfMatCsrReal const *matCsrReal) {
  return &matCsrReal->super;
}

/** Downcasting: */

BfMatCsrReal *bfMatToMatCsrReal(BfMat *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_CSR_REAL)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatCsrReal *)mat;
  }
}

BfMatCsrReal const *bfMatConstToMatCsrRealConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_CSR_REAL)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatCsrReal const *)mat;
  }
}

/** Implementation: MatCsrReal */

BfMatCsrReal *bfMatCsrRealNew() {
  BF_ERROR_BEGIN();

  BfMatCsrReal *mat = bfMemAlloc(1, sizeof(BfMatCsrReal));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BF_ERROR_END() {}

  return mat;
}

void bfMatCsrRealInit(BfMatCsrReal *mat, BfSize numRows, BfSize numCols,
                      BfSize const *rowptr, BfSize const *colind,
                      BfReal const *data) {
  BF_ERROR_BEGIN();

  BfSize nnz = rowptr[numRows];

  bfMatInit(&mat->super, &MAT_VTABLE, numRows, numCols);
  HANDLE_ERROR();

  mat->rowptr = bfMemAlloc(numRows + 1, sizeof(BfSize));
  if (mat->rowptr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  mat->colind = bfMemAlloc(nnz, sizeof(BfSize));
  if (mat->colind == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  mat->data = bfMemAlloc(nnz, sizeof(BfReal));
  if (mat->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  bfMemCopy(rowptr, numRows + 1, sizeof(BfSize), mat->rowptr);
  bfMemCopy(colind, nnz, sizeof(BfSize), mat->colind);
  bfMemCopy(data, nnz, sizeof(BfSize), mat->data);

  BF_ERROR_END()
    bfMatCsrRealDeinit(mat);
}

void bfMatCsrRealDeinit(BfMatCsrReal *mat) {
  if (!(mat->super.props & BF_MAT_PROPS_VIEW)) {
    bfMemFree(mat->rowptr);
    bfMemFree(mat->colind);
    bfMemFree(mat->data);
  }

  mat->rowptr = NULL;
  mat->colind = NULL;
  mat->data = NULL;

  bfMatDeinit(&mat->super);
}

void bfMatCsrRealDealloc(BfMatCsrReal **mat) {
  bfMemFree(*mat);
  *mat = NULL;
}

void bfMatCsrRealDeinitAndDealloc(BfMatCsrReal **mat) {
  bfMatCsrRealDeinit(*mat);
  bfMatCsrRealDealloc(mat);
}

void bfMatCsrRealDump(BfMatCsrReal const *matCsrReal, char const *rowptrPath,
                      char const *colindPath, char const *dataPath) {
  BF_ERROR_BEGIN();

  BfSize numRows = matCsrReal->super.numRows;
  BfSize nnz = matCsrReal->rowptr[numRows];

  FILE *fp;

  /* Save `rowptr` to disk: */

  fp = fopen(rowptrPath, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);
  fwrite(matCsrReal->rowptr, sizeof(BfSize), numRows + 1, fp);
  fclose(fp);

  /* Save `colind` to disk: */

  fp = fopen(colindPath, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);
  fwrite(matCsrReal->colind, sizeof(BfSize), nnz, fp);
  fclose(fp);

  /* Save `data` to disk: */

  fp = fopen(dataPath, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);
  fwrite(matCsrReal->data, sizeof(BfReal), nnz, fp);
  fclose(fp);

  BF_ERROR_END() {
    // TODO: delete all files
    fclose(fp);
  }
}

bool bfMatCsrRealHasSameSparsityPattern(BfMatCsrReal const *matCsrReal, BfMatCsrReal const *otherMatCsrReal) {
  BfMat const *mat = &matCsrReal->super;
  BfMat const *otherMat = &otherMatCsrReal->super;

  BfSize m = bfMatGetNumRows(mat);
  if (m != bfMatGetNumRows(otherMat))
    return false;

  if (bfMatGetNumCols(mat) != bfMatGetNumCols(otherMat))
    return false;

  for (BfSize i = 0; i < m; ++i)
    if (matCsrReal->rowptr[i] != otherMatCsrReal->rowptr[i])
      return false;

  BfSize nnz = matCsrReal->rowptr[m];
  if (nnz != otherMatCsrReal->rowptr[m])
    return false;

  for (BfSize i = 0; i < nnz; ++i)
    if (matCsrReal->colind[i] != otherMatCsrReal->colind[i])
      return false;

  return true;
}
