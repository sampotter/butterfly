#include <bf/mat_dense_complex.h>

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <bf/blas.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_diag_real.h>

BF_DEFINE_MAT_VTABLE(MatDenseComplex);

static enum CBLAS_TRANSPOSE getCblasTranspose(BfMatDenseComplex const *mat) {
  BfMat const *super = bfMatDenseComplexGetMatConstPtr(mat);

  if (super->props & (BF_MAT_PROPS_TRANS | BF_MAT_PROPS_CONJ))
    return CblasConjTrans;
  else if (super->props & BF_MAT_PROPS_TRANS)
    return CblasTrans;
  else
    return CblasNoTrans;
}

static BfSize getLeadingDimension(BfMatDenseComplex const *mat) {
  return mat->rowStride;
}

BfMatDenseComplex *bfMatDenseComplexNew() {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *mat = malloc(sizeof(BfMatDenseComplex));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return mat;
}

BfMatDenseComplex *bfMatDenseComplexZeros(BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *zeros = NULL;

  zeros = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(zeros, numRows, numCols);
  HANDLE_ERROR();

  for (BfSize i = 0; i < numRows*numCols; ++i)
    zeros->data[i++] = 0;

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDealloc(&zeros);

  return zeros;
}

void bfMatDenseComplexInit(BfMatDenseComplex *mat,
                           BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  bfMatInit(&mat->super, &matVtbl, numRows, numCols);
  HANDLE_ERROR();

  mat->rowStride = numCols;
  mat->colStride = 1;

  mat->data = malloc(numRows*numCols*sizeof(BfComplex));
  if (mat->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING()
    bfMatDeinit(&mat->super);
}

void bfMatDenseComplexDelete(BfMatDenseComplex **mat) {
  bfMatDenseComplexDeinitAndDealloc(mat);
}

BfMatDenseComplex *
bfMatDenseComplexEmptyLike(BfMatDenseComplex const *mat,
                           BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *result = NULL;

  if (numRows == BF_SIZE_BAD_VALUE)
    numRows = bfMatDenseComplexGetNumRows(mat);

  if (numCols == BF_SIZE_BAD_VALUE)
    numCols = bfMatDenseComplexGetNumCols(mat);

  result = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(result, numRows, numCols);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDealloc(&result);

  return result;
}

BfMatDenseComplex *
bfMatDenseComplexZerosLike(BfMatDenseComplex const *mat,
                           BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *zeros = NULL;

  if (numRows == BF_SIZE_BAD_VALUE)
    numRows = bfMatDenseComplexGetNumRows(mat);

  if (numCols == BF_SIZE_BAD_VALUE)
    numCols = bfMatDenseComplexGetNumCols(mat);

  zeros = bfMatDenseComplexZeros(numRows, numCols);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDealloc(&zeros);

  return zeros;
}

void bfMatDenseComplexDeinit(BfMatDenseComplex *mat) {
  if (!(mat->super.props & BF_MAT_PROPS_VIEW))
    free(mat->data);

  mat->data = NULL;
}

void bfMatDenseComplexDealloc(BfMatDenseComplex **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatDenseComplexDeinitAndDealloc(BfMatDenseComplex **mat) {
  bfMatDenseComplexDeinit(*mat);
  bfMatDenseComplexDealloc(mat);
}

BfMatDenseComplex *bfMatDenseComplexFromMatPtr(BfMat *mat) {
  BEGIN_ERROR_HANDLING();

  if (bfMatGetType(mat) != BF_MAT_TYPE_DENSE_COMPLEX)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  END_ERROR_HANDLING() {}

  return (BfMatDenseComplex *)mat;
}

BfMatDenseComplex const *bfMatDenseComplexFromMatConstPtr(BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  if (bfMatGetType(mat) != BF_MAT_TYPE_DENSE_COMPLEX)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  END_ERROR_HANDLING() {}

  return (BfMatDenseComplex const *)mat;
}

BfMatType bfMatDenseComplexGetType(BfMatDenseComplex const *mat) {
  return BF_MAT_TYPE_DENSE_COMPLEX;
}

bool bfMatDenseComplexInstanceOf(BfMatDenseComplex const *mat, BfMatType matType) {
  BfMat const *parent = bfMatDenseComplexGetMatConstPtr(mat);
  return bfMatTypeDerivedFrom(bfMatGetType(parent), matType);
}

BfSize bfMatDenseComplexNumBytes(BfMatDenseComplex const *mat) {
  (void)mat;
  assert(false);
  return BF_SIZE_BAD_VALUE;
}

void bfMatDenseComplexSave(BfMatDenseComplex const *mat, char const *path) {
  BEGIN_ERROR_HANDLING();

  FILE *fp = fopen(path, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  BfMat const *super = bfMatDenseComplexGetMatConstPtr(mat);

  BfSize numElts = super->numRows*super->numCols;

  fwrite(mat->data, sizeof(BfComplex), numElts, fp);
  /* TODO: error-handling */

  END_ERROR_HANDLING() {}

  fclose(fp);
}

BfSize bfMatDenseComplexGetNumRows(BfMatDenseComplex const *mat) {
  BfMat const *super = bfMatDenseComplexGetMatConstPtr(mat);
  return bfMatIsTransposed(super) ? mat->super.numCols : mat->super.numRows;
}

BfSize bfMatDenseComplexGetNumCols(BfMatDenseComplex const *mat) {
  BfMat const *super = bfMatDenseComplexGetMatConstPtr(mat);
  return bfMatIsTransposed(super) ? mat->super.numRows : mat->super.numCols;
}

BfMatDenseComplex *bfMatDenseComplexNewView(BfMatDenseComplex *mat) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *view = malloc(sizeof(BfMatDenseComplex));
  if (view == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  *view = *mat;

  bfMatDenseComplexGetMatPtr(view)->props |= BF_MAT_PROPS_VIEW;

  END_ERROR_HANDLING() {
    free(view);
    view = NULL;
  }

  return view;
}

BfMat *bfMatDenseComplexGetMatPtr(BfMatDenseComplex *mat) {
  return &mat->super;
}

BfMat const *bfMatDenseComplexGetMatConstPtr(BfMatDenseComplex const *mat) {
  return &mat->super;
}

BfSize bfMatDenseComplexGetRowStride(BfMatDenseComplex const *mat) {
  return bfMatIsTransposed(bfMatDenseComplexGetMatConstPtr(mat)) ?
    mat->colStride :
    mat->rowStride;
}

BfSize bfMatDenseComplexGetColStride(BfMatDenseComplex const *mat) {
  return bfMatIsTransposed(bfMatDenseComplexGetMatConstPtr(mat)) ?
    mat->rowStride :
    mat->colStride;
}

void bfMatDenseComplexCopy(BfMatDenseComplex *dst,
                           BfMatDenseComplex const *src) {
  BEGIN_ERROR_HANDLING();

  BfSize dstRows = bfMatDenseComplexGetNumRows(dst);
  if (dstRows != bfMatDenseComplexGetNumRows(src))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize dstCols = bfMatDenseComplexGetNumCols(dst);
  if (dstCols != bfMatDenseComplexGetNumCols(src))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize dstRowStride = bfMatDenseComplexGetRowStride(dst);
  BfSize dstColStride = bfMatDenseComplexGetColStride(dst);

  BfSize srcRowStride = bfMatDenseComplexGetRowStride(src);
  BfSize srcColStride = bfMatDenseComplexGetColStride(src);

  for (BfSize i = 0; i < dstRows; ++i) {
    BfComplex *dstData = dst->data + i*dstRowStride;
    BfComplex const *srcData = src->data + i*srcRowStride;
    for (BfSize j = 0; j < dstCols; ++j)
      *(dstData + j*dstColStride) = *(srcData + j*srcColStride);
  }

  END_ERROR_HANDLING() {}
}

BfMatDenseComplex *
bfMatDenseComplexGetColRange(BfMatDenseComplex *mat, BfSize j0, BfSize j1) {
  BfMat *super = bfMatDenseComplexGetMatPtr(mat);
  assert(!bfMatIsTransposed(super)); // TODO: implement

  BfSize numCols = super->numCols;

  assert(j0 < j1);
  assert(j1 <= numCols);

  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *submat = bfMatDenseComplexNewView(mat);
  HANDLE_ERROR();

  if (j1 - j0 == numCols)
    return submat;

  submat->super.numCols = j1 - j0;
  submat->data += submat->colStride*j0;

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDealloc(&submat);

  return submat;
}

BfMatDenseComplex *
bfMatDenseComplexGetRowRange(BfMatDenseComplex *mat, BfSize i0, BfSize i1) {
  BfMat *super = bfMatDenseComplexGetMatPtr(mat);
  BfSize numRows = super->numRows;

  assert(i0 < i1);
  assert(i1 <= numRows);
  assert(!bfMatIsTransposed(super)); // TODO: implement

  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *submat = bfMatDenseComplexNewView(mat);
  HANDLE_ERROR();

  if (i1 - i0 == numRows)
    return submat;

  bfMatDenseComplexGetMatPtr(submat)->numRows = i1 - i0;
  submat->data += submat ->rowStride*i0;

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDealloc(&submat);

  return submat;
}

void
bfMatDenseComplexSetRowRange(BfMatDenseComplex *mat,
                             BfSize i0, BfSize i1, BfMat const *rows) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *matRows = bfMatDenseComplexGetRowRange(mat, i0, i1);
  HANDLE_ERROR();

  switch (bfMatGetType(rows)) {
  case BF_MAT_TYPE_DENSE_COMPLEX:
    bfMatDenseComplexCopy(matRows, bfMatDenseComplexFromMatConstPtr(rows));
    HANDLE_ERROR();
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING() {}
}

void bfMatDenseComplexDenseComplexAddInplace(BfMatDenseComplex *op1,
                                             BfMatDenseComplex const *op2) {
  BEGIN_ERROR_HANDLING();

  BfSize numRows = bfMatDenseComplexGetNumRows(op1);
  if (numRows != bfMatDenseComplexGetNumRows(op2))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numCols = bfMatDenseComplexGetNumCols(op1);
  if (numCols != bfMatDenseComplexGetNumCols(op2))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  for (BfSize i = 0; i < numRows*numCols; ++i)
    op1->data[i] += op2->data[i];

  END_ERROR_HANDLING() {}
}

void bfMatDenseComplexAddInplace(BfMatDenseComplex *op1, BfMat const *op2) {
  BEGIN_ERROR_HANDLING();

  switch (bfMatGetType(op2)) {
  case BF_MAT_TYPE_DENSE_COMPLEX:
    bfMatDenseComplexDenseComplexAddInplace(
      op1, bfMatDenseComplexFromMatConstPtr(op2));
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING() {}
}

BfMat *bfMatDenseComplexMul(BfMatDenseComplex const *op1, BfMat const *op2) {
  BEGIN_ERROR_HANDLING();

  BfMat *result;

  switch (bfMatGetType(op2)) {
  case BF_MAT_TYPE_DENSE_COMPLEX:
    result = bfMatDenseComplexGetMatPtr(
      bfMatDenseComplexDenseComplexMul(
        op1, bfMatDenseComplexFromMatConstPtr(op2)));
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING()
    bfMatDelete(&result);

  return result;
}

BfMatDenseComplex *
bfMatDenseComplexDenseComplexMul(BfMatDenseComplex const *op1,
                                 BfMatDenseComplex const *op2)
{
  BEGIN_ERROR_HANDLING();

  enum CBLAS_TRANSPOSE transa = getCblasTranspose(op1);
  enum CBLAS_TRANSPOSE transb = getCblasTranspose(op2);

  BfMat const *super1 = bfMatDenseComplexGetMatConstPtr(op1);
  BfMat const *super2 = bfMatDenseComplexGetMatConstPtr(op2);

  BfMatDenseComplex *result = NULL;

  BfSize k = bfMatGetNumCols(super1);
  if (k != bfMatGetNumRows(super2))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  result = bfMatDenseComplexNew();
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(super1);
  BfSize n = bfMatGetNumCols(super2);
  bfMatDenseComplexInit(result, m, n);
  HANDLE_ERROR();

  BfComplex alpha = 1, beta = 0;

  BfSize lda = getLeadingDimension(op1);
  BfSize ldb = getLeadingDimension(op2);
  BfSize ldc = getLeadingDimension(result);

  assert(m > 0 && n > 0 && k > 0);

  if (transa == CblasNoTrans) {
    assert(lda >= k);
  } else {
    assert(transa == CblasTrans || transa == CblasConjTrans);
    assert(lda >= m);
  }

  if (transb == CblasNoTrans) {
    assert(ldb >= n);
  } else {
    assert(transb == CblasTrans || transb == CblasConjTrans);
    assert(ldb >= k);
  }

  assert(ldc >= n);

  cblas_zgemm(CblasRowMajor, transa, transb, m, n, k,
              &alpha, op1->data, lda, op2->data, ldb, &beta, result->data, ldc);

  /* TODO: handle cblas errors  */

  END_ERROR_HANDLING() {
    if (result != NULL)
      bfMatDenseComplexDeinitAndDealloc(&result);
  }

  return result;
}

BfMat *bfMatDenseComplexLstSq(BfMatDenseComplex const *lhs, BfMat const *rhs) {
  BEGIN_ERROR_HANDLING();

  BfMat *result;

  switch (bfMatGetType(rhs)) {
  case BF_MAT_TYPE_DENSE_COMPLEX:
    result = bfMatDenseComplexGetMatPtr(
      bfMatDenseComplexDenseComplexLstSq(
        lhs, bfMatDenseComplexFromMatConstPtr(rhs)));
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING()
    bfMatDelete(&result);

  return result;
}

BfMatDenseComplex *
bfMatDenseComplexDenseComplexLstSq(BfMatDenseComplex const *lhs,
                                   BfMatDenseComplex const *rhs)
{
  BfMat const *lhsSuper = bfMatDenseComplexGetMatConstPtr(lhs);
  assert(!bfMatIsTransposed(lhsSuper));

  BfSize m = lhsSuper->numRows;
  BfSize n = lhsSuper->numCols;
  BfSize p = m < n ? m : n;

  BEGIN_ERROR_HANDLING();

  /* compute SVD of A */

  BfMatDenseComplex *U = bfMatDenseComplexNew();
  bfMatDenseComplexInit(U, m, p);

  BfMatDiagReal *S = bfMatDiagRealNew();
  bfMatDiagRealInit(S, p, p);

  BfMatDenseComplex *VH = bfMatDenseComplexNew();
  bfMatDenseComplexInit(VH, p, n);

  bfMatDenseComplexSvd(lhs, U, S, VH);
  HANDLE_ERROR();

  /* compute tolerance and compute number of terms to
   * retain in pseudoinverse */

  BfReal const atol = BF_EPS_MACH;
  BfReal const rtol = (m > n ? m : n)*BF_EPS_MACH;

  BfReal const *sigma = S->data;

  BfReal tol = rtol*sigma[0] + atol;

  BfSize k;
  for (k = 0; k < p; ++k)
    if (sigma[k] < tol)
      break;

  /* get subblocks of truncated SVD */

  BfMatDenseComplex *UkH = bfMatDenseComplexGetColRange(U, 0, k);
  bfMatConjTrans(bfMatDenseComplexGetMatPtr(UkH));

  BfMatDiagReal *Sk = bfMatDiagRealGetDiagBlock(S, 0, k);

  BfMatDenseComplex *Vk = bfMatDenseComplexGetRowRange(VH, 0, k);
  bfMatConjTrans(bfMatDenseComplexGetMatPtr(Vk));

  /* solve least squares problem */

  BfMatDenseComplex *tmp1 = bfMatDenseComplexDenseComplexMul(UkH, rhs);
  HANDLE_ERROR();

  BfMatDenseComplex *tmp2 = bfMatDiagRealDenseComplexSolve(Sk, tmp1);
  HANDLE_ERROR();

  BfMatDenseComplex *result = bfMatDenseComplexDenseComplexMul(Vk, tmp2);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDealloc(&result);

  bfMatDenseComplexDeinitAndDealloc(&U);
  bfMatDiagRealDeinitAndDealloc(&S);
  bfMatDenseComplexDeinitAndDealloc(&VH);

  bfMatDenseComplexDeinitAndDealloc(&UkH);
  bfMatDiagRealDeinitAndDealloc(&Sk);
  bfMatDenseComplexDeinitAndDealloc(&Vk);

  bfMatDenseComplexDeinitAndDealloc(&tmp1);
  bfMatDenseComplexDeinitAndDealloc(&tmp2);

  return result;
}

void bfMatDenseComplexSvd(BfMatDenseComplex const *mat, BfMatDenseComplex *U,
                          BfMatDiagReal *S, BfMatDenseComplex *VH) {
  BEGIN_ERROR_HANDLING();

  BfMat const *super = bfMatDenseComplexGetMatConstPtr(mat);

  BfSize m = super->numRows;
  BfSize n = super->numCols;

  BfReal *superb = NULL;

  /* zgesvd will overwrite A, so allocate space for a copy */
  BfComplex *dataCopy = malloc(m*n*sizeof(BfComplex));
  if (dataCopy == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* copy contents of A */
  memcpy(dataCopy, mat->data, m*n*sizeof(BfComplex));
  /* TODO: error-handling */

  /* output array which contains information about superdiagonal
   * elements which didn't converge
   *
   * more info here: tinyurl.com/2p8f5ev3 */
  superb = malloc((((m < n) ? m : n) - 1)*sizeof(BfReal));
  if (superb == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* compute the SVD */
  lapack_int info = LAPACKE_zgesvd(
    LAPACK_ROW_MAJOR, 'S', 'S', m, n, dataCopy, m, S->data, U->data, m,
    VH->data, n, superb);

  /* check for invalid arguments */
  if (info < 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* check for errors */
  if (info > 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  END_ERROR_HANDLING() {}

  bfMatDenseComplexGetMatPtr(U)->props |= BF_MAT_PROPS_ORTHO;
  bfMatDenseComplexGetMatPtr(VH)->props |= BF_MAT_PROPS_ORTHO;

  free(dataCopy);
  free(superb);
}

bool bfMatDenseComplexIsFinite(BfMatDenseComplex const *mat) {
  for (BfSize i = 0; i < mat->super.numRows; ++i) {
    BfComplex *rowPtr = mat->data + i*mat->rowStride;
    for (BfSize j = 0; j < mat->super.numCols; ++j) {
      BfComplex elt = *(rowPtr + j*mat->colStride);
      if (!isfinite(creal(elt)) || !isfinite(cimag(elt)))
        return false;
    }
  }
  return true;
}
