#include <bf/mat_dense_complex.h>

#include <assert.h>
#include <math.h>
#include <openblas/lapacke.h>
#include <stdlib.h>
#include <string.h>

#include <bf/blas.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_diag_real.h>

/** Static functions: */

static enum CBLAS_TRANSPOSE getCblasTranspose(BfMatDenseComplex const *mat) {
  BfMat const *super = bfMatDenseComplexConstToMatConst(mat);
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

/** Interface: Mat */

#define INTERFACE BF_INTERFACE_Mat
BF_DEFINE_VTABLE(Mat, MatDenseComplex)
#undef INTERFACE

BfMat *bfMatDenseComplexEmptyLike(BfMat const *mat, BfSize numRows, BfSize numCols) {
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

  return bfMatDenseComplexToMat(result);
}

BfMat *bfMatDenseComplexZerosLike(BfMat const *mat, BfSize numRows, BfSize numCols) {
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

  return bfMatDenseComplexToMat(zeros);
}

BfMatType bfMatDenseComplexGetType(BfMat const *mat) {
  (void)mat;
  return BF_MAT_TYPE_DENSE_COMPLEX;
}

bool bfMatDenseComplexInstanceOf(BfMat const *mat, BfMatType matType) {
  return bfMatTypeDerivedFrom(bfMatGetType(mat), matType);
}

BfSize bfMatDenseComplexNumBytes(BfMat const *mat) {
  (void)mat;
  assert(false);
  return BF_SIZE_BAD_VALUE;
}

void bfMatDenseComplexSave(BfMat const *mat, char const *path) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex const *matDenseComplex = bfMatConstToMatDenseComplexConst(mat);

  FILE *fp = fopen(path, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  BfSize numElts = mat->numRows*mat->numCols;

  fwrite(matDenseComplex->data, sizeof(BfComplex), numElts, fp);
  /* TODO: error-handling */

  END_ERROR_HANDLING() {}

  fclose(fp);
}

BF_STUB(void, MatDenseComplexPrint, FILE *, BfMat const *)

BfSize bfMatDenseComplexGetNumRows(BfMat const *mat) {
  return bfMatIsTransposed(mat) ? mat->numCols : mat->numRows;
}

BfSize bfMatDenseComplexGetNumCols(BfMat const *mat) {
  return bfMatIsTransposed(mat) ? mat->numRows : mat->numCols;
}

BfMat *bfMatDenseComplexGetRowRange(BfMat *mat, BfSize i0, BfSize i1) {
  BfSize numRows = mat->numRows;

  assert(i0 < i1);
  assert(i1 <= numRows);
  assert(!bfMatIsTransposed(mat)); // TODO: implement

  BEGIN_ERROR_HANDLING();

  BfMat *matView = bfMatGetView(mat);

  BfMatDenseComplex *submat = bfMatToMatDenseComplex(matView);
  HANDLE_ERROR();

  if (i1 - i0 != numRows) {
    bfMatDenseComplexToMat(submat)->numRows = i1 - i0;
    submat->data += submat ->rowStride*i0;
  }

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDealloc(&submat);

  return bfMatDenseComplexToMat(submat);
}

BfMat *bfMatDenseComplexGetColRange(BfMat *mat, BfSize j0, BfSize j1) {
  BEGIN_ERROR_HANDLING();

  if (bfMatIsTransposed(mat))
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfSize numCols = mat->numCols;
  if (!(j0 <= j1 && j1 <= numCols))
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  BfMat *matView = bfMatGetView(mat);

  BfMatDenseComplex *matDenseComplexView = bfMatToMatDenseComplex(matView);
  HANDLE_ERROR();

  if (j1 - j0 < numCols) {
    matDenseComplexView->super.numCols = j1 - j0;
    matDenseComplexView->data += matDenseComplexView->colStride*j0;
  }

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDealloc(&matDenseComplexView);

  return bfMatDenseComplexToMat(matDenseComplexView);
}

void bfMatDenseComplexSetRowRange(BfMat *mat, BfSize i0, BfSize i1, BfMat const *rows) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *matRows =
    bfMatToMatDenseComplex(bfMatDenseComplexGetRowRange(mat, i0, i1));
  HANDLE_ERROR();

  switch (bfMatGetType(rows)) {
  case BF_MAT_TYPE_DENSE_COMPLEX:
    bfMatDenseComplexCopy(matRows, bfMatConstToMatDenseComplexConst(rows));
    HANDLE_ERROR();
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING() {}
}

BF_STUB(BfMat *, MatDenseComplexRowDists, BfMat const *, BfMat const *)
BF_STUB(void, MatDenseComplexScaleCols, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatDenseComplexSumCols, BfMat const *)

static void bfMatDenseComplexDenseComplexAddInplace(
  BfMatDenseComplex *op1, BfMatDenseComplex const *op2)
{
  BEGIN_ERROR_HANDLING();

  BfMat *mat = bfMatDenseComplexToMat(op1);
  BfMat const *otherMat = bfMatDenseComplexConstToMatConst(op2);

  BfSize numRows = bfMatDenseComplexGetNumRows(mat);
  if (numRows != bfMatDenseComplexGetNumRows(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numCols = bfMatDenseComplexGetNumCols(mat);
  if (numCols != bfMatDenseComplexGetNumCols(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  for (BfSize i = 0; i < numRows*numCols; ++i)
    op1->data[i] += op2->data[i];

  END_ERROR_HANDLING() {}
}

void bfMatDenseComplexAddInplace(BfMat *mat, BfMat const *otherMat) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *matDenseComplex = bfMatToMatDenseComplex(mat);

  switch (bfMatGetType(otherMat)) {
  case BF_MAT_TYPE_DENSE_COMPLEX:
    bfMatDenseComplexDenseComplexAddInplace(
      matDenseComplex, bfMatConstToMatDenseComplexConst(otherMat));
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING() {}
}

BF_STUB(void, MatDenseComplexAddDiag, BfMat *, BfMat const *)

BfMat *bfMatDenseComplexMul(BfMat const *op1, BfMat const *op2) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex const *matDenseComplex = bfMatConstToMatDenseComplexConst(op1);

  BfMat *result;

  switch (bfMatGetType(op2)) {
  case BF_MAT_TYPE_DENSE_COMPLEX:
    result = bfMatDenseComplexToMat(
      bfMatDenseComplexDenseComplexMul(
        matDenseComplex, bfMatConstToMatDenseComplexConst(op2)));
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

  BfMat const *super1 = bfMatDenseComplexConstToMatConst(op1);
  BfMat const *super2 = bfMatDenseComplexConstToMatConst(op2);

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

BF_STUB(void, MatDenseComplexMulInplace, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatDenseComplexSolve, BfMat const *, BfMat const *)

BfMat *bfMatDenseComplexLstSq(BfMat const *lhs, BfMat const *rhs) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex const *matDenseComplex = bfMatConstToMatDenseComplexConst(lhs);

  BfMat *result;

  switch (bfMatGetType(rhs)) {
  case BF_MAT_TYPE_DENSE_COMPLEX:
    result = bfMatDenseComplexToMat(
      bfMatDenseComplexDenseComplexLstSq(
        matDenseComplex, bfMatConstToMatDenseComplexConst(rhs)));
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
  BfMat const *lhsSuper = bfMatDenseComplexConstToMatConst(lhs);
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

  /* TODO: the casting below got kind of out of hand... clean it up */

  /* get subblocks of truncated SVD */

  BfMat *UkH = bfMatDenseComplexGetColRange(bfMatDenseComplexToMat(U), 0, k);
  bfMatConjTrans(UkH);

  BfMatDiagReal *Sk = bfMatDiagRealGetDiagBlock(S, 0, k);

  BfMat *Vk = bfMatGetRowRange(bfMatDenseComplexToMat(VH), 0, k);
  bfMatConjTrans(Vk);

  /* solve least squares problem */

  BfMat *tmp1 = bfMatMul(UkH, bfMatDenseComplexConstToMatConst(rhs));
  HANDLE_ERROR();

  BfMat *tmp2 = bfMatDenseComplexToMat(
    bfMatDiagRealDenseComplexSolve(Sk, bfMatToMatDenseComplex(tmp1)));
  HANDLE_ERROR();

  BfMat *result = bfMatMul(Vk, tmp2);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfMatDenseComplexDelete(&result);

  bfMatDenseComplexDeinitAndDealloc(&U);
  bfMatDiagRealDeinitAndDealloc(&S);
  bfMatDenseComplexDeinitAndDealloc(&VH);

  bfMatDenseComplexDelete(&UkH);
  bfMatDiagRealDeinitAndDealloc(&Sk);
  bfMatDenseComplexDelete(&Vk);

  bfMatDelete(&tmp1);
  bfMatDelete(&tmp2);

  return bfMatToMatDenseComplex(result);
}

BfMatDenseComplex *bfMatDenseComplexNewView(BfMatDenseComplex *matDenseComplex) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *matDenseComplexView = malloc(sizeof(BfMatDenseComplex));
  if (matDenseComplexView == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  *matDenseComplexView = *matDenseComplex;

  bfMatDenseComplexToMat(matDenseComplexView)->props |= BF_MAT_PROPS_VIEW;

  END_ERROR_HANDLING() {
    free(matDenseComplexView);
    matDenseComplexView = NULL;
  }

  return matDenseComplexView;
}

BfSize bfMatDenseComplexGetRowStride(BfMatDenseComplex const *mat) {
  return bfMatIsTransposed(bfMatDenseComplexConstToMatConst(mat)) ?
    mat->colStride :
    mat->rowStride;
}

BfSize bfMatDenseComplexGetColStride(BfMatDenseComplex const *mat) {
  return bfMatIsTransposed(bfMatDenseComplexConstToMatConst(mat)) ?
    mat->rowStride :
    mat->colStride;
}

void bfMatDenseComplexCopy(BfMatDenseComplex *dst,
                           BfMatDenseComplex const *src) {
  BEGIN_ERROR_HANDLING();

  BfMat *mat = bfMatDenseComplexToMat(dst);
  BfMat const *matOther = bfMatDenseComplexConstToMatConst(src);

  BfSize dstRows = bfMatDenseComplexGetNumRows(mat);
  if (dstRows != bfMatDenseComplexGetNumRows(matOther))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize dstCols = bfMatDenseComplexGetNumCols(mat);
  if (dstCols != bfMatDenseComplexGetNumCols(matOther))
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

void bfMatDenseComplexSvd(BfMatDenseComplex const *mat, BfMatDenseComplex *U,
                          BfMatDiagReal *S, BfMatDenseComplex *VH) {
  BEGIN_ERROR_HANDLING();

  BfMat const *super = bfMatDenseComplexConstToMatConst(mat);

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

  bfMatDenseComplexToMat(U)->props |= BF_MAT_PROPS_ORTHO;
  bfMatDenseComplexToMat(VH)->props |= BF_MAT_PROPS_ORTHO;

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

void bf_zmat_add_diag(BfMatDenseComplex *mat, BfReal value) {
  BfSize m = mat->super.numRows;
  BfSize n = mat->super.numCols;
  BfSize p = m < n ? m : n;
  for (BfSize i = 0; i < p; ++i)
    *(mat->data + i*mat->rowStride + i*mat->colStride) += value;
}

BfMatDenseComplex *bf_zmat_solve(BfMatDenseComplex *A, BfMatDenseComplex *B) {
  BfSize m = A->super.numRows;
  if (B->super.numRows != m) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return NULL;
  }

  BfSize n = A->super.numCols;
  if (m != n) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return NULL;
  }

  BfSize p = B->super.numCols; /* == number of RHSs */

  BEGIN_ERROR_HANDLING();

  int *ipiv = malloc(m*sizeof(int));

  lapack_int error = LAPACKE_zgesv(
    LAPACK_ROW_MAJOR, m, p, A->data, A->rowStride,
    ipiv, B->data, B->rowStride);

  if (error == LAPACK_WORK_MEMORY_ERROR)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  if (error != 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR)

  END_ERROR_HANDLING()
    free(ipiv);

  return B;
}

/** Upcasting: */

BfMat *bfMatDenseComplexToMat(BfMatDenseComplex *matDenseComplex) {
  return &matDenseComplex->super;
}

BfMat const *bfMatDenseComplexConstToMatConst(BfMatDenseComplex const *matDenseComplex) {
  return &matDenseComplex->super;
}

/** Downcasting: */

BfMatDenseComplex *bfMatToMatDenseComplex(BfMat *mat) {
  if (!bfMatInstanceOf(mat, BF_MAT_TYPE_DENSE_COMPLEX)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatDenseComplex *)mat;
  }
}

BfMatDenseComplex const *bfMatConstToMatDenseComplexConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_MAT_TYPE_DENSE_COMPLEX)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatDenseComplex const *)mat;
  }
}

/** Implementation: MatDenseComplex */

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

BfMatDenseComplex *bfMatDenseComplexFromFile(char const *path, BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *matDenseComplex = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(matDenseComplex, numRows, numCols);
  HANDLE_ERROR();

  FILE *fp = fopen(path, "r");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  BfSize size = numRows*numCols;
  if (fread(matDenseComplex->data, sizeof(BfComplex), size, fp) != size)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDealloc(&matDenseComplex);

  return matDenseComplex;
}

void bfMatDenseComplexInit(BfMatDenseComplex *mat,
                           BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  bfMatInit(&mat->super, &MatVtbl, numRows, numCols);
  HANDLE_ERROR();

  mat->rowStride = numCols;
  mat->colStride = 1;

  mat->data = malloc(numRows*numCols*sizeof(BfComplex));
  if (mat->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING()
    bfMatDeinit(&mat->super);
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

BF_STUB(BfMat *, MatDenseComplexGetView, BfMat *)

void bfMatDenseComplexDelete(BfMat **mat) {
  bfMatDenseComplexDeinitAndDealloc((BfMatDenseComplex **)mat);
}
