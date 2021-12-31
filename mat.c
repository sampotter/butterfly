#include "mat.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <cblas.h>
#include <lapacke.h>

#include "error_macros.h"
#include "rand.h"

BfMat bfGetUninitializedMat() {
  return (BfMat) {
    .dtype = BF_DTYPE_VOID,
    .props = BF_MAT_PROP_NONE,
    .numRows = BF_SIZE_BAD_VALUE,
    .numCols = BF_SIZE_BAD_VALUE,
    .rowStride = BF_SIZE_BAD_VALUE,
    .colStride = BF_SIZE_BAD_VALUE,
    .data = NULL
  };
}

bool bfMatIsInitialized(BfMat const *A) {
  return A->data != NULL;
}

BfSize bfMatSize(BfMat const *A)
{
  assert(A->dtype != BF_DTYPE_MAT);

  return A->numRows*A->numCols;
}

void bfMatNumBytes(BfMat const *A, BfSize *nbytes)
{
  // TODO: implement
  assert(A->dtype != BF_DTYPE_MAT);

  BfSize dtype_size = bfDtypeSize(A->dtype);

  BfSize m = A->numRows, n = A->numCols;

  BfSize num_entries = A->props & BF_MAT_PROP_DIAGONAL ?
    (m < n ? m : n) :
    m*n;

  *nbytes = dtype_size*num_entries;
}

BfSize bfMatNumRows(BfMat const *A) {
  return bfMatIsTransposed(A) ? A->numCols : A->numRows;
}

BfSize bfMatNumCols(BfMat const *A) {
  return bfMatIsTransposed(A) ? A->numRows : A->numCols;
}

bool bfMatIsAligned(BfMat const *A) {
  BfSize dtype_size = bfDtypeSize(A->dtype);

  return ((BfSize)A->data % dtype_size == 0)
      && (A->rowStride % dtype_size == 0)
      && (A->colStride % dtype_size == 0);
}

BfSize bfMatRowStride(BfMat const *A) {
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  return bfMatIsTransposed(A) ? A->colStride : A->rowStride;
}

BfSize bfMatColStride(BfMat const *A) {
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  return bfMatIsTransposed(A) ? A->rowStride : A->colStride;
}

void bfFreeMat(BfMat *A)
{
  assert(A->dtype != BF_DTYPE_MAT);

  free(A->data);
}

void
bfInitEmptyMat(BfMat *A, enum BfDtypes dtype, enum BfMatProps props,
               BfSize numRows, BfSize numCols)
{
  BEGIN_ERROR_HANDLING();

  A->dtype = dtype;
  A->props = props;
  A->numRows = numRows;
  A->numCols = numCols;

  BfSize dtype_size = bfDtypeSize(A->dtype);

  /* when we set the row stride, if this is a diagonal matrix, we just
   * set the row stride to zero, which makes indexing simpler since we
   * just store each diagonal entry contiguously in memory */
  A->rowStride = A->props & BF_MAT_PROP_DIAGONAL ? 0 : dtype_size*numCols;

  A->colStride = dtype_size;

  A->data = malloc(dtype_size*numRows*numCols);
  if (A->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {
    free(A->data);
  }
}

void
bfMatZeros(BfMat *A, enum BfDtypes dtype, BfSize numRows, BfSize numCols)
{
  BEGIN_ERROR_HANDLING();

  assert(dtype == BF_DTYPE_REAL || dtype == BF_DTYPE_COMPLEX);

  bfInitEmptyMat(A, dtype, BF_MAT_PROP_NONE, numRows, numCols);
  HANDLE_ERROR();

  if (dtype == BF_DTYPE_REAL) {
    BfByte *row_ptr = A->data;
    for (BfSize i = 0; i < numRows; ++i) {
      BfByte *col_ptr = row_ptr;
      for (BfSize j = 0; j < numCols; ++j) {
        *(BfReal *)col_ptr = 0;
        col_ptr += A->colStride;
      }
      row_ptr += A->rowStride;
    }
  }

  if (dtype == BF_DTYPE_COMPLEX) {
    BfByte *row_ptr = A->data;
    for (BfSize i = 0; i < numRows; ++i) {
      BfByte *col_ptr = row_ptr;
      for (BfSize j = 0; j < numCols; ++j) {
        *(BfComplex *)col_ptr = 0;
        col_ptr += A->colStride;
      }
      row_ptr += A->rowStride;
    }
  }

  END_ERROR_HANDLING() {
    bfFreeMat(A);
  }
}

void
bfSaveMat(BfMat const *A, char const *path)
{
  assert(A->dtype != BF_DTYPE_MAT);

  FILE *fp = fopen(path, "w");
  if (fp == NULL) {
    bfSetError(BF_ERROR_RUNTIME_ERROR);
    goto cleanup;
  }

  BfSize dtype_size = bfDtypeSize(A->dtype);
  BfSize m = A->numRows, n = A->numCols;
  BfSize num_entries = A->props & BF_MAT_PROP_DIAGONAL ?
    (m < n ? m : n) :
    m*n;

  fwrite(A->data, dtype_size, num_entries, fp);
  /* TODO: error-handling */

cleanup:
  fclose(fp);
}

void
bfFillMatRandn(BfMat *A)
{
  assert(A->dtype != BF_DTYPE_MAT);

  if (A->dtype == BF_DTYPE_REAL) {
    bfRandn(bfMatSize(A), A->data);
    return;
  }

  if (A->dtype == BF_DTYPE_COMPLEX) {
    bfRandn(2*bfMatSize(A), A->data);
    return;
  }

  if (!bfDtypeIsValid(A->dtype))
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
}

static bool validIndex(BfMat const *A, BfSize i, BfSize j) {
  return i < A->numRows && j < A->numCols;
}

static BfSize offset(BfMat const *A, BfSize i, BfSize j) {
  if (A->props & BF_MAT_PROP_DIAGONAL)
    assert(i == j);

  return i*A->rowStride + j*A->colStride;
}

static void
getDiagonalMatElt(BfMat const *A, BfSize i, BfSize j, BfPtr ptr) {
  assert(A->dtype != BF_DTYPE_MAT);

  if (i == j) {
    BfSize size = bfDtypeSize(A->dtype);
    memcpy(ptr, (BfByte *)A->data + i*size, size);
  } else {
    // TODO: we should make sure to return a zero matrix with a
    // conformable shape here
    bfGetDtypeZero(A->dtype, ptr);
  }
}

void
bfGetMatElt(BfMat const *A, BfSize i, BfSize j, BfPtr ptr) {
  assert(A->dtype != BF_DTYPE_MAT);

  if (!validIndex(A, i, j)) {
    bfSetError(BF_ERROR_OUT_OF_RANGE);
    return;
  }

  if (A->props & BF_MAT_PROP_DIAGONAL) {
    getDiagonalMatElt(A, i, j, ptr);
  } else {
    BfSize dtype_size = bfDtypeSize(A->dtype);
    memcpy(ptr, (BfByte *)A->data + offset(A, i, j), dtype_size);
  }
}

static void
getDiagonalMatEltPtr(BfMat const *A, BfSize i, BfSize j, BfPtr *ptr) {
  if (i != j) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  BfSize size = bfDtypeSize(A->dtype);

  *ptr = (BfByte *)A->data + i*size;
}

void bfGetMatEltPtr(BfMat const *A, BfSize i, BfSize j, BfPtr *ptr) {
  if (!validIndex(A, i, j)) {
    bfSetError(BF_ERROR_OUT_OF_RANGE);
    return;
  }

  if (A->props & BF_MAT_PROP_DIAGONAL) {
    getDiagonalMatEltPtr(A, i, j, ptr);
    return;
  }

  *ptr = (BfByte *)A->data + offset(A, i, j);
}

BfVec bfGetMatRow(BfMat const *A, BfSize i) {
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  return (BfVec) {
    .dtype = A->dtype,
    .props = BF_VEC_PROP_VIEW,
    .size = bfMatNumCols(A),
    .stride = bfMatColStride(A),
    .data = (BfByte *)A->data + i*bfMatRowStride(A)
  };
}

void bfGetRowPtr(BfMat const *A, BfSize i, BfPtr *ptr)
{
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  if (i >= bfMatNumRows(A)) {
    bfSetError(BF_ERROR_OUT_OF_RANGE);
    return;
  }

  *ptr = (BfByte *)A->data + i*bfMatRowStride(A);
}

void bfSetMatRow(BfMat *A, BfSize i, void const *data)
{
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  BfSize m = A->numRows, n = A->numCols;

  if (i >= m) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  BfSize nbytes = bfDtypeSize(A->dtype);

  BfSize row_nbytes = nbytes*n;

  memcpy((BfByte *)A->data + row_nbytes*i, data, row_nbytes);
  /* TODO: error-handling */
}

void
bfCopyMatRow(BfMat const *A, BfSize i, BfMat *B, BfSize j)
{
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  BfSize m = bfMatNumCols(A);

  if (m != bfMatNumCols(B)) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  if (i >= bfMatNumRows(A)) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  if (j >= bfMatNumRows(B)) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  if (A->dtype != B->dtype) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  BfByte *A_row_ptr, *B_row_ptr;
  bfGetRowPtr(A, i, (BfPtr *)&A_row_ptr);
  bfGetRowPtr(B, j, (BfPtr *)&B_row_ptr);

  BfSize A_col_stride = bfMatColStride(A);
  BfSize B_col_stride = bfMatColStride(B);

  // TODO: this sucks
  if (A->dtype == BF_DTYPE_COMPLEX) {
    bool shouldConj = (A->props & BF_MAT_PROP_CONJ_TRANS)
                    ^ (B->props & BF_MAT_PROP_CONJ_TRANS);

    for (BfSize k = 0; k < m; ++k) {
      *(BfComplex *)B_row_ptr = shouldConj ?
        conj(*(BfComplex *)A_row_ptr) : *(BfComplex *)A_row_ptr;

      A_row_ptr += A_col_stride;
      B_row_ptr += B_col_stride;
    }
  } else {
    assert(false); // TODO: not implemented
  }
}

BfMat bfGetMatRowRange(BfMat const *A, BfSize i0, BfSize i1) {
  assert(i0 < i1);
  assert(i1 <= A->numRows);

  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  assert(!bfMatIsTransposed(A));

  BfMat A_rows = *A;
  if (i1 - i0 == A->numRows)
    return A_rows;

  A_rows.props |= BF_MAT_PROP_VIEW;

  if (A_rows.props & BF_MAT_PROP_UNITARY) {
    A_rows.props ^= BF_MAT_PROP_UNITARY;
    A_rows.props ^= BF_MAT_PROP_SEMI_UNITARY;
  }

  A_rows.numRows = i1 - i0;
  A_rows.data = (BfByte *)A_rows.data + A->rowStride*i0;

  return A_rows;
}

BfMat bfGetMatColRange(BfMat const *A, BfSize j0, BfSize j1) {
  assert(j0 < j1);
  assert(j1 <= A->numCols);

  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  assert(!bfMatIsTransposed(A));

  BfMat A_cols = *A;

  A_cols.props |= BF_MAT_PROP_VIEW;
  if (j1 - j0 == A->numCols)
    return A_cols;

  if (A_cols.props & BF_MAT_PROP_UNITARY) {
    A_cols.props ^= BF_MAT_PROP_UNITARY;
    A_cols.props ^= BF_MAT_PROP_SEMI_UNITARY;
  }

  A_cols.numCols = j1 - j0;
  A_cols.data = (BfByte *)A_cols.data + A->colStride*j0;

  return A_cols;
}

BfMat bfGetMatContSubblock(BfMat const *A, BfSize i0, BfSize i1,
                           BfSize j0, BfSize j1) {
  assert(A->dtype != BF_DTYPE_MAT);

  // TODO: if i0 != j0, we need to return a block matrix with zero
  // subblocks except for a single diagonal subblock... this is
  // tedious so just disallow this for now
  if (A->dtype & BF_MAT_PROP_DIAGONAL)
    assert(i0 == j0);

  BfMat A_subblock = *A;

  A_subblock.props |= BF_MAT_PROP_VIEW;

  A_subblock.numRows = i1 - i0;
  A_subblock.numCols = j1 - j0;

  A_subblock.data = (BfByte *)A_subblock.data + offset(A, i0, j0);

  return A_subblock;
}

bool bfMatIsTransposed(BfMat const *A) {
  return A->props & (BF_MAT_PROP_TRANS | BF_MAT_PROP_CONJ_TRANS);
}

BfMat bfConjTrans(BfMat const *A) {
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  BfMat AH = *A;

  /* make AH a view of A */
  AH.props |= BF_MAT_PROP_VIEW;

  /* make AH the Hermitian transpose of A */
  AH.props ^= BF_MAT_PROP_CONJ_TRANS;

  return AH;
}

static enum CBLAS_TRANSPOSE getCblasTranspose(BfMat const *A) {
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  if (A->props & BF_MAT_PROP_CONJ_TRANS)
    return CblasConjTrans;
  else if (A->props & BF_MAT_PROP_TRANS)
    return CblasTrans;
  else
    return CblasNoTrans;
}

static BfSize getLeadingDimension(BfMat const *A) {
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  BfSize dtype_size = bfDtypeSize(A->dtype);
  return bfMatIsTransposed(A) ?
    bfMatColStride(A)/dtype_size :
    bfMatRowStride(A)/dtype_size;
}

void bfMatAddInplace(BfMat *A, BfMat const *B) {
  if (A->dtype != B->dtype) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  if (A->dtype != BF_DTYPE_REAL && A->dtype != BF_DTYPE_COMPLEX) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  BfSize m = A->numRows;
  if (m != B->numRows) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  BfSize n = A->numCols;
  if (n != B->numCols) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  if (A->dtype == BF_DTYPE_REAL) {
    BfByte *A_row_ptr = A->data, *B_row_ptr = B->data;
    for (BfSize i = 0; i < m; ++i) {
      BfByte *A_col_ptr = A_row_ptr, *B_col_ptr = B_row_ptr;
      for (BfSize j = 0; j < n; ++j) {
        *(BfReal *)A_col_ptr += *(BfReal *)B_col_ptr;
        A_col_ptr += A->colStride;
        B_col_ptr += B->colStride;
      }
      A_row_ptr += A->rowStride;
      B_row_ptr += B->rowStride;
    }
  }

  if (A->dtype == BF_DTYPE_COMPLEX) {
    BfByte *A_row_ptr = A->data, *B_row_ptr = B->data;
    for (BfSize i = 0; i < m; ++i) {
      BfByte *A_col_ptr = A_row_ptr, *B_col_ptr = B_row_ptr;
      for (BfSize j = 0; j < n; ++j) {
        *(BfComplex *)A_col_ptr += *(BfComplex *)B_col_ptr;
        A_col_ptr += A->colStride;
        B_col_ptr += B->colStride;
      }
      A_row_ptr += A->rowStride;
      B_row_ptr += B->rowStride;
    }
  }
}

static void
complexComplexMatMul(BfMat const *A, BfMat const *B, BfMat *C)
{
  BEGIN_ERROR_HANDLING();

  assert(A->dtype == BF_DTYPE_COMPLEX);
  assert(B->dtype == BF_DTYPE_COMPLEX);

  assert(!(A->props & BF_MAT_PROP_DIAGONAL));
  assert(!(B->props & BF_MAT_PROP_DIAGONAL));

  assert(!bfMatIsInitialized(C));

  BfComplex alpha = 1, beta = 0;

  enum CBLAS_TRANSPOSE A_trans = getCblasTranspose(A);
  enum CBLAS_TRANSPOSE B_trans = getCblasTranspose(B);

  BfSize m = bfMatNumRows(A);
  BfSize n = bfMatNumCols(B);
  BfSize k = bfMatNumCols(A);
  assert(k == bfMatNumRows(B));
  bfInitEmptyMat(C, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE, m, n);
  HANDLE_ERROR();

  BfSize lda = getLeadingDimension(A);
  BfSize ldb = getLeadingDimension(B);
  BfSize ldc = getLeadingDimension(C);

  cblas_zgemm(CblasRowMajor, A_trans, B_trans, m, n, k,
              &alpha, A->data, lda, B->data, ldb, &beta, C->data, ldc);
  /* TODO: error-handling */

  END_ERROR_HANDLING() {
    bfFreeMat(C);
  }
}

void
bfMatMul(BfMat const *A, BfMat const *B, BfMat *C)
{
  assert(A->dtype != BF_DTYPE_MAT);
  assert(B->dtype != BF_DTYPE_MAT);

  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  assert(bfMatIsAligned(A));
  assert(bfMatIsAligned(B));

  if (bfMatNumCols(A) != bfMatNumRows(B)) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  if (A->dtype == BF_DTYPE_COMPLEX && B->dtype == BF_DTYPE_COMPLEX) {
    complexComplexMatMul(A, B, C);
    enum BfError error = bfGetError();
    if (error)
      bfSetError(error);
    return;
  }

  bfSetError(BF_ERROR_NOT_IMPLEMENTED);
}

void
realComplexMatSolve(BfMat const *A, BfMat const *B, BfMat *C)
{
  BEGIN_ERROR_HANDLING();

  assert(A->dtype == BF_DTYPE_REAL);
  assert(B->dtype == BF_DTYPE_COMPLEX);

  BfSize m = bfMatNumCols(A);
  BfSize k = bfMatNumRows(B);
  BfSize n = bfMatNumCols(B);

  bfInitEmptyMat(C, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE, m, n);
  HANDLE_ERROR();

  if (A->props & BF_MAT_PROP_DIAGONAL) {
    BfReal scale;
    BfVec row;
    for (BfSize i = 0; i < k; ++i) {
      bfCopyMatRow(B, i, C, i);
      row = bfGetMatRow(C, i);
      bfGetMatElt(A, i, i, &scale);
      HANDLE_ERROR();
      bfVecScaleByReal(&row, 1/scale);
    }
  }

  END_ERROR_HANDLING() {
    bfFreeMat(C);
  }
}

void
complexComplexMatSolve(BfMat const *A, BfMat const *B, BfMat *C)
{
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  if (A->props & BF_MAT_PROP_UNITARY ||
      (A->props & BF_MAT_PROP_SEMI_UNITARY
       && bfMatNumRows(A) > bfMatNumCols(A))) {
    BfMat AH = *A;
    AH.props ^= (AH.props & BF_MAT_PROP_CONJ_TRANS);
    complexComplexMatMul(&AH, B, C);
    enum BfError error = bfGetError();
    if (error) {
      bfSetError(error);
      return;
    }
  }
}

/* Solve the linear system B = A*C for C. Analogous to A\B = C in
 * MATLAB. */
void
bfMatSolve(BfMat const *A, BfMat const *B, BfMat *C)
{
  assert(A->dtype != BF_DTYPE_MAT);
  assert(B->dtype != BF_DTYPE_MAT);

  if (bfMatNumRows(A) != bfMatNumRows(B)) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  if (A->dtype == BF_DTYPE_REAL && B->dtype == BF_DTYPE_COMPLEX) {
    realComplexMatSolve(A, B, C);
    enum BfError error = bfGetError();
    if (error)
      bfSetError(error);
    return;
  }

  if (A->dtype == BF_DTYPE_COMPLEX && B->dtype == BF_DTYPE_COMPLEX) {
    complexComplexMatSolve(A, B, C);
    enum BfError error = bfGetError();
    if (error)
      bfSetError(error);
    return;
  }

  bfSetError(BF_ERROR_NOT_IMPLEMENTED);
}

void
initEmptySvdMatsComplex(BfMat const *A, BfMat *U, BfMat *S, BfMat *Vt)
{
  BEGIN_ERROR_HANDLING();

  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  BfSize m = A->numRows;
  BfSize n = A->numCols;
  BfSize p = m < n ? m : n;

  /* init U */
  bfInitEmptyMat(U, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE, m, p);
  HANDLE_ERROR();

  /* init S */
  bfInitEmptyMat(S, BF_DTYPE_REAL, BF_MAT_PROP_DIAGONAL, p, p);
  HANDLE_ERROR();

  /* init Vt */
  bfInitEmptyMat(Vt, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE, p, n);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    bfFreeMat(U);
    bfFreeMat(S);
    bfFreeMat(Vt);
  }
}

void
bfInitEmptySvdMats(BfMat const *A, BfMat *U, BfMat *S, BfMat *Vt)
{
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  if (A->dtype == BF_DTYPE_COMPLEX) {
    initEmptySvdMatsComplex(A, U, S, Vt);
    enum BfError error = bfGetError();
    if (error)
      bfSetError(error);
    return;
  }

  bfSetError(BF_ERROR_NOT_IMPLEMENTED);
}

void
computeMatSvdComplex(BfMat const *A, BfMat *U, BfMat *S, BfMat *Vt)
{
  BEGIN_ERROR_HANDLING();

  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  lapack_int m = A->numRows;
  lapack_int n = A->numCols;

  BfReal *superb = NULL;

  /* zgesvd will overwrite A, so allocate space for a copy */
  BfComplex *A_data_copy = malloc(m*n*sizeof(BfComplex));
  if (A_data_copy == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* copy contents of A */
  memcpy(A_data_copy, A->data, m*n*sizeof(BfComplex));
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
    LAPACK_ROW_MAJOR, 'S', 'S', m, n, A_data_copy, m, S->data, U->data, m,
    Vt->data, n, superb);

  /* check for invalid arguments */
  if (info < 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* check for errors */
  if (info > 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  U->props |= BF_MAT_PROP_UNITARY;
  S->props |= BF_MAT_PROP_DIAGONAL;
  Vt->props |= BF_MAT_PROP_UNITARY;

  END_ERROR_HANDLING() {}

  free(A_data_copy);
  free(superb);
}

void
bfComputeMatSvd(BfMat const *A, BfMat *U, BfMat *S, BfMat *Vt)
{
  BEGIN_ERROR_HANDLING();

  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  bfInitEmptySvdMats(A, U, S, Vt);
  HANDLE_ERROR();

  if (A->dtype == BF_DTYPE_COMPLEX) {
    computeMatSvdComplex(A, U, S, Vt);
    HANDLE_ERROR();
  } else {
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING() {
    bfFreeMat(U);
    bfFreeMat(S);
    bfFreeMat(Vt);
  }
}

void
bfMatLstSq(BfMat const *A, BfMat const *B, BfMat *C)
{
  BEGIN_ERROR_HANDLING();

  assert(!(A->props & BF_MAT_PROP_DIAGONAL));
  assert(A->dtype != BF_DTYPE_MAT);

  BfReal const atol = BF_EPS_MACH;

  BfSize m = A->numRows, n = A->numCols;

  BfReal const rtol = (m > n ? m : n)*BF_EPS_MACH;

  /* compute SVD of A */

  BfMat U = bfGetUninitializedMat();
  BfMat S = bfGetUninitializedMat();
  BfMat VH = bfGetUninitializedMat();
  bfComputeMatSvd(A, &U, &S, &VH);
  HANDLE_ERROR();

  /* compute tolerance and compute number of terms to
   * retain in pseudoinverse */

  BfSize k_max = m < n ? m : n;

  BfReal *sigma = S.data;

  BfReal tol = rtol*sigma[0] + atol;

  BfSize k;
  for (k = 0; k < k_max; ++k)
    if (sigma[k] < tol)
      break;

  /* get subblocks of truncated SVD */

  BfMat UkH = bfGetMatColRange(&U, 0, k);
  UkH = bfConjTrans(&UkH);

  BfMat Sk = bfGetMatContSubblock(&S, 0, k, 0, k);

  BfMat Vk = bfGetMatRowRange(&VH, 0, k);
  Vk = bfConjTrans(&Vk);

  /* solve least squares problem */

  BfMat tmp1 = bfGetUninitializedMat();
  bfMatMul(&UkH, B, &tmp1);
  HANDLE_ERROR();

  BfMat tmp2 = bfGetUninitializedMat();
  bfMatSolve(&Sk, &tmp1, &tmp2);
  HANDLE_ERROR();

  bfMatMul(&Vk, &tmp2, C);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    bfFreeMat(C);
  }

  bfFreeMat(&U);
  bfFreeMat(&S);
  bfFreeMat(&VH);
  bfFreeMat(&tmp1);
  bfFreeMat(&tmp2);
}
