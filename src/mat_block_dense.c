#include <bf/mat_block_dense.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>

#define INTERFACE BF_INTERFACE_Mat
BF_DEFINE_VTABLE(Mat, MatBlockDense)
#undef INTERFACE

#define INTERFACE BF_INTERFACE_MatBlock
BF_DEFINE_VTABLE(MatBlock, MatBlockDense)
#undef INTERFACE

/** Interface: Mat */

BF_STUB(BfMat *, MatBlockDenseCopy, BfMat const *)
BF_STUB(BfMat *, MatBlockDenseGetView, BfMat *)
BF_STUB(BfVec *, MatBlockDenseGetRowView, BfMat *, BfSize)
BF_STUB(BfVec *, MatBlockDenseGetColView, BfMat *, BfSize)
BF_STUB(BfVec *, MatBlockDenseGetColRangeView, BfMat *, BfSize, BfSize, BfSize)

void bfMatBlockDenseDelete(BfMat **mat) {
  bfMatBlockDenseDeinitAndDealloc((BfMatBlockDense **)mat);
}

BF_STUB(BfMat *, MatBlockDenseEmptyLike, BfMat const *, BfSize, BfSize)
BF_STUB(BfMat *, MatBlockDenseZerosLike, BfMat const *, BfSize, BfSize)

BfType bfMatBlockDenseGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_BLOCK_DENSE;
}

bool bfMatBlockDenseInstanceOf(BfMat const *mat, BfType type) {
  return bfTypeDerivedFrom(bfMatGetType(mat), type);
}

BfSize bfMatBlockDenseNumBytes(BfMat const *mat) {
  BfSize num_bytes = 0;

  BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);

  /* memory occupied by mat itself */
  num_bytes += sizeof(BfMatBlockDense);

  /* memory occupied by blocks */
  BfSize numBlocks = bfMatBlockDenseNumBlocks(matBlock);
  for (BfSize i = 0; i < numBlocks; ++i)
    num_bytes += bfMatNumBytes(matBlock->block[i]);

  /* mat->super.rowOffset and mat->super.colOffset */
  num_bytes += (mat->numRows + 1)*sizeof(BfSize);
  num_bytes += (mat->numCols + 1)*sizeof(BfSize);

  /* mat->rowInd and mat->colInd */
  num_bytes += 2*numBlocks*sizeof(BfSize);

  return num_bytes;
}

BF_STUB(void, MatBlockDenseSave, BfMat const *, char const *)
BF_STUB(void, MatBlockDensePrint, BfMat const *, FILE *)

BfSize bfMatBlockDenseGetNumRows(BfMat const *mat) {
  BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
  return bfMatIsTransposed(mat) ?
    matBlock->colOffset[mat->numCols] :
    matBlock->rowOffset[mat->numRows];
}

BfSize bfMatBlockDenseGetNumCols(BfMat const *mat) {
  BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
  return bfMatIsTransposed(mat) ?
    matBlock->rowOffset[mat->numRows] :
    matBlock->colOffset[mat->numCols];
}

BF_STUB(void, MatBlockDenseSetRow, BfMat *, BfSize, BfVec const *)
BF_STUB(void, MatBlockDenseSetCol, BfMat *, BfSize, BfVec const *)
BF_STUB(void, MatBlockDenseSetColRange, BfMat *, BfSize, BfSize, BfSize, BfVec const *)

BfMat *bfMatBlockDenseGetRowRange(BfMat *mat, BfSize i0, BfSize i1) {
  BEGIN_ERROR_HANDLING();

  BfMatBlock *matBlock = bfMatToMatBlock(mat);

  BfMatBlockDense *view = NULL;

  if (i0 >= i1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numRows = bfMatBlockDenseGetNumRows(mat);

  if (i0 >= numRows)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  if (i1 >= numRows)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  /* Find the offset of each row of blocks. For now, we require that
   * i0 and i1 exactly correspond to the start of two distinct block
   * rows.
   *
   * In the future, we could relax this assumption and allow block
   * matrices to get partitioned when i0 and i1 don't exactly align
   * with the block partitions. But probably it would be better to
   * just do that in a separate function. */

  BfSize p0 = bfMatBlockFindRowBlock(matBlock, i0);
  if (p0 == BF_SIZE_BAD_VALUE || i0 != matBlock->rowOffset[p0])
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize p1 = bfMatBlockFindRowBlock(matBlock, i1);
  if (p1 == BF_SIZE_BAD_VALUE || i1 != matBlock->rowOffset[p1])
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numColBlocks = bfMatBlockGetNumColBlocks(matBlock);

  view = bfMatBlockDenseNew();
  HANDLE_ERROR();

  bfMatBlockDenseInit(view, p1 - p0, numColBlocks);

  /* the row range is a view of the original matrix */
  bfMatBlockDenseToMat(view)->props |= BF_MAT_PROPS_VIEW;

  /* point to the subrange of blocks */
  view->super.block = matBlock->block + numColBlocks*p0;

  END_ERROR_HANDLING() {}

  return bfMatBlockDenseToMat(view);
}

BF_STUB(BfMat *, MatBlockDenseGetColRange, BfMat *, BfSize, BfSize)
BF_STUB(BfMat *, MatBlockDenseGetColRangeCopy, BfMat const *, BfSize, BfSize)
BF_STUB(BfMat *, MatBlockDenseGetRowRangeCopy, BfMat const *, BfSize, BfSize)
BF_STUB(void, MatBlockDenseSetRowRange, BfMat *, BfSize, BfSize, BfMat const *)
BF_STUB(void, MatBlockDensePermuteRows, BfMat *, BfPerm const *)
BF_STUB(void, MatBlockDensePermuteCols, BfMat *, BfPerm const *)
BF_STUB(BfVec *, MatBlockDenseRowDists, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatBlockDenseColDists, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatBlockDenseColDots, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatBlockDenseColNorms, BfMat const *)

void bfMatBlockDenseScaleCols(BfMat *mat, BfVec const *vec) {
  BEGIN_ERROR_HANDLING();

  BfMatBlock *matBlock = NULL;
  BfMatBlockDense *matBlockDense = NULL;

  matBlock = bfMatToMatBlock(mat);
  HANDLE_ERROR();

  matBlockDense = bfMatToMatBlockDense(mat);
  HANDLE_ERROR();

  if (bfMatGetNumCols(mat) != vec->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize m_blk = bfMatBlockDenseGetNumRowBlocks(matBlock);
  BfSize n_blk = bfMatBlockDenseGetNumColBlocks(matBlock);

  for (BfSize i_blk = 0; i_blk < m_blk; ++i_blk) {
    for (BfSize j_blk = 0; j_blk < n_blk; ++j_blk) {
      BfSize j0 = matBlockDense->super.colOffset[j_blk];
      BfSize j1 = matBlockDense->super.colOffset[j_blk + 1];
      BfVec *subvec = bfVecGetSubvecCopy(vec, j0, j1);
      BfMat *block = bfMatBlockDenseGetBlock(matBlockDense, i_blk, j_blk);
      bfMatScaleCols(block, subvec);
      bfVecDelete(&subvec);
    }
  }

  END_ERROR_HANDLING() {}
}

BF_STUB(BfVec *, MatBlockDenseSumCols, BfMat const *)

void bfMatBlockDenseAddInplace(BfMat *mat, BfMat const *otherMat) {
  BEGIN_ERROR_HANDLING();

  BfMatBlock *matBlock = NULL;
  BfMatBlockDense *matBlockDense = NULL;

  matBlock = bfMatToMatBlock(mat);
  HANDLE_ERROR();

  matBlockDense = bfMatToMatBlockDense(mat);
  HANDLE_ERROR();

  if (bfMatGetNumRows(mat) != bfMatGetNumRows(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (bfMatGetNumCols(mat) != bfMatGetNumCols(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize m_blk = bfMatBlockDenseGetNumRowBlocks(matBlock);
  BfSize n_blk = bfMatBlockDenseGetNumColBlocks(matBlock);

  for (BfSize i_blk = 0; i_blk < m_blk; ++i_blk) {
    BfSize i0 = matBlockDense->super.rowOffset[i_blk];
    BfSize i1 = matBlockDense->super.rowOffset[i_blk + 1];
    BfMat const *otherRows = bfMatGetRowRangeCopy(otherMat, i0, i1);

    for (BfSize j_blk = 0; j_blk < n_blk; ++j_blk) {
      BfSize j0 = matBlockDense->super.colOffset[j_blk];
      BfSize j1 = matBlockDense->super.colOffset[j_blk + 1];
      BfMat const *otherBlock = bfMatGetColRangeCopy(otherRows, j0, j1);

      BfMat *block = bfMatBlockDenseGetBlock(matBlockDense, i_blk, j_blk);

      bfMatAddInplace(block, otherBlock);
    }
  }

  END_ERROR_HANDLING() {}
}

BF_STUB(void, MatBlockDenseAddDiag, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatBlockDenseSub, BfMat const *, BfMat const *)
BF_STUB(void, MatBlockDenseSubInplace, BfMat *, BfMat const *)

BfMat *bfMatBlockDenseMul(BfMat const *mat, BfMat const *otherMat) {
  BEGIN_ERROR_HANDLING();

  BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
  BfMatBlockDense const *matBlockDense = bfMatConstToMatBlockDenseConst(mat);

  BfSize numRows, numCols, numRowBlocks, numColBlocks;
  BfMat *block = NULL, *op2Rows = NULL;
  BfMat *result = NULL, *resultRows = NULL, *tmp = NULL;

  numRowBlocks = bfMatBlockGetNumRowBlocks(matBlock);
  numColBlocks = bfMatBlockGetNumColBlocks(matBlock);

  numRows = bfMatGetNumRows(mat);
  numCols = bfMatGetNumCols(otherMat);

  result = bfMatZerosLike(otherMat, numRows, numCols);
  if (result == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0, i0, i1; i < numRowBlocks; ++i) {
    i0 = matBlock->rowOffset[i];
    i1 = matBlock->rowOffset[i + 1];
    resultRows = bfMatGetRowRange(result, i0, i1);
    for (BfSize j = 0, j0, j1; j < numColBlocks; ++j) {
      j0 = matBlock->colOffset[j];
      j1 = matBlock->colOffset[j + 1];
      op2Rows = bfMatGetRowRange((BfMat *)otherMat, j0, j1);
      block = bfMatBlockDenseGetBlock((BfMatBlockDense *)matBlockDense, i, j);
      assert(bfMatGetNumRows(block) == i1 - i0);
      assert(bfMatGetNumCols(block) == j1 - j0);
      tmp = bfMatMul(block, op2Rows);
      bfMatAddInplace(resultRows, tmp);
      bfMatDelete(&tmp);
    }
  }

  END_ERROR_HANDLING()
    bfMatDelete(&result);

  return result;
}

BF_STUB(BfVec *, MatBlockDenseMulVec, BfMat const *, BfVec const *)
BF_STUB(void, MatBlockDenseMulInplace, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatBlockDenseSolveLU, BfMat const *, BfMat const *)
BF_STUB(BfMat *, MatBlockDenseLstSq, BfMat const *, BfMat const *)
BF_STUB(BfMat *, MatBlockDenseGetGivensRotation, BfVec const *, BfSize, BfSize)
BF_STUB(bool, MatBlockDenseIsUpperTri, BfMat const *)
BF_STUB(BfVec *, MatBlockDenseBackwardSolveVec, BfMat const *, BfVec const *)
BF_STUB(bool, MatBlockDenseIsZero, BfMat const *)

/** Interface: MatBlock */

BfSize bfMatBlockDenseNumBlocks(BfMatBlock const *matBlock) {
  BfMat const *mat = bfMatBlockConstToMatConst(matBlock);
  return mat->numRows*mat->numCols;
}

BfSize bfMatBlockDenseGetNumRowBlocks(BfMatBlock const *matBlock) {
  BfMat const *mat = bfMatBlockConstToMatConst(matBlock);
  if (bfMatGetType(mat) != BF_TYPE_MAT_BLOCK_DENSE) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return BF_SIZE_BAD_VALUE;
  } else {
    return mat->numRows;
  }
}

BfSize bfMatBlockDenseGetNumColBlocks(BfMatBlock const *matBlock) {
  BfMat const *mat = bfMatBlockConstToMatConst(matBlock);
  if (bfMatGetType(mat) != BF_TYPE_MAT_BLOCK_DENSE) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return BF_SIZE_BAD_VALUE;
  } else {
    return mat->numCols;
  }
}

BF_STUB(BfSize, MatBlockDenseGetNumBlockRows, BfMatBlock const *, BfSize)
BF_STUB(BfSize, MatBlockDenseGetNumBlockCols, BfMatBlock const *, BfSize)

/** Upcasting: */

BfMat *bfMatBlockDenseToMat(BfMatBlockDense *matBlockDense) {
  return &matBlockDense->super.super;
}

BfMat const *bfMatBlockDenseConstToMatConst(BfMatBlockDense const *matBlockDense) {
  return &matBlockDense->super.super;
}

BfMatBlock *bfMatBlockDenseToMatBlock(BfMatBlockDense *matBlockDense) {
  return &matBlockDense->super;
}

/** Downcasting: */

BfMatBlockDense *bfMatToMatBlockDense(BfMat *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_BLOCK_DENSE)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatBlockDense *)mat;
  }
}

BfMatBlockDense const *bfMatConstToMatBlockDenseConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_BLOCK_DENSE)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatBlockDense const *)mat;
  }
}

/** Implementation: */

BfMatBlockDense *bfMatBlockDenseNew() {
  return malloc(sizeof(BfMatBlockDense));
}

void bfMatBlockDenseInit(BfMatBlockDense *mat,
                         BfSize numBlockRows, BfSize numBlockCols) {
  BEGIN_ERROR_HANDLING();

  BfSize numBlocks = numBlockRows*numBlockCols;

  bfMatBlockInit(&mat->super,
                 &MatVtbl, &MatBlockVtbl,
                 numBlocks, numBlockRows, numBlockCols);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfMatBlockDeinit(&mat->super);
}

void bfMatBlockDenseDeinit(BfMatBlockDense *mat) {
  bfMatBlockDeinit(&mat->super);
}

void bfMatBlockDenseDealloc(BfMatBlockDense **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatBlockDenseDeinitAndDealloc(BfMatBlockDense **mat) {
  bfMatBlockDenseDeinit(*mat);
  bfMatBlockDenseDealloc(mat);
}

BfMat *bfMatBlockDenseGetBlock(BfMatBlockDense *mat, BfSize i, BfSize j) {
  BEGIN_ERROR_HANDLING();

  BfMat *block = NULL;

  BfSize numBlockRows = mat->super.super.numRows;
  if (i >= numBlockRows)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  BfSize numBlockCols = mat->super.super.numCols;
  if (j >= numBlockCols)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  block = mat->super.block[numBlockCols*i + j];

  END_ERROR_HANDLING() {}

  return block;
}

void bfMatBlockDenseSetBlock(BfMatBlockDense *mat, BfSize i, BfSize j,
                             BfMat *block) {
  BEGIN_ERROR_HANDLING();

  BfSize numBlockRows = mat->super.super.numRows;
  if (i >= numBlockRows)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  BfSize numBlockCols = mat->super.super.numCols;
  if (j >= numBlockCols)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  mat->super.block[numBlockCols*i + j] = block;

  END_ERROR_HANDLING() {}
}
