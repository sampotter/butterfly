#include <bf/mat_block_dense.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_coo_complex.h>
#include <bf/mat_sum.h>
#include <bf/mat_zero.h>
#include <bf/mem.h>
#include <bf/size_array.h>
#include <bf/util.h>
#include <bf/vec_complex.h>
#include <bf/vec_real.h>

/** Helper macros: */

#define NUM_ROW_BLOCKS(mat) mat->super.super.numRows
#define NUM_COL_BLOCKS(mat) mat->super.super.numCols
#define NUM_BLOCKS(mat) NUM_ROW_BLOCKS(mat)*NUM_COL_BLOCKS(mat)
#define ROW_OFFSET(mat, i) mat->super.rowOffset[i]
#define COL_OFFSET(mat, j) mat->super.colOffset[j]
#define BLOCK_ROW_OFFSET(mat, k) ROW_OFFSET(mat, k/NUM_COL_BLOCKS(mat))
#define BLOCK_COL_OFFSET(mat, k) COL_OFFSET(mat, k % NUM_COL_BLOCKS(mat))
#define NUM_BLOCK_ROWS(mat, i) ROW_OFFSET(mat, i + 1) - ROW_OFFSET(mat, i)
#define NUM_BLOCK_COLS(mat, j) COL_OFFSET(mat, j + 1) - COL_OFFSET(mat, j)
#define BLOCK(mat, i, j) mat->super.block[i*NUM_COL_BLOCKS(mat) + j]
#define LINEAR_BLOCK(mat, k) mat->super.block[k]

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .Copy = (__typeof__(&bfMatBlockDenseCopy))bfMatBlockDenseCopy,
  .Steal = (__typeof__(&bfMatSteal))bfMatBlockDenseSteal,
  .GetRowCopy = (__typeof__(&bfMatBlockDenseGetRowCopy))bfMatBlockDenseGetRowCopy,
  .Delete = (__typeof__(&bfMatBlockDenseDelete))bfMatBlockDenseDelete,
  .GetType = (__typeof__(&bfMatBlockDenseGetType))bfMatBlockDenseGetType,
  .NumBytes = (__typeof__(&bfMatBlockDenseNumBytes))bfMatBlockDenseNumBytes,
  .Save = (__typeof__(&bfMatSave))bfMatBlockDenseSave,
  .Dump = (__typeof__(&bfMatDump))bfMatBlockDenseDump,
  .GetNumRows = (__typeof__(&bfMatBlockDenseGetNumRows))bfMatBlockDenseGetNumRows,
  .GetNumCols = (__typeof__(&bfMatBlockDenseGetNumCols))bfMatBlockDenseGetNumCols,
  .GetRowRange = (__typeof__(&bfMatBlockDenseGetRowRange))bfMatBlockDenseGetRowRange,
  .GetRowRangeCopy = (__typeof__(&bfMatGetRowRangeCopy))bfMatBlockDenseGetRowRangeCopy,
  .ScaleCols = (__typeof__(&bfMatScaleCols))bfMatBlockDenseScaleCols,
  .AddInplace = (__typeof__(&bfMatAddInplace))bfMatBlockDenseAddInplace,
  .Mul = (__typeof__(&bfMatBlockDenseMul))bfMatBlockDenseMul,
  .MulVec = (__typeof__(&bfMatMulVec))bfMatBlockDenseMulVec,
  .RmulVec = (__typeof__(&bfMatRmulVec))bfMatBlockDenseRmulVec,
  .ToType = (__typeof__(&bfMatBlockDenseToType))bfMatBlockDenseToType,
  .GetNonzeroColumnRanges = (__typeof__(&bfMatGetNonzeroColumnRanges))bfMatBlockDenseGetNonzeroColumnRanges,
  .PrintBlocksDeep = (__typeof__(&bfMatPrintBlocksDeep))bfMatBlockDensePrintBlocksDeep
};

BfMat *bfMatBlockDenseCopy(BfMat const *mat) {
  BF_ERROR_BEGIN();

  BfMatBlock const *matBlock = NULL;
  BfMatBlockDense const *matBlockDense = NULL;
  BfMatBlockDense *matBlockDenseCopy = NULL;
  BfMat const *block = NULL;
  BfMat *blockCopy = NULL;

  matBlock = bfMatConstToMatBlockConst(mat);
  HANDLE_ERROR();

  matBlockDense = bfMatConstToMatBlockDenseConst(mat);
  HANDLE_ERROR();

  BfSize numRowBlocks = bfMatBlockGetNumRowBlocks(matBlock);
  BfSize numColBlocks = bfMatBlockGetNumColBlocks(matBlock);

  matBlockDenseCopy = bfMatBlockDenseNew();
  HANDLE_ERROR();

  bfMatBlockDenseInit(matBlockDenseCopy, numRowBlocks, numColBlocks);
  HANDLE_ERROR();

  for (BfSize i = 0; i <= numRowBlocks; ++i)
    matBlockDenseCopy->super.rowOffset[i] = matBlock->rowOffset[i];

  for (BfSize j = 0; j <= numColBlocks; ++j)
    matBlockDenseCopy->super.colOffset[j] = matBlock->colOffset[j];

  for (BfSize i = 0; i < numRowBlocks; ++i) {
    for (BfSize j = 0; j < numColBlocks; ++j) {
      block = bfMatBlockDenseGetBlockConst(matBlockDense, i, j);
      blockCopy = bfMatCopy(block);
      HANDLE_ERROR();
      bfMatBlockDenseSetBlock(matBlockDenseCopy, i, j, blockCopy);
    }
  }

  BF_ERROR_END() {
    bfMatDelete(&blockCopy);
    bfMatBlockDenseDeinitAndDealloc(&matBlockDenseCopy);
  }

  return bfMatBlockDenseToMat(matBlockDenseCopy);
}

BfMat *bfMatBlockDenseSteal(BfMatBlockDense *matBlockDense) {
  BF_ERROR_BEGIN();

  BfMat *mat = bfMatBlockDenseToMat(matBlockDense);

  if (bfMatIsView(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatBlockDense *matBlockDenseNew = bfMatBlockDenseNew();
  HANDLE_ERROR();

  *matBlockDenseNew = *matBlockDense;

  mat->props |= BF_MAT_PROPS_VIEW;

  BF_ERROR_END() {
    BF_DIE();
  }

  return bfMatBlockDenseToMat(matBlockDenseNew);
}

BfVec *bfMatBlockDenseGetRowCopy(BfMat const *mat, BfSize i) {
  BF_ERROR_BEGIN();

  BfVec *rowCopy = NULL;
  BfMatBlock const *matBlock = NULL;
  BfMatBlockDense const *matBlockDense = NULL;

  matBlock = bfMatConstToMatBlockConst(mat);
  HANDLE_ERROR();

  matBlockDense = bfMatConstToMatBlockDenseConst(mat);
  HANDLE_ERROR();

  BfSize numRows = bfMatGetNumRows(mat);
  if (i >= numRows)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numRowBlocks = bfMatBlockGetNumRowBlocks(matBlock);
  if (numRowBlocks == 0)
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfSize ib = 0, i0, i1;
  for (; ib < numRowBlocks; ++ib) {
    i0 = ROW_OFFSET(matBlockDense, ib);
    i1 = ROW_OFFSET(matBlockDense, ib + 1);
    if (i0 <= i && i < i1)
      break;
  }

  BfSize numColBlocks = bfMatBlockGetNumColBlocks(matBlock);

  for (BfSize jb = 0; jb < numColBlocks; ++jb) {
    BfMat *block = BLOCK(matBlockDense, ib, jb);
    BfVec *blockRowCopy = bfMatGetRowCopy(block, i - i0);
    BF_ASSERT(blockRowCopy->size == NUM_BLOCK_COLS(matBlockDense, jb));
    if (rowCopy == NULL) {
      rowCopy = blockRowCopy;
    } else {
      BfVec *cat = bfVecConcat(rowCopy, blockRowCopy);
      HANDLE_ERROR();
      bfVecDelete(&rowCopy);
      bfVecDelete(&blockRowCopy);
      rowCopy = cat;
    }
    BF_ASSERT(rowCopy->size == matBlockDense->super.colOffset[jb + 1]);
  }

  BF_ASSERT(rowCopy->size == bfMatGetNumCols(mat));

  BF_ERROR_END() {}

  return rowCopy;
}

void bfMatBlockDenseDelete(BfMat **mat) {
  bfMatBlockDenseDeinitAndDealloc((BfMatBlockDense **)mat);
}

BfType bfMatBlockDenseGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_BLOCK_DENSE;
}

BfSize bfMatBlockDenseNumBytes(BfMat const *mat) {
  BfSize num_bytes = 0;

  BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
  BfMatBlockDense const *matBlockDense = bfMatConstToMatBlockDenseConst(mat);

  /* memory occupied by mat itself */
  num_bytes += sizeof(BfMatBlockDense);

  /* memory occupied by blocks */
  BfSize numBlocks = bfMatBlockDenseNumBlocks(matBlockDense);
  for (BfSize i = 0; i < numBlocks; ++i)
    num_bytes += bfMatNumBytes(matBlock->block[i]);

  /* mat->super.rowOffset and mat->super.colOffset */
  num_bytes += (mat->numRows + 1)*sizeof(BfSize);
  num_bytes += (mat->numCols + 1)*sizeof(BfSize);

  /* mat->rowInd and mat->colInd */
  num_bytes += 2*numBlocks*sizeof(BfSize);

  return num_bytes;
}

void bfMatBlockDenseSave(BfMatBlockDense const *matBlockDense, char const *path) {
  BF_ERROR_BEGIN();

  FILE *fp = fopen(path, "w");

  bfMatBlockDenseDump(matBlockDense, fp);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  fclose(fp);
}

void bfMatBlockDenseDump(BfMatBlockDense const *matBlockDense, FILE *fp) {
  BF_ERROR_BEGIN();

  /* Serialize the number of row and column blocks: */
  fwrite(&NUM_ROW_BLOCKS(matBlockDense), sizeof(BfSize), 1, fp);
  fwrite(&NUM_COL_BLOCKS(matBlockDense), sizeof(BfSize), 1, fp);

  /* Serialize the row and column offsets: */
  fwrite(&ROW_OFFSET(matBlockDense, 0), sizeof(BfSize), NUM_ROW_BLOCKS(matBlockDense) + 1, fp);
  fwrite(&COL_OFFSET(matBlockDense, 0), sizeof(BfSize), NUM_COL_BLOCKS(matBlockDense) + 1, fp);

  /* Recursively write each of the blocks: */
  for (BfSize k = 0; k < NUM_BLOCKS(matBlockDense); ++k) {
    BfMat const *block = LINEAR_BLOCK(matBlockDense, k);
    bfMatDump(block, fp);
    HANDLE_ERROR();
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

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

BfMat *bfMatBlockDenseGetRowRange(BfMat *mat, BfSize i0, BfSize i1) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END() {}

  return bfMatBlockDenseToMat(view);
}

BfMat *bfMatBlockDenseGetRowRangeCopy(BfMatBlockDense const *matBlockDense,
                                      BfSize i0, BfSize i1) {
  BF_ERROR_BEGIN();

  if (i0 > i1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMat const *mat = bfMatBlockDenseConstToMatConst(matBlockDense);

  BfSize numRows = bfMatBlockDenseGetNumRows(mat);

  /* Check whether i0 <= i1 index a valid row range: */
  if (i0 == i1 && i0 > numRows)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);
  if (i1 > numRows)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  /** Seek through the row offsets to find i0 and i1. Since there
   * could be repeat row offsets, we need to search from above and
   * below here since we want to make sure to capture *all* row
   * offsets bracketed by i0 and i1.
   *
   * Afterwards, [p0, p1) indexed the block range corresponding to the
   * row range [i0, i1). */

  BfSize p0 = 0;
  while (p0 < NUM_ROW_BLOCKS(matBlockDense) &&
         ROW_OFFSET(matBlockDense, p0) < i0)
    ++p0;

  BfSize p1 = NUM_ROW_BLOCKS(matBlockDense);
  while (p1 > p0 && ROW_OFFSET(matBlockDense, p1 - 1) > i0)
    --p1;

  BF_ASSERT(p0 <= p1);

  if (p0 == p1) {
    --p0;
    BF_ASSERT(ROW_OFFSET(matBlockDense, p0) <= i0
              && i1 <= ROW_OFFSET(matBlockDense, p1));
  }

  /* The range [i0_, i1) is the row index range of every block in the
   * row with respect to `matBlockDense`'s index system. */
  BfSize i0_ = ROW_OFFSET(matBlockDense, p0);
  BfSize i1_ = ROW_OFFSET(matBlockDense, p1);

  if (i0_ == i1_) {
    BF_ASSERT(i0 <= i0_ && i0_ <= i1);
  } else {
    BF_ASSERT(i1_ >= i0 || i0_ <= i1);
  }

  BfSize m_ = i1_ - i0_;

  /* The range of rows which should be copied from each block in the
   * row with respect to the block row index system. */
  BfSize i0__ = i0_ < i0 ? i0 - i0_ : 0;
  BfSize i1__ = m_ - (i1 < i1_ ? i1_ - i1 : 0);
  BF_ASSERT(i0__ <= i1__ && i1__ <= m_);

  /** Copy all blocks in the row range and assemble the resulting
   * dense block matrix containing the row range. */

  /* Array holding all `(p1 - p0)*NUM_COL_BLOCKS(matBlockDense)`
   * blocks in row major order. */
  BfPtrArray *blocks = bfPtrArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  for (BfSize p = p0; p < p1; ++p) {
    for (BfSize q = 0; q < NUM_COL_BLOCKS(matBlockDense); ++q) {
      BfMat const *block = BLOCK(matBlockDense, p, q);

      BF_ASSERT(bfMatGetNumRows(block) == m_);

      BfMat *blockRowRangeCopy = bfMatGetRowRangeCopy(block, i0__, i1__);
      HANDLE_ERROR();

      bfPtrArrayAppend(blocks, blockRowRangeCopy);
      HANDLE_ERROR();
    }
  }

  BfMatBlockDense *result = bfMatBlockDenseNewFromBlocks(
    p1 - p0, NUM_COL_BLOCKS(matBlockDense), blocks, BF_POLICY_STEAL);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  bfPtrArrayDelete(&blocks);

  return bfMatBlockDenseToMat(result);
}

void bfMatBlockDenseScaleCols(BfMatBlockDense *matBlockDense, BfVec const *vec) {
  BF_ERROR_BEGIN();

  BfMat *mat = bfMatBlockDenseToMat(matBlockDense);

  if (bfMatGetNumCols(mat) != vec->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize m_blk = bfMatBlockDenseGetNumRowBlocks(matBlockDense);
  BfSize n_blk = bfMatBlockDenseGetNumColBlocks(matBlockDense);

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

  BF_ERROR_END() {}
}

void bfMatBlockDenseAddInplace(BfMatBlockDense *matBlockDense, BfMat const *otherMat) {
  BF_ERROR_BEGIN();

  BfMat *mat = bfMatBlockDenseToMat(matBlockDense);

  if (bfMatGetNumRows(mat) != bfMatGetNumRows(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (bfMatGetNumCols(mat) != bfMatGetNumCols(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize m_blk = bfMatBlockDenseGetNumRowBlocks(matBlockDense);
  BfSize n_blk = bfMatBlockDenseGetNumColBlocks(matBlockDense);

  for (BfSize i_blk = 0; i_blk < m_blk; ++i_blk) {
    BfSize i0 = matBlockDense->super.rowOffset[i_blk];
    BfSize i1 = matBlockDense->super.rowOffset[i_blk + 1];
    BfMat const *otherRows = bfMatGetRowRangeCopy(otherMat, i0, i1);

    for (BfSize j_blk = 0; j_blk < n_blk; ++j_blk) {
      BfSize j0 = matBlockDense->super.colOffset[j_blk];
      BfSize j1 = matBlockDense->super.colOffset[j_blk + 1];
      BfMat *otherBlock = bfMatGetColRangeCopy(otherRows, j0, j1);

      if (bfMatIsZero(otherBlock))
        continue;

      BfMat *block = bfMatBlockDenseGetBlock(matBlockDense, i_blk, j_blk);

      BfMat *newBlock = NULL;

      if (bfMatInstanceOf(block, BF_TYPE_MAT_PRODUCT)) {
        BfMatSum *matSum = bfMatSumNew();
        bfMatSumInit(matSum);
        bfMatSumAddTerm(matSum, block);
        bfMatSumAddTerm(matSum, otherBlock);
        newBlock = bfMatSumToMat(matSum);
      }

      if (newBlock == NULL)
        bfMatAddInplace(block, otherBlock);
      else
        bfMatBlockDenseSetBlock(matBlockDense, i_blk, j_blk, newBlock);
    }
  }

  BF_ERROR_END() {}
}

BfMat *bfMatBlockDenseMul(BfMat const *mat, BfMat const *otherMat) {
  BF_ERROR_BEGIN();

  BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
  BfMatBlockDense const *matBlockDense = bfMatConstToMatBlockDenseConst(mat);

  BfSize numRows, numCols, numRowBlocks, numColBlocks;
  BfMat const *block = NULL, *op2Rows = NULL;
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
      block = bfMatBlockDenseGetBlockConst(matBlockDense, i, j);
      BF_ASSERT(bfMatGetNumRows(block) == i1 - i0);
      BF_ASSERT(bfMatGetNumCols(block) == j1 - j0);
      tmp = bfMatMul(block, op2Rows);
      bfMatAddInplace(resultRows, tmp);
      bfMatDelete(&tmp);
    }
  }

  BF_ERROR_END()
    bfMatDelete(&result);

  return result;
}

BfVec *bfMatBlockDenseMulVec(BfMatBlockDense const *matBlockDense,
                             BfVec const *vec) {
  BF_ERROR_BEGIN();

  BfMat const *mat = bfMatBlockDenseConstToMatConst(matBlockDense);

  BfSize numRows = bfMatGetNumRows(mat);
  BfSize numRowBlocks = bfMatBlockDenseGetNumRowBlocks(matBlockDense);
  BfSize numColBlocks = bfMatBlockDenseGetNumColBlocks(matBlockDense);

  BfVec *result = NULL;
  switch (bfVecGetType(vec)) {
  case BF_TYPE_VEC_REAL:
    result = bfVecRealToVec(bfVecRealNewWithValue(numRows, 0));
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  for (BfSize i = 0; i < numRowBlocks; ++i) {
    BfSize i0 = ROW_OFFSET(matBlockDense, i);
    BfSize i1 = ROW_OFFSET(matBlockDense, i + 1);

    /* If this is a degenerate block row, skip it */
    if (i0 == i1)
      continue;

    BfVec *resultSubvecView = bfVecGetSubvecView(result, i0, i1);
    HANDLE_ERROR();

    for (BfSize j = 0; j < numColBlocks; ++j) {
      BfMat const *block = BLOCK(matBlockDense, i, j);

      BfSize j0 = COL_OFFSET(matBlockDense, j);
      BfSize j1 = COL_OFFSET(matBlockDense, j + 1);

      /* Skip degenerate block columns */
      if (j0 == j1)
        continue;

      /* Sanity check... basically unnecessary, but whatever */
      BF_ASSERT(bfMatGetNumRows(block) == i1 - i0);
      BF_ASSERT(bfMatGetNumCols(block) == j1 - j0);

      BfVec const *subvecView = bfVecGetSubvecViewConst(vec, j0, j1);
      HANDLE_ERROR();

      BfVec *tmp = bfMatMulVec(block, subvecView);
      HANDLE_ERROR();

      bfVecAddInplace(resultSubvecView, tmp);

      bfVecDelete(&tmp);
    }

    bfVecDelete(&resultSubvecView);
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return result;
}

BfVec *bfMatBlockDenseRmulVec(BfMatBlockDense const *matBlockDense, BfVec const *vec) {
  BF_ERROR_BEGIN();

  BfMat const *mat = bfMatBlockDenseConstToMatConst(matBlockDense);

  BfSize numCols = bfMatGetNumCols(mat);
  BfSize numRowBlocks = bfMatBlockDenseGetNumRowBlocks(matBlockDense);
  BfSize numColBlocks = bfMatBlockDenseGetNumColBlocks(matBlockDense);

  BfVec *result = NULL;
  switch (bfVecGetType(vec)) {
  case BF_TYPE_VEC_REAL:
    result = bfVecRealToVec(bfVecRealNewWithValue(numCols, 0));
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  for (BfSize j = 0; j < numColBlocks; ++j) {
    BfSize j0 = COL_OFFSET(matBlockDense, j);
    BfSize j1 = COL_OFFSET(matBlockDense, j + 1);

    /* If this is a degenerate block row, skip it */
    if (j0 == j1)
      continue;

    BfVec *resultSubvecView = bfVecGetSubvecView(result, j0, j1);
    HANDLE_ERROR();

    for (BfSize i = 0; i < numRowBlocks; ++i) {
      BfMat const *block = BLOCK(matBlockDense, i, j);

      BfSize i0 = ROW_OFFSET(matBlockDense, i);
      BfSize i1 = ROW_OFFSET(matBlockDense, i + 1);

      /* Skip degenerate block columns */
      if (i0 == i1)
        continue;

      /* Sanity check... basically unnecessary, but whatever */
      BF_ASSERT(bfMatGetNumRows(block) == i1 - i0);
      BF_ASSERT(bfMatGetNumCols(block) == j1 - j0);

      BfVec const *subvecView = bfVecGetSubvecViewConst(vec, i0, i1);
      HANDLE_ERROR();

      BfVec *tmp = bfMatRmulVec(block, subvecView);
      HANDLE_ERROR();

      bfVecAddInplace(resultSubvecView, tmp);

      bfVecDelete(&tmp);
    }

    bfVecDelete(&resultSubvecView);
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return result;
}

static BfMat *toType_cooComplex(BfMat const *mat) {
  BF_ERROR_BEGIN();

  BfMatCooComplex *matCooComplex = NULL;

  BfSize numRows = bfMatGetNumRows(mat);
  BfSize numCols = bfMatGetNumCols(mat);

  BfSize *rowInd = NULL;
  BfSize *colInd = NULL;
  BfComplex *value = NULL;

  BfSize numElts = 0, capacity = 16;

  rowInd = bfMemAlloc(capacity, sizeof(BfSize));
  if (rowInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  colInd = bfMemAlloc(capacity, sizeof(BfSize));
  if (colInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  value = bfMemAlloc(capacity, sizeof(BfComplex));
  if (value == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0; i < numRows; ++i) {
    BfVec *vec = bfMatGetRowCopy(mat, i);
    HANDLE_ERROR();

    BfVecComplex *vecComplex = bfVecToVecComplex(vec);
    HANDLE_ERROR();

    BfComplex *ptr = vecComplex->data;

    /* Decrement the pointer here to simplify the logic of the
     * following loop. */
    --ptr;

    /* Accumulate nonzero values and their row and column indices into
     * our dynamic arrays */
    for (BfSize j = 0; j < vec->size; ++j) {
      if (*++ptr == 0)
        continue;

      /* If we've met the capacity of our dynamic arrays, increase
       * their capacity and reallocate */
      if (numElts == capacity) {
        capacity *= 2;

        BfSize *newRowInd = bfMemRealloc(rowInd, capacity, sizeof(BfSize));
        if (newRowInd == NULL)
          RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
        rowInd = newRowInd;

        BfSize *newColInd = bfMemRealloc(colInd, capacity, sizeof(BfSize));
        if (newColInd == NULL)
          RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
        colInd = newColInd;

        BfComplex *newValue = bfMemRealloc(value, capacity, sizeof(BfComplex));
        if (newValue == NULL)
          RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
        value = newValue;
      }

      rowInd[numElts] = i;
      colInd[numElts] = j;
      value[numElts++] = *ptr;
    }

    bfVecDelete(&vec);
  }

  matCooComplex = bfMatCooComplexNew();
  HANDLE_ERROR();

  bfMatCooComplexInitEmpty(matCooComplex, numRows, numCols, numElts);
  HANDLE_ERROR();

  bfMemCopy(rowInd, numElts, sizeof(BfSize), matCooComplex->rowInd);
  bfMemCopy(colInd, numElts, sizeof(BfSize), matCooComplex->colInd);
  bfMemCopy(value, numElts, sizeof(BfComplex), matCooComplex->value);

  BF_ERROR_END() {}

  bfMemFree(rowInd);
  bfMemFree(colInd);
  bfMemFree(value);

  return bfMatCooComplexToMat(matCooComplex);
}

BfMat *bfMatBlockDenseToType(BfMat const *mat, BfType type) {
  switch (type) {
  case BF_TYPE_MAT_COO_COMPLEX:
    return toType_cooComplex(mat);
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

BfSizeArray *bfMatBlockDenseGetNonzeroColumnRanges(BfMatBlockDense const *matBlockDense) {
  BF_ERROR_BEGIN();

  BfSizeArray *nonzeroColumnRanges = bfSizeArrayNew();
  HANDLE_ERROR();

  bfSizeArrayInitWithDefaultCapacity(nonzeroColumnRanges);
  HANDLE_ERROR();

  /* TODO: for now just assume there are *no* zero blocks. We can fix
   * this later. */
  for (BfSize p = 0; p < NUM_ROW_BLOCKS(matBlockDense); ++p)
    for (BfSize q = 0; q < NUM_COL_BLOCKS(matBlockDense); ++q)
      BF_ASSERT(bfMatGetType(BLOCK(matBlockDense, p, q)) != BF_TYPE_MAT_ZERO);

  BfMat const *mat = bfMatBlockDenseConstToMatConst(matBlockDense);

  bfSizeArrayAppend(nonzeroColumnRanges, 0);
  bfSizeArrayAppend(nonzeroColumnRanges, bfMatGetNumCols(mat));

  BF_ERROR_END() {
    BF_DIE();
  }

  return nonzeroColumnRanges;
}

void bfMatBlockDensePrintBlocksDeep(BfMatBlockDense const *matBlockDense,
                                    FILE *fp,
                                    BfSize i0, BfSize j0, BfSize depth) {
  BfMat const *mat = bfMatBlockDenseConstToMatConst(matBlockDense);

  BfSize i1 = i0 + bfMatGetNumRows(mat);
  BfSize j1 = j0 + bfMatGetNumCols(mat);

  fprintf(fp, "%u %lu %lu %lu %lu %lu\n",
          BF_TYPE_MAT_BLOCK_DENSE, i0, i1, j0, j1, depth);

  BfSize numBlocks = bfMatBlockDenseNumBlocks(matBlockDense);
  for (BfSize k = 0; k < numBlocks; ++k) {
    BfSize i0_ = i0 + BLOCK_ROW_OFFSET(matBlockDense, k);
    BfSize j0_ = j0 + BLOCK_COL_OFFSET(matBlockDense, k);
    BfMat const *block = LINEAR_BLOCK(matBlockDense, k);
    bfMatPrintBlocksDeep(block, fp, i0_, j0_, depth + 1);
  }
}

/** Interface: MatBlock */

static BfMatBlockVtable MAT_BLOCK_VTABLE = {
  .NumBlocks = (__typeof__(&bfMatBlockNumBlocks))bfMatBlockDenseNumBlocks,
  .GetNumRowBlocks = (__typeof__(&bfMatBlockGetNumRowBlocks))bfMatBlockDenseGetNumRowBlocks,
  .GetNumColBlocks = (__typeof__(&bfMatBlockGetNumColBlocks))bfMatBlockDenseGetNumColBlocks,
  .GetRowOffset = (__typeof__(&bfMatBlockGetRowOffset))bfMatBlockDenseGetRowOffset,
  .GetColOffset = (__typeof__(&bfMatBlockGetColOffset))bfMatBlockDenseGetColOffset,
  .GetBlockConst = (__typeof__(&bfMatBlockGetBlockConst))bfMatBlockDenseGetBlockConst,
};

BfSize bfMatBlockDenseNumBlocks(BfMatBlockDense const *matBlock) {
  BfMat const *mat = bfMatBlockDenseConstToMatConst(matBlock);
  return mat->numRows*mat->numCols;
}

BfSize bfMatBlockDenseGetNumRowBlocks(BfMatBlockDense const *matBlockDense) {
  BfMat const *mat = bfMatBlockDenseConstToMatConst(matBlockDense);
  if (bfMatGetType(mat) != BF_TYPE_MAT_BLOCK_DENSE) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return BF_SIZE_BAD_VALUE;
  } else {
    return mat->numRows;
  }
}

BfSize bfMatBlockDenseGetNumColBlocks(BfMatBlockDense const *matBlockDense) {
  BfMat const *mat = bfMatBlockDenseConstToMatConst(matBlockDense);
  if (bfMatGetType(mat) != BF_TYPE_MAT_BLOCK_DENSE) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return BF_SIZE_BAD_VALUE;
  } else {
    return mat->numCols;
  }
}

BfSize bfMatBlockDenseGetRowOffset(BfMatBlockDense const *matBlockDense, BfSize i) {
  return ROW_OFFSET(matBlockDense, i);
}

BfSize bfMatBlockDenseGetColOffset(BfMatBlockDense const *matBlockDense, BfSize j) {
  return COL_OFFSET(matBlockDense, j);
}

BfMat const *bfMatBlockDenseGetBlockConst(BfMatBlockDense const *mat, BfSize i, BfSize j) {
  BF_ERROR_BEGIN();

  BfMat *block = NULL;

  if (i >= NUM_ROW_BLOCKS(mat))
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  if (j >= NUM_COL_BLOCKS(mat))
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  block = mat->super.block[i*NUM_COL_BLOCKS(mat) + j];

  BF_ERROR_END() {}

  return block;
}

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
  BF_ERROR_BEGIN();

  BfMatBlockDense *matBlockDense = bfMemAlloc(1, sizeof(BfMatBlockDense));
  if (matBlockDense == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BF_ERROR_END()
    matBlockDense = NULL;

  return matBlockDense;
}

BfMatBlockDense *bfMatBlockDenseNewFromBlocks(BfSize numRowBlocks, BfSize numColBlocks, BfPtrArray *blocks, BfPolicy policy) {
  BF_ERROR_BEGIN();

  if (bfPtrArrayIsEmpty(blocks))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatBlockDense *matBlockDense = bfMatBlockDenseNew();
  HANDLE_ERROR();

  BfSize numBlocks = bfPtrArraySize(blocks);
  BF_ASSERT(numBlocks == numRowBlocks*numColBlocks);

  bfMatBlockInit(&matBlockDense->super, &MAT_VTABLE, &MAT_BLOCK_VTABLE,
                 numBlocks, numRowBlocks, numColBlocks);
  HANDLE_ERROR();

  /* Compute the row and column offsets */

  ROW_OFFSET(matBlockDense, 0) = 0;
  for (BfSize p = 0; p < numRowBlocks; ++p)
    ROW_OFFSET(matBlockDense, p + 1) = BF_SIZE_BAD_VALUE;

  COL_OFFSET(matBlockDense, 0) = 0;
  for (BfSize q = 0; q < numColBlocks; ++q)
    COL_OFFSET(matBlockDense, q + 1) = BF_SIZE_BAD_VALUE;

  for (BfSize p = 0; p < numRowBlocks; ++p) {
    for (BfSize q = 0; q < numColBlocks; ++q) {
      BfSize k = numColBlocks*p + q;

      BfMat *block = bfPtrArrayGet(blocks, k);

      BLOCK(matBlockDense, p, q) = bfMatGet(block, policy);

      BfSize numRows = bfMatGetNumRows(block);
      BfSize numCols = bfMatGetNumCols(block);

      if (ROW_OFFSET(matBlockDense, p + 1) != BF_SIZE_BAD_VALUE)
        BF_ASSERT(ROW_OFFSET(matBlockDense, p + 1) == numRows);
      ROW_OFFSET(matBlockDense, p + 1) = numRows;

      if (COL_OFFSET(matBlockDense, q + 1) != BF_SIZE_BAD_VALUE)
        BF_ASSERT(COL_OFFSET(matBlockDense, q + 1) == numCols);
      COL_OFFSET(matBlockDense, q + 1) = numCols;
    }
  }

  bfSizeRunningSum(numRowBlocks + 1, &ROW_OFFSET(matBlockDense, 0));
  bfSizeRunningSum(numColBlocks + 1, &COL_OFFSET(matBlockDense, 0));

  BF_ERROR_END() {
    BF_DIE();
  }

  return matBlockDense;
}

BfMatBlockDense *bfMatBlockDenseNewColFromBlocks(BfPtrArray *blocks, BfPolicy policy) {
  BF_ERROR_BEGIN();

  if (bfPtrArrayIsEmpty(blocks))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatBlockDense *matBlockDense = bfMatBlockDenseNew();
  HANDLE_ERROR();

  BfSize numBlocks = bfPtrArraySize(blocks);

  bfMatBlockInit(&matBlockDense->super, &MAT_VTABLE, &MAT_BLOCK_VTABLE,
                 numBlocks, numBlocks, 1);
  HANDLE_ERROR();

  /* Compute the row and column offsets */

  BfMat *block = bfPtrArrayGet(blocks, 0);

  BfSize n = bfMatGetNumCols(block);

  COL_OFFSET(matBlockDense, 0) = 0;
  COL_OFFSET(matBlockDense, 1) = n;

  ROW_OFFSET(matBlockDense, 0) = 0;
  for (BfSize k = 0; k < numBlocks; ++k) {
    block = bfPtrArrayGet(blocks, k);
    BLOCK(matBlockDense, k, 0) = bfMatGet(block, policy);
    ROW_OFFSET(matBlockDense, k + 1) = bfMatGetNumRows(block);
  }
  bfSizeRunningSum(numBlocks + 1, &ROW_OFFSET(matBlockDense, 0));

  BF_ERROR_END() {
    BF_DIE();
  }

  return matBlockDense;
}

BfMatBlockDense *bfMatBlockDenseNewRowFromBlocks(BfPtrArray *blocks, BfPolicy policy) {
  BF_ERROR_BEGIN();

  if (bfPtrArrayIsEmpty(blocks))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatBlockDense *matBlockDense = bfMatBlockDenseNew();
  HANDLE_ERROR();

  BfSize numBlocks = bfPtrArraySize(blocks);

  bfMatBlockInit(&matBlockDense->super, &MAT_VTABLE, &MAT_BLOCK_VTABLE,
                 numBlocks, 1, numBlocks);
  HANDLE_ERROR();

  /* Compute the row and column offsets */

  BfMat *block = bfPtrArrayGet(blocks, 0);

  BfSize m = bfMatGetNumRows(block);

  ROW_OFFSET(matBlockDense, 0) = 0;
  ROW_OFFSET(matBlockDense, 1) = m;

  COL_OFFSET(matBlockDense, 0) = 0;
  for (BfSize k = 0; k < numBlocks; ++k) {
    block = bfPtrArrayGet(blocks, k);
    BLOCK(matBlockDense, 0, k) = bfMatGet(block, policy);
    COL_OFFSET(matBlockDense, k + 1) = bfMatGetNumCols(block);
  }
  bfSizeRunningSum(numBlocks + 1, &COL_OFFSET(matBlockDense, 0));

  BF_ERROR_END() {
    BF_DIE();
  }

  return matBlockDense;
}

void bfMatBlockDenseInit(BfMatBlockDense *mat,
                         BfSize numRowBlocks, BfSize numColBlocks) {
  BF_ERROR_BEGIN();

  BfSize numBlocks = numRowBlocks*numColBlocks;

  bfMatBlockInit(&mat->super,
                 &MAT_VTABLE, &MAT_BLOCK_VTABLE,
                 numBlocks, numRowBlocks, numColBlocks);
  HANDLE_ERROR();

  BF_ERROR_END()
    bfMatBlockDeinit(&mat->super);
}

void bfMatBlockDenseDeinit(BfMatBlockDense *mat) {
  bfMatBlockDeinit(&mat->super);
}

void bfMatBlockDenseDealloc(BfMatBlockDense **mat) {
  bfMemFree(*mat);
  *mat = NULL;
}

void bfMatBlockDenseDeinitAndDealloc(BfMatBlockDense **mat) {
  bfMatBlockDenseDeinit(*mat);
  bfMatBlockDenseDealloc(mat);
}

BfMat *bfMatBlockDenseGetBlock(BfMatBlockDense *mat, BfSize i, BfSize j) {
  BF_ERROR_BEGIN();

  BfMat *block = NULL;

  if (i >= NUM_ROW_BLOCKS(mat))
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  if (j >= NUM_COL_BLOCKS(mat))
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  block = mat->super.block[i*NUM_COL_BLOCKS(mat) + j];

  BF_ERROR_END() {}

  return block;
}

void bfMatBlockDenseSetBlock(BfMatBlockDense *mat, BfSize i, BfSize j,
                             BfMat *block) {
  BF_ERROR_BEGIN();

  if (i >= NUM_ROW_BLOCKS(mat))
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  if (j >= NUM_COL_BLOCKS(mat))
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  if (bfMatGetNumRows(block) != NUM_BLOCK_ROWS(mat, i))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (bfMatGetNumCols(block) != NUM_BLOCK_COLS(mat, j))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BLOCK(mat, i, j) = block;

  BF_ERROR_END() {}
}

void bfMatBlockDenseAppendColumn(BfMatBlockDense *mat, BfSize numBlockCols) {
  BF_ERROR_BEGIN();

  BfSize oldNumColBlocks = BF_SIZE_BAD_VALUE;

  BfMatBlock *matBlock = NULL;

  matBlock = bfMatBlockDenseToMatBlock(mat);
  HANDLE_ERROR();

  /** Update the number of block columns and their offsets: */

  oldNumColBlocks = NUM_COL_BLOCKS(mat)++;

  /* Allocate more space for column offsets */
  BfSize *newColOffset = bfMemRealloc(
    matBlock->colOffset, NUM_COL_BLOCKS(mat) + 1, sizeof(BfSize));
  if (newColOffset == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  matBlock->colOffset = newColOffset;

  /* Compute new final column offset */
  matBlock->colOffset[NUM_COL_BLOCKS(mat)] =
    matBlock->colOffset[oldNumColBlocks] + numBlockCols;

  /** Allocate space for a new column of blocks, arrange them
   * correctly, and insert zero blocks */

  /* Allocate more space for the new block column */
  BfMat **newBlockPtr = bfMemRealloc(
    mat->super.block, NUM_BLOCKS(mat), sizeof(BfMat **));
  if (newBlockPtr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  mat->super.block = newBlockPtr;

  /* Adjust the position of each of the old rows */
  for (BfSize i = NUM_ROW_BLOCKS(mat); i > 0; --i)
    bfMemMove(mat->super.block + (i - 1)*oldNumColBlocks,
              oldNumColBlocks, sizeof(BfMat *),
              mat->super.block + (i - 1)*NUM_COL_BLOCKS(mat));

  /* Set new blocks to zero */
  for (BfSize i = 0; i < NUM_ROW_BLOCKS(mat); ++i) {
    BfMatZero *block = bfMatZeroNew();
    HANDLE_ERROR();

    bfMatZeroInit(block, NUM_BLOCK_ROWS(mat, i), numBlockCols);
    HANDLE_ERROR();

    bfMatBlockDenseSetBlock(
      mat, i, NUM_COL_BLOCKS(mat) - 1, bfMatZeroToMat(block));
  }

  BF_ERROR_END() {
    /* Reset the number of block columns */
    if (NUM_COL_BLOCKS(mat) != oldNumColBlocks)
      NUM_COL_BLOCKS(mat) = oldNumColBlocks;
  }
}

void bfMatBlockDenseAppendRow(BfMatBlockDense *mat, BfSize numBlockRows) {
  BF_ERROR_BEGIN();

  BfMatBlock *matBlock = NULL;

  matBlock = bfMatBlockDenseToMatBlock(mat);
  HANDLE_ERROR();

  /** Update the number of block rows and their offsets: */

  BfSize oldNumRowBlocks = NUM_ROW_BLOCKS(mat)++;

  /* Allocate more space for row offsets */
  BfSize *newRowOffset = bfMemRealloc(
    matBlock->rowOffset, NUM_ROW_BLOCKS(mat) + 1, sizeof(BfSize));
  if (newRowOffset == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  matBlock->rowOffset = newRowOffset;

  /* Compute new final row offset */
  matBlock->rowOffset[NUM_ROW_BLOCKS(mat)] =
    matBlock->rowOffset[oldNumRowBlocks] + numBlockRows;

  /** Allocate space for a new row of blocks, insert zeros */

  /* Allocate more space for the new block row */
  BfMat **newBlockPtr = bfMemRealloc(
    mat->super.block, NUM_BLOCKS(mat), sizeof(BfMat **));
  if (newBlockPtr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  mat->super.block = newBlockPtr;

  /* Set new blocks to zero */
  for (BfSize j = 0; j < NUM_COL_BLOCKS(mat); ++j) {
    BfMatZero *block = bfMatZeroNew();
    HANDLE_ERROR();

    bfMatZeroInit(block, numBlockRows, NUM_BLOCK_COLS(mat, j));
    HANDLE_ERROR();

    bfMatBlockDenseSetBlock(
      mat, NUM_ROW_BLOCKS(mat) - 1, j, bfMatZeroToMat(block));
  }

  BF_ERROR_END() {}
}
