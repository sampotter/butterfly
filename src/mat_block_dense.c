#include <bf/mat_block_dense.h>

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_coo_complex.h>
#include <bf/mat_sum.h>
#include <bf/mat_zero.h>
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
  .GetRowCopy = (__typeof__(&bfMatBlockDenseGetRowCopy))bfMatBlockDenseGetRowCopy,
  .Delete = (__typeof__(&bfMatBlockDenseDelete))bfMatBlockDenseDelete,
  .GetType = (__typeof__(&bfMatBlockDenseGetType))bfMatBlockDenseGetType,
  .NumBytes = (__typeof__(&bfMatBlockDenseNumBytes))bfMatBlockDenseNumBytes,
  .GetNumRows = (__typeof__(&bfMatBlockDenseGetNumRows))bfMatBlockDenseGetNumRows,
  .GetNumCols = (__typeof__(&bfMatBlockDenseGetNumCols))bfMatBlockDenseGetNumCols,
  .GetRowRange = (__typeof__(&bfMatBlockDenseGetRowRange))bfMatBlockDenseGetRowRange,
  .GetRowRangeCopy = (__typeof__(&bfMatGetRowRangeCopy))bfMatBlockDenseGetRowRangeCopy,
  .ScaleCols = (__typeof__(&bfMatScaleCols))bfMatBlockDenseScaleCols,
  .AddInplace = (__typeof__(&bfMatAddInplace))bfMatBlockDenseAddInplace,
  .Mul = (__typeof__(&bfMatBlockDenseMul))bfMatBlockDenseMul,
  .MulVec = (__typeof__(&bfMatMulVec))bfMatBlockDenseMulVec,
  .ToType = (__typeof__(&bfMatBlockDenseToType))bfMatBlockDenseToType,
  .PrintBlocksDeep = (__typeof__(&bfMatPrintBlocksDeep))bfMatBlockDensePrintBlocksDeep
};

BfMat *bfMatBlockDenseCopy(BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

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

  END_ERROR_HANDLING() {
    bfMatDelete(&blockCopy);
    bfMatBlockDenseDeinitAndDealloc(&matBlockDenseCopy);
  }

  return bfMatBlockDenseToMat(matBlockDenseCopy);
}

BfVec *bfMatBlockDenseGetRowCopy(BfMat const *mat, BfSize i) {
  BEGIN_ERROR_HANDLING();

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
  BfSize numColBlocks = bfMatBlockGetNumColBlocks(matBlock);

  BfSize ib = 0, i0, i1;
  for (; ib < numRowBlocks; ++ib) {
    i0 = ROW_OFFSET(matBlockDense, ib);
    i1 = ROW_OFFSET(matBlockDense, ib + 1);
    if (i0 <= i && i < i1)
      break;
  }

  for (BfSize jb = 0; jb < numColBlocks; ++jb) {
    BfMat *block = BLOCK(matBlockDense, ib, jb);
    BfVec *blockRowCopy = bfMatGetRowCopy(block, i - i0);
    assert(blockRowCopy->size == NUM_BLOCK_COLS(matBlockDense, jb));
    if (rowCopy == NULL) {
      rowCopy = blockRowCopy;
    } else {
      BfVec *cat = bfVecConcat(rowCopy, blockRowCopy);
      HANDLE_ERROR();
      bfVecDelete(&rowCopy);
      bfVecDelete(&blockRowCopy);
      rowCopy = cat;
    }
    assert(rowCopy->size == matBlockDense->super.colOffset[jb + 1]);
  }

  assert(rowCopy->size == bfMatGetNumCols(mat));

  END_ERROR_HANDLING() {}

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

BfMat *bfMatBlockDenseGetRowRangeCopy(BfMatBlockDense const *matBlockDense,
                                      BfSize i0, BfSize i1) {
  BEGIN_ERROR_HANDLING();

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
  while (p1 > 0 && ROW_OFFSET(matBlockDense, p1 - 1) > i0)
    --p1;

  assert(p0 < p1);
  assert(ROW_OFFSET(matBlockDense, p0) == i0);
  assert(ROW_OFFSET(matBlockDense, p1) == i1);

  /** Copy all blocks in the row range and assemble the resulting
   * dense block matrix containing the row range. */

  /* Array holding all `(p1 - p0)*NUM_COL_BLOCKS(matBlockDense)`
   * blocks in row major order. */
  BfPtrArray *blocks = bfPtrArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  for (BfSize p = p0; p < p1; ++p) {
    for (BfSize q = 0; q < NUM_COL_BLOCKS(matBlockDense); ++q) {
      BfMat const *block = BLOCK(matBlockDense, p, q);

      BfMat *blockCopy = bfMatCopy(block);
      HANDLE_ERROR();

      bfPtrArrayAppend(blocks, blockCopy);
      HANDLE_ERROR();
    }
  }

  BfMatBlockDense *result = bfMatBlockDenseNewFromBlocks(
    p1 - p0, NUM_COL_BLOCKS(matBlockDense), blocks);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    assert(false);
  }

  bfPtrArrayDelete(&blocks);

  return bfMatBlockDenseToMat(result);
}

void bfMatBlockDenseScaleCols(BfMatBlockDense *matBlockDense, BfVec const *vec) {
  BEGIN_ERROR_HANDLING();

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

  END_ERROR_HANDLING() {}
}

void bfMatBlockDenseAddInplace(BfMatBlockDense *matBlockDense, BfMat const *otherMat) {
  BEGIN_ERROR_HANDLING();

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

  END_ERROR_HANDLING() {}
}

BfMat *bfMatBlockDenseMul(BfMat const *mat, BfMat const *otherMat) {
  BEGIN_ERROR_HANDLING();

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

BfVec *bfMatBlockDenseMulVec(BfMatBlockDense const *matBlockDense,
                             BfVec const *vec) {
  BEGIN_ERROR_HANDLING();

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

    BfVec *resultSubvecView = bfVecGetSubvecView(result, i0, i1);
    HANDLE_ERROR();

    for (BfSize j = 0; j < numColBlocks; ++j) {
      BfMat const *block = BLOCK(matBlockDense, i, j);

      BfSize j0 = COL_OFFSET(matBlockDense, j);
      BfSize j1 = COL_OFFSET(matBlockDense, j + 1);

      /* Sanity check... basically unnecessary, but whatever */
      assert(bfMatGetNumRows(block) == i1 - i0);
      assert(bfMatGetNumCols(block) == j1 - j0);

      BfVec const *subvecView = bfVecGetSubvecViewConst(vec, j0, j1);
      HANDLE_ERROR();

      BfVec *tmp = bfMatMulVec(block, subvecView);
      HANDLE_ERROR();

      bfVecAddInplace(resultSubvecView, tmp);

      bfVecDelete(&tmp);
    }

    bfVecDelete(&resultSubvecView);
  }

  END_ERROR_HANDLING() {
    assert(false);
  }

  return result;
}

static BfMat *toType_cooComplex(BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfMatCooComplex *matCooComplex = NULL;

  BfSize numRows = bfMatGetNumRows(mat);
  BfSize numCols = bfMatGetNumCols(mat);

  BfSize *rowInd = NULL;
  BfSize *colInd = NULL;
  BfComplex *value = NULL;

  BfSize numElts = 0, capacity = 16;

  rowInd = malloc(capacity*sizeof(BfSize));
  if (rowInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  colInd = malloc(capacity*sizeof(BfSize));
  if (colInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  value = malloc(capacity*sizeof(BfComplex));
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

        BfSize *newRowInd = realloc(rowInd, capacity*sizeof(BfSize));
        if (newRowInd == NULL)
          RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
        rowInd = newRowInd;

        BfSize *newColInd = realloc(colInd, capacity*sizeof(BfSize));
        if (newColInd == NULL)
          RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
        colInd = newColInd;

        BfComplex *newValue = realloc(value, capacity*sizeof(BfComplex));
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

  memcpy(matCooComplex->rowInd, rowInd, numElts*sizeof(BfSize));
  memcpy(matCooComplex->colInd, colInd, numElts*sizeof(BfSize));
  memcpy(matCooComplex->value, value, numElts*sizeof(BfComplex));

  END_ERROR_HANDLING() {}

  free(rowInd);
  free(colInd);
  free(value);

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
  BEGIN_ERROR_HANDLING();

  BfMatBlockDense *matBlockDense = malloc(sizeof(BfMatBlockDense));
  if (matBlockDense == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING()
    matBlockDense = NULL;

  return matBlockDense;
}

BfMatBlockDense *bfMatBlockDenseNewFromBlocks(BfSize numRowBlocks,
                                              BfSize numColBlocks,
                                              BfPtrArray *blocks) {
  BEGIN_ERROR_HANDLING();

  if (bfPtrArrayIsEmpty(blocks))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatBlockDense *matBlockDense = bfMatBlockDenseNew();
  HANDLE_ERROR();

  BfSize numBlocks = bfPtrArraySize(blocks);
  assert(numBlocks == numRowBlocks*numColBlocks);

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

      BLOCK(matBlockDense, p, q) = block;

      BfSize numRows = bfMatGetNumRows(block);
      BfSize numCols = bfMatGetNumCols(block);

      if (ROW_OFFSET(matBlockDense, p + 1) != BF_SIZE_BAD_VALUE)
        assert(ROW_OFFSET(matBlockDense, p + 1) == numRows);
      ROW_OFFSET(matBlockDense, p + 1) = numRows;

      if (COL_OFFSET(matBlockDense, q + 1) != BF_SIZE_BAD_VALUE)
        assert(COL_OFFSET(matBlockDense, q + 1) == numCols);
      COL_OFFSET(matBlockDense, q + 1) = numCols;
    }
  }

  bfSizeRunningSum(numRowBlocks + 1, &ROW_OFFSET(matBlockDense, 0));
  bfSizeRunningSum(numColBlocks + 1, &COL_OFFSET(matBlockDense, 0));

  END_ERROR_HANDLING() {
    assert(false);
  }

  return matBlockDense;
}

BfMatBlockDense *bfMatBlockDenseNewColFromBlocks(BfPtrArray *blocks) {
  BEGIN_ERROR_HANDLING();

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
    BLOCK(matBlockDense, k, 0) = block;
    ROW_OFFSET(matBlockDense, k + 1) = bfMatGetNumRows(block);
  }
  bfSizeRunningSum(numBlocks + 1, &ROW_OFFSET(matBlockDense, 0));

  END_ERROR_HANDLING() {
    assert(false);
  }

  return matBlockDense;
}

BfMatBlockDense *bfMatBlockDenseNewRowFromBlocks(BfPtrArray *blocks) {
  BEGIN_ERROR_HANDLING();

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
    BLOCK(matBlockDense, 0, k) = block;
    COL_OFFSET(matBlockDense, k + 1) = bfMatGetNumCols(block);
  }
  bfSizeRunningSum(numBlocks + 1, &COL_OFFSET(matBlockDense, 0));

  END_ERROR_HANDLING() {
    assert(false);
  }

  return matBlockDense;
}

void bfMatBlockDenseInit(BfMatBlockDense *mat,
                         BfSize numRowBlocks, BfSize numColBlocks) {
  BEGIN_ERROR_HANDLING();

  BfSize numBlocks = numRowBlocks*numColBlocks;

  bfMatBlockInit(&mat->super,
                 &MAT_VTABLE, &MAT_BLOCK_VTABLE,
                 numBlocks, numRowBlocks, numColBlocks);
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

  if (i >= NUM_ROW_BLOCKS(mat))
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  if (j >= NUM_COL_BLOCKS(mat))
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  block = mat->super.block[i*NUM_COL_BLOCKS(mat) + j];

  END_ERROR_HANDLING() {}

  return block;
}

BfMat const *bfMatBlockDenseGetBlockConst(BfMatBlockDense const *mat, BfSize i, BfSize j) {
  BEGIN_ERROR_HANDLING();

  BfMat *block = NULL;

  if (i >= NUM_ROW_BLOCKS(mat))
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  if (j >= NUM_COL_BLOCKS(mat))
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  block = mat->super.block[i*NUM_COL_BLOCKS(mat) + j];

  END_ERROR_HANDLING() {}

  return block;
}

void bfMatBlockDenseSetBlock(BfMatBlockDense *mat, BfSize i, BfSize j,
                             BfMat *block) {
  BEGIN_ERROR_HANDLING();

  if (i >= NUM_ROW_BLOCKS(mat))
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  if (j >= NUM_COL_BLOCKS(mat))
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  if (bfMatGetNumRows(block) != NUM_BLOCK_ROWS(mat, i))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (bfMatGetNumCols(block) != NUM_BLOCK_COLS(mat, j))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BLOCK(mat, i, j) = block;

  END_ERROR_HANDLING() {}
}

void bfMatBlockDenseAppendColumn(BfMatBlockDense *mat, BfSize numBlockCols) {
  BEGIN_ERROR_HANDLING();

  BfMatBlock *matBlock = NULL;

  matBlock = bfMatBlockDenseToMatBlock(mat);
  HANDLE_ERROR();

  /** Update the number of block columns and their offsets: */

  BfSize oldNumColBlocks = NUM_COL_BLOCKS(mat)++;

  /* Allocate more space for column offsets */
  BfSize *newColOffset = realloc(
    matBlock->colOffset, (NUM_COL_BLOCKS(mat) + 1)*sizeof(BfSize));
  if (newColOffset == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  matBlock->colOffset = newColOffset;

  /* Compute new final column offset */
  matBlock->colOffset[NUM_COL_BLOCKS(mat)] =
    matBlock->colOffset[oldNumColBlocks] + numBlockCols;

  /** Allocate space for a new column of blocks, arrange them
   * correctly, and insert zero blocks */

  /* Allocate more space for the new block column */
  BfMat **newBlockPtr = realloc(
    mat->super.block, NUM_BLOCKS(mat)*sizeof(BfMat **));
  if (newBlockPtr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  mat->super.block = newBlockPtr;

  /* Adjust the position of each of the old rows */
  for (BfSize i = NUM_ROW_BLOCKS(mat); i > 0; --i)
    memmove(mat->super.block + (i - 1)*NUM_COL_BLOCKS(mat),
            mat->super.block + (i - 1)*oldNumColBlocks,
            oldNumColBlocks*sizeof(BfMat *));

  /* Set new blocks to zero */
  for (BfSize i = 0; i < NUM_ROW_BLOCKS(mat); ++i) {
    BfMatZero *block = bfMatZeroNew();
    HANDLE_ERROR();

    bfMatZeroInit(block, NUM_BLOCK_ROWS(mat, i), numBlockCols);
    HANDLE_ERROR();

    bfMatBlockDenseSetBlock(
      mat, i, NUM_COL_BLOCKS(mat) - 1, bfMatZeroToMat(block));
  }

  END_ERROR_HANDLING() {
    /* Reset the number of block columns */
    if (NUM_COL_BLOCKS(mat) != oldNumColBlocks)
      NUM_COL_BLOCKS(mat) = oldNumColBlocks;
  }
}

void bfMatBlockDenseAppendRow(BfMatBlockDense *mat, BfSize numBlockRows) {
  BEGIN_ERROR_HANDLING();

  BfMatBlock *matBlock = NULL;

  matBlock = bfMatBlockDenseToMatBlock(mat);
  HANDLE_ERROR();

  /** Update the number of block rows and their offsets: */

  BfSize oldNumRowBlocks = NUM_ROW_BLOCKS(mat)++;

  /* Allocate more space for row offsets */
  BfSize *newRowOffset = realloc(
    matBlock->rowOffset, (NUM_ROW_BLOCKS(mat) + 1)*sizeof(BfSize));
  if (newRowOffset == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  matBlock->rowOffset = newRowOffset;

  /* Compute new final row offset */
  matBlock->rowOffset[NUM_ROW_BLOCKS(mat)] =
    matBlock->rowOffset[oldNumRowBlocks] + numBlockRows;

  /** Allocate space for a new row of blocks, insert zeros */

  /* Allocate more space for the new block row */
  BfMat **newBlockPtr = realloc(
    mat->super.block, NUM_BLOCKS(mat)*sizeof(BfMat **));
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

  END_ERROR_HANDLING() {}
}
