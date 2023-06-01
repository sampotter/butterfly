#include <bf/mat_block_coo.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/indexed_mat.h>
#include <bf/mat_zero.h>
#include <bf/mem.h>
#include <bf/ptr_array.h>
#include <bf/size_array.h>
#include <bf/util.h>
#include <bf/vec_real.h>

/** Helper macros: */

#define NUM_ROW_BLOCKS(mat) mat->super.super.numRows
#define NUM_COL_BLOCKS(mat) mat->super.super.numCols
#define ROW_IND(mat, k) mat->rowInd[k]
#define COL_IND(mat, k) mat->colInd[k]
#define ROW_OFFSET(mat, i) mat->super.rowOffset[i]
#define COL_OFFSET(mat, j) mat->super.colOffset[j]
#define BLOCK_ROW_OFFSET(mat, k) ROW_OFFSET(mat, mat->rowInd[k])
#define BLOCK_COL_OFFSET(mat, k) COL_OFFSET(mat, mat->colInd[k])
#define BLOCK(mat, k) mat->super.block[k]

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .GetView = (__typeof__(&bfMatGetView))bfMatBlockCooGetView,
  .Copy = (__typeof__(&bfMatBlockCooCopy))bfMatBlockCooCopy,
  .Steal = (__typeof__(&bfMatSteal))bfMatBlockCooSteal,
  .GetRowCopy = (__typeof__(&bfMatBlockCooGetRowCopy))bfMatBlockCooGetRowCopy,
  .Delete = (__typeof__(&bfMatBlockCooDelete))bfMatBlockCooDelete,
  .GetType = (__typeof__(&bfMatBlockCooGetType))bfMatBlockCooGetType,
  .NumBytes = (__typeof__(&bfMatBlockCooNumBytes))bfMatBlockCooNumBytes,
  .Dump = (__typeof__(&bfMatDump))bfMatBlockCooDump,
  .GetNumRows = (__typeof__(&bfMatBlockCooGetNumRows))bfMatBlockCooGetNumRows,
  .GetNumCols = (__typeof__(&bfMatBlockCooGetNumCols))bfMatBlockCooGetNumCols,
  .GetRowRangeCopy = (__typeof__(&bfMatGetRowRangeCopy))bfMatBlockCooGetRowRangeCopy,
  .Mul = (__typeof__(&bfMatBlockCooMul))bfMatBlockCooMul,
  .MulVec = (__typeof__(&bfMatMulVec))bfMatBlockCooMulVec,
  .RmulVec = (__typeof__(&bfMatRmulVec))bfMatBlockCooRmulVec,
  .Negate = (__typeof__(&bfMatBlockCooNegate))bfMatBlockCooNegate,
  .GetNonzeroColumnRanges = (__typeof__(&bfMatGetNonzeroColumnRanges))bfMatBlockCooGetNonzeroColumnRanges,
  .PrintBlocksDeep = (__typeof__(&bfMatPrintBlocksDeep))bfMatBlockCooPrintBlocksDeep,
};

BfMat *bfMatBlockCooGetView(BfMatBlockCoo const *matBlockCoo) {
  BF_ERROR_BEGIN();

  BfMatBlockCoo *matBlockCooView = bfMatBlockCooNew();
  HANDLE_ERROR();

  *matBlockCooView = *matBlockCoo;

  BfMat *matView = bfMatBlockCooToMat(matBlockCooView);

  matView->props |= BF_MAT_PROPS_VIEW;

  BF_ERROR_END() {
    BF_DIE();
  }

  return matView;
}

BfMat *bfMatBlockCooCopy(BfMat const *mat) {
  BF_ERROR_BEGIN();

  BfMatBlock const *matBlock = NULL;
  BfMatBlockCoo const *matBlockCoo = NULL;
  BfMatBlockCoo *matBlockCooCopy = NULL;
  BfMat *blockCopy = NULL;

  matBlock = bfMatConstToMatBlockConst(mat);
  HANDLE_ERROR();

  matBlockCoo = bfMatConstToMatBlockCooConst(mat);
  HANDLE_ERROR();

  BfSize numRowBlocks = bfMatBlockGetNumRowBlocks(matBlock);
  BfSize numColBlocks = bfMatBlockGetNumColBlocks(matBlock);
  BfSize numBlocks = bfMatBlockNumBlocks(matBlock);

  matBlockCooCopy = bfMatBlockCooNew();
  HANDLE_ERROR();

  bfMatBlockCooInit(matBlockCooCopy, numRowBlocks, numColBlocks, numBlocks);
  HANDLE_ERROR();

  for (BfSize k = 0; k < numBlocks; ++k) {
    matBlockCooCopy->rowInd[k] = matBlockCoo->rowInd[k];
    matBlockCooCopy->colInd[k] = matBlockCoo->colInd[k];
  }

  for (BfSize i = 0; i <= numRowBlocks; ++i)
    matBlockCooCopy->super.rowOffset[i] = matBlock->rowOffset[i];

  for (BfSize j = 0; j <= numColBlocks; ++j)
    matBlockCooCopy->super.colOffset[j] = matBlock->colOffset[j];

  for (BfSize k = 0; k < numBlocks; ++k) {
    blockCopy = bfMatCopy(matBlock->block[k]);
    HANDLE_ERROR();
    matBlockCooCopy->super.block[k] = blockCopy;
  }

  BF_ERROR_END() {
    bfMatDelete(&blockCopy);
    bfMatBlockCooDeinitAndDealloc(&matBlockCooCopy);
  }

  return bfMatBlockCooToMat(matBlockCooCopy);
}

BfMat *bfMatBlockCooSteal(BfMatBlockCoo *matBlockCoo) {
  BF_ERROR_BEGIN();

  BfMat *mat = bfMatBlockCooToMat(matBlockCoo);

  if (bfMatIsView(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatBlockCoo *matBlockCooNew = bfMatBlockCooNew();
  HANDLE_ERROR();

  *matBlockCooNew = *matBlockCoo;

  mat->props |= BF_MAT_PROPS_VIEW;

  BF_ERROR_END() {
    BF_DIE();
  }

  return bfMatBlockCooToMat(matBlockCooNew);
}

typedef struct {
  BfSize i0, j0;
  BfMat const *block;
} BfMatBlockCooEntry;

static int cmpColumnOrder(BfMatBlockCooEntry const **arg1,
                          BfMatBlockCooEntry const **arg2) {
  int64_t arg1_j0 = (*arg1)->j0;
  int64_t arg2_j0 = (*arg2)->j0;
  return arg1_j0 - arg2_j0;
}

BfVec *bfMatBlockCooGetRowCopy(BfMat const *mat, BfSize i) {
  BF_ERROR_BEGIN();

  BfVec *rowCopy = NULL;
  BfMatBlock const *matBlock = NULL;
  BfMatBlockCoo const *matBlockCoo = NULL;

  matBlock = bfMatConstToMatBlockConst(mat);
  HANDLE_ERROR();

  matBlockCoo = bfMatConstToMatBlockCooConst(mat);
  HANDLE_ERROR();

  BfSize numRows = bfMatGetNumRows(mat);
  if (i >= numRows)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numBlocks = bfMatBlockNumBlocks(matBlock);
  BfSize numRowBlocks = bfMatBlockGetNumRowBlocks(matBlock);

  BfSize rowInd = 0;
  for (; rowInd < numRowBlocks; ++rowInd) {
    BfSize i0 = ROW_OFFSET(matBlockCoo, rowInd);
    BfSize i1 = ROW_OFFSET(matBlockCoo, rowInd + 1);
    if (i0 <= i && i < i1)
      break;
  }

  BfPtrArray blocks;
  bfInitPtrArrayWithDefaultCapacity(&blocks);

  for (BfSize k = 0; k < numBlocks; ++k) {
    if (ROW_IND(matBlockCoo, k) != rowInd)
      continue;
    BfMatBlockCooEntry *entry = bfMemAlloc(1, sizeof(BfMatBlockCooEntry));
    entry->i0 = BLOCK_ROW_OFFSET(matBlockCoo, k);
    entry->j0 = BLOCK_COL_OFFSET(matBlockCoo, k);
    entry->block = BLOCK(matBlockCoo, k);
    bfPtrArrayAppend(&blocks, entry);
  }

  /* Sort the array of blocks into dictionary order */
  bfPtrArraySort(&blocks, (BfPtrCmp)cmpColumnOrder, NULL);

  /* Verify that none of the collected blocks are overlapping */
#if BF_DEBUG
  BfSize j1 = 0;
  for (BfSize k = 0; k < bfPtrArraySize(&blocks); ++k) {
    BfMatBlockCooEntry const *entry = bfPtrArrayGet(&blocks, k);
    BF_ASSERT(j1 <= entry->j0);
    j1 = entry->j0;
  }
#endif

  for (BfSize k = 0; k < bfPtrArraySize(&blocks); ++k) {
    BfMatBlockCooEntry const *entry = bfPtrArrayGet(&blocks, k);
    BfVec *blockRowCopy = bfMatGetRowCopy(entry->block, i - entry->i0);
    if (rowCopy == NULL) {
      rowCopy = blockRowCopy;
    } else {
      BfVec *cat = bfVecConcat(rowCopy, blockRowCopy);
      HANDLE_ERROR();
      bfVecDelete(&rowCopy);
      bfVecDelete(&blockRowCopy);
      rowCopy = cat;
    }
  }

  BF_ASSERT(rowCopy->size == bfMatGetNumCols(mat));

  BF_ERROR_END() {}

  for (BfSize k = 0; k < bfPtrArraySize(&blocks); ++k)
    bfMemFree(bfPtrArrayGet(&blocks, k));
  bfPtrArrayDeinit(&blocks);

  return rowCopy;
}

void bfMatBlockCooDelete(BfMat **mat) {
  bfMatBlockCooDeinitAndDealloc((BfMatBlockCoo **)mat);
}

BfType bfMatBlockCooGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_BLOCK_COO;
}

BfSize bfMatBlockCooNumBytes(BfMat const *mat) {
  BfMatBlockCoo const *matBlockCoo = bfMatConstToMatBlockCooConst(mat);

  BfSize num_bytes = 0;

  /* memory occupied by mat itself */
  num_bytes += sizeof(BfMatBlockCoo);

  /* memory occupied by blocks */
  for (BfSize i = 0; i < matBlockCoo->numBlocks; ++i)
    num_bytes += bfMatNumBytes(matBlockCoo->super.block[i]);

  /* mat->super.rowOffset and mat->super.colOffset */
  num_bytes += (mat->numRows + 1)*sizeof(BfSize);
  num_bytes += (mat->numCols + 1)*sizeof(BfSize);

  /* mat->rowInd and mat->colInd */
  num_bytes += 2*matBlockCoo->numBlocks*sizeof(BfSize);

  return num_bytes;
}

void bfMatBlockCooDump(BfMatBlockCoo const *matBlockCoo, FILE *fp) {
  BF_ERROR_BEGIN();

  /* Write number of row and column blocks: */
  fwrite(&NUM_ROW_BLOCKS(matBlockCoo), sizeof(BfSize), 1, fp);
  fwrite(&NUM_COL_BLOCKS(matBlockCoo), sizeof(BfSize), 1, fp);

  /* Write number of blocks: */
  BfSize numBlocks = bfMatBlockCooNumBlocks(matBlockCoo);
  fwrite(&numBlocks, sizeof(BfSize), 1, fp);

  /* Write the row and column indices: */
  fwrite(&ROW_IND(matBlockCoo, 0), sizeof(BfSize), numBlocks, fp);
  fwrite(&COL_IND(matBlockCoo, 0), sizeof(BfSize), numBlocks, fp);

  /* Write the row and column offsets: */
  fwrite(&ROW_OFFSET(matBlockCoo, 0), sizeof(BfSize), NUM_ROW_BLOCKS(matBlockCoo) + 1, fp);
  fwrite(&COL_OFFSET(matBlockCoo, 0), sizeof(BfSize), NUM_COL_BLOCKS(matBlockCoo) + 1, fp);

  /* Recursively dump the blocks: */
  for (BfSize k = 0; k < numBlocks; ++k) {
    BfMat const *mat = BLOCK(matBlockCoo, k);
    bfMatDump(mat, fp);
    HANDLE_ERROR();
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

BfSize bfMatBlockCooGetNumRows(BfMat const *mat) {
  BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
  return bfMatIsTransposed(mat) ?
    matBlock->colOffset[bfMatBlockGetNumColBlocks(matBlock)] :
    matBlock->rowOffset[bfMatBlockGetNumRowBlocks(matBlock)];
}

BfSize bfMatBlockCooGetNumCols(BfMat const *mat) {
  BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
  return bfMatIsTransposed(mat) ?
    matBlock->rowOffset[bfMatBlockGetNumRowBlocks(matBlock)] :
    matBlock->colOffset[bfMatBlockGetNumColBlocks(matBlock)];
}

BfMat *bfMatBlockCooGetRowRangeCopy(BfMatBlockCoo const *matBlockCoo, BfSize i0, BfSize i1) {
  BF_ERROR_BEGIN();

  if (i0 > i1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMat const *mat = bfMatBlockCooConstToMatConst(matBlockCoo);
  BfMatBlock const *matBlock = bfMatBlockCooConstToMatBlockConst(matBlockCoo);

  BfMatBlockCoo *blockRow = NULL;

  BfPtrArray indexedRowBlocks;
  bfInitPtrArrayWithDefaultCapacity(&indexedRowBlocks);
  HANDLE_ERROR();

  BfSize numBlocks = bfMatBlockNumBlocks(matBlock);

  for (BfSize k = 0; k < numBlocks; ++k) {
    BfMat const *block = BLOCK(matBlockCoo, k);

    BfSize m_ = bfMatGetNumRows(block);

    /* The range [i0_, i1_) is the row index range of `block` with
     * respect to `matBlockCoo`'s index system. */
    BfSize i0_ = BLOCK_ROW_OFFSET(matBlockCoo, k);
    BfSize i1_ = i0_ + m_;

    /* Skip blocks which don't overlap with the row range [i0, i1). */
    if (i1_ <= i0 || i1 <= i0_)
      continue;

    /* The starting column index of `block` with respect to
     * `matBlockCoo`'s index system. */
    BfSize j0_ = BLOCK_COL_OFFSET(matBlockCoo, k);

    /* The range of rows which should be copied from `block` with
     * respect to `block`'s index system. */
    BfSize i0__ = i0_ < i0 ? i0 - i0_ : 0;
    BfSize i1__ = m_ - (i1 < i1_ ? i1_ - i1 : 0);

    /* Copy the row range---the whole point here is that this can
     * easily be a subset of the rows of `block`! */
    BfMat *blockRowRange = bfMatGetRowRangeCopy(block, i0__, i1__);
    HANDLE_ERROR();

    // TODO: explain how the indexing for i0 works
    BfIndexedMat *indexedRowBlock = bfIndexedMatNewFromMat(
      i0_ < i0 ? 0 : i0_ - i0,
      j0_,
      blockRowRange,
      BF_POLICY_STEAL);

    bfPtrArrayAppend(&indexedRowBlocks, indexedRowBlock);
    HANDLE_ERROR();

    bfMatDelete(&blockRowRange);
  }

  BfSize m = i1 - i0;
  BfSize n = bfMatGetNumCols(mat);

  blockRow = bfMatBlockCooNewFromIndexedBlocks(m, n, &indexedRowBlocks, BF_POLICY_STEAL);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  for (BfSize k = 0; k < bfPtrArraySize(&indexedRowBlocks); ++k) {
    BfIndexedMat *indexedRowBlock = bfPtrArrayGet(&indexedRowBlocks, k);
    bfIndexedMatDelete(&indexedRowBlock);
  }
  bfPtrArrayDeinit(&indexedRowBlocks);

  return bfMatBlockCooToMat(blockRow);
}

BfMat *bfMatBlockCooMul(BfMat const *op1, BfMat const *op2) {
  BF_ERROR_BEGIN();

  BfMatBlock const *matBlock1 = bfMatConstToMatBlockConst(op1);
  BfMatBlockCoo const *matBlockCoo1 = bfMatConstToMatBlockCooConst(op1);

  BfSize numRows = bfMatGetNumRows(op1);
  BfSize numCols = bfMatGetNumCols(op2);

  if (bfMatGetNumCols(op1) != bfMatGetNumRows(op2))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numBlocks = bfMatBlockCooNumBlocks(matBlockCoo1);

  BfMat *result = NULL;
  BfMat *op2Rows = NULL;
  BfMat *tmp = NULL;
  BfMat *resultRows = NULL;

  result = bfMatZerosLike(op2, numRows, numCols);
  HANDLE_ERROR();

  for (BfSize k = 0, i0, i1, j0, j1; k < numBlocks; ++k) {
    i0 = matBlock1->rowOffset[matBlockCoo1->rowInd[k]];
    i1 = matBlock1->rowOffset[matBlockCoo1->rowInd[k] + 1];
    j0 = matBlock1->colOffset[matBlockCoo1->colInd[k]];
    j1 = matBlock1->colOffset[matBlockCoo1->colInd[k] + 1];

    op2Rows = bfMatGetRowRange((BfMat *)op2, j0, j1);
    tmp = bfMatMul(matBlock1->block[k], op2Rows);
    resultRows = bfMatGetRowRange(result, i0, i1);
    bfMatAddInplace(resultRows, tmp);

    bfMatDelete(&resultRows);
    bfMatDelete(&tmp);
    bfMatDelete(&op2Rows);
  }

  BF_ERROR_END()
    bfMatDelete(&result);

  return result;

}

BfVec *bfMatBlockCooMulVec(BfMatBlockCoo const *matBlockCoo, BfVec const *vec) {
  BF_ERROR_BEGIN();

  BfMat const *mat = bfMatBlockCooConstToMatConst(matBlockCoo);

  if (bfMatGetNumCols(mat) != vec->size)
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfSize resultSize = bfMatGetNumRows(mat);

  BfVec *result = NULL;
  switch (bfVecGetType(vec)) {
  case BF_TYPE_VEC_REAL:
    result = bfVecRealToVec(bfVecRealNewWithValue(resultSize, 0));
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  BfSize numBlocks = bfMatBlockCooNumBlocks(matBlockCoo);
  for (BfSize k = 0; k < numBlocks; ++k) {
    BfMat const *block = BLOCK(matBlockCoo, k);

    BfSize i0 = BLOCK_ROW_OFFSET(matBlockCoo, k);
    BfSize i1 = i0 + bfMatGetNumRows(block);

    BfSize j0 = BLOCK_COL_OFFSET(matBlockCoo, k);
    BfSize j1 = j0 + bfMatGetNumCols(block);

    /* Multiply the current block with the relevant part of `vec`: */
    BfVec const *subvecView = bfVecGetSubvecViewConst(vec, j0, j1);
    BfVec *tmp = bfMatMulVec(block, subvecView);
    bfVecDelete((BfVec **)&subvecView);

    /* Assign the result to the relevant part of `result`: */
    BfVec *resultSubvecView = bfVecGetSubvecView(result, i0, i1);
    bfVecAddInplace(resultSubvecView, tmp);
    bfVecDelete(&resultSubvecView);

    bfVecDelete(&tmp);
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return result;
}

BfVec *bfMatBlockCooRmulVec(BfMatBlockCoo const *matBlockCoo, BfVec const *vec) {
  BF_ERROR_BEGIN();

  BfMat const *mat = bfMatBlockCooConstToMatConst(matBlockCoo);

  if (bfMatGetNumRows(mat) != vec->size)
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfSize resultSize = bfMatGetNumCols(mat);
  BfVec *result = NULL;
  switch (bfVecGetType(vec)) {
  case BF_TYPE_VEC_REAL:
    result = bfVecRealToVec(bfVecRealNewWithValue(resultSize, 0));
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  BfSize numBlocks = bfMatBlockCooNumBlocks(matBlockCoo);
  for (BfSize k = 0; k < numBlocks; ++k) {
    BfMat const *block = BLOCK(matBlockCoo, k);

    BfSize i0 = BLOCK_ROW_OFFSET(matBlockCoo, k);
    BfSize i1 = i0 + bfMatGetNumRows(block);

    BfSize j0 = BLOCK_COL_OFFSET(matBlockCoo, k);
    BfSize j1 = j0 + bfMatGetNumCols(block);

    /* Multiply the current block with the relevant part of `vec`: */
    BfVec const *subvecView = bfVecGetSubvecViewConst(vec, i0, i1);
    BfVec *tmp = bfMatRmulVec(block, subvecView);
    bfVecDelete((BfVec **)&subvecView);

    /* Assign the result to the relevant part of `result`: */
    BfVec *resultSubvecView = bfVecGetSubvecView(result, j0, j1);
    bfVecAddInplace(resultSubvecView, tmp);
    bfVecDelete(&resultSubvecView);
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return result;
}

void bfMatBlockCooNegate(BfMat *mat) {
  BF_ERROR_BEGIN();

  BfMatBlock *matBlock = NULL;

  matBlock = bfMatToMatBlock(mat);
  HANDLE_ERROR();

  BfSize numBlocks = bfMatBlockNumBlocks(matBlock);

  for (BfSize i = 0; i < numBlocks; ++i)
    bfMatNegate(matBlock->block[i]);

  BF_ERROR_END() {}
}

static void setColumnNonzerosForBlock(BfMat const *mat, bool *nonzero) {
  BfSize m = bfMatGetNumCols(mat);
  BfSize n = bfMatGetNumCols(mat);
  BfType type = bfMatGetType(mat);
  switch (type) {
  case BF_TYPE_MAT_BLOCK_COO: {
    BfMatBlockCoo const *matBlockCoo = bfMatConstToMatBlockCooConst(mat);
    BfSize numBlocks = bfMatBlockCooNumBlocks(matBlockCoo);
    for (BfSize k = 0; k < numBlocks; ++k) {
      BfMat const *block = BLOCK(matBlockCoo, k);
      BfSize j0 = BLOCK_COL_OFFSET(matBlockCoo, k);
      setColumnNonzerosForBlock(block, &nonzero[j0]);
    }
    break;
  }
  case BF_TYPE_MAT_IDENTITY: {
    BF_ASSERT(m == n);
    for (BfSize j = 0; j < n; ++j) nonzero[j] = true;
    break;
  }
  case BF_TYPE_MAT_DENSE_REAL: {
    for (BfSize j = 0; j < n; ++j) nonzero[j] = true;
    break;
  }
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
  }
}


BfSizeArray *bfMatBlockCooGetNonzeroColumnRanges(BfMatBlockCoo const *matBlockCoo) {
  BF_ERROR_BEGIN();

  BfSizeArray *nonzeroColumnRanges = bfSizeArrayNew();
  HANDLE_ERROR();

  bfSizeArrayInitWithDefaultCapacity(nonzeroColumnRanges);
  HANDLE_ERROR();

  BfSize n = bfMatGetNumCols(bfMatBlockCooConstToMatConst(matBlockCoo));

  bool *nonzero = bfMemAllocAndZero(n, sizeof(bool));
  if (nonzero == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfSize numBlocks = bfMatBlockCooNumBlocks(matBlockCoo);
  for (BfSize k = 0; k < numBlocks; ++k) {
    BfMat const *block = BLOCK(matBlockCoo, k);
    BfSize j0 = BLOCK_COL_OFFSET(matBlockCoo, k);
    setColumnNonzerosForBlock(block, &nonzero[j0]);
  }

  BfSize j0Next = 0;

  while (j0Next < n) {
    BfSize j0 = j0Next;
    while (j0 < n && !nonzero[j0]) ++j0;

    BfSize j1 = j0;
    while (j1 < n && nonzero[j1]) ++j1;

    if (j0 < j1) {
      bfSizeArrayAppend(nonzeroColumnRanges, j0);
      bfSizeArrayAppend(nonzeroColumnRanges, j1);
    }

    j0Next = j1;
  }

  BF_ASSERT(j0Next == n);

  BF_ERROR_END() {
    BF_DIE();
  }

  return nonzeroColumnRanges;
}

void bfMatBlockCooPrintBlocksDeep(BfMatBlockCoo const *matBlockCoo, FILE *fp, BfSize i0, BfSize j0, BfSize depth) {
  BfMat const *mat = bfMatBlockCooConstToMatConst(matBlockCoo);

  BfSize i1 = i0 + bfMatGetNumRows(mat);
  BfSize j1 = j0 + bfMatGetNumCols(mat);

  fprintf(fp, "%u %lu %lu %lu %lu %lu\n", BF_TYPE_MAT_BLOCK_COO, i0, i1, j0, j1, depth);

  BfSize numBlocks = bfMatBlockCooNumBlocks(matBlockCoo);
  for (BfSize k = 0; k < numBlocks; ++k) {
    BfSize i0_ = i0 + BLOCK_ROW_OFFSET(matBlockCoo, k);
    BfSize j0_ = j0 + BLOCK_COL_OFFSET(matBlockCoo, k);
    BfMat const *block = BLOCK(matBlockCoo, k);
    bfMatPrintBlocksDeep(block, fp, i0_, j0_, depth + 1);
  }
}

/** Interface: MatBlock */

static BfMatBlockVtable MAT_BLOCK_VTABLE = {
  .NumBlocks = (__typeof__(&bfMatBlockNumBlocks))bfMatBlockCooNumBlocks,
  .GetNumRowBlocks = (__typeof__(&bfMatBlockGetNumRowBlocks))bfMatBlockCooGetNumRowBlocks,
  .GetNumColBlocks = (__typeof__(&bfMatBlockGetNumRowBlocks))bfMatBlockCooGetNumColBlocks,
  .GetRowOffset = (__typeof__(&bfMatBlockGetRowOffset))bfMatBlockCooGetRowOffset,
  .GetColOffset = (__typeof__(&bfMatBlockGetColOffset))bfMatBlockCooGetColOffset,
  .GetBlockConst = (__typeof__(&bfMatBlockGetBlockConst))bfMatBlockCooGetBlockConst,
};

BfSize bfMatBlockCooNumBlocks(BfMatBlockCoo const *matBlockCoo) {
  return matBlockCoo->numBlocks;
}

BfSize bfMatBlockCooGetNumRowBlocks(BfMatBlockCoo const *matBlockCoo) {
  return NUM_ROW_BLOCKS(matBlockCoo);
}

BfSize bfMatBlockCooGetNumColBlocks(BfMatBlockCoo const *matBlockCoo) {
  return NUM_COL_BLOCKS(matBlockCoo);
}

BfSize bfMatBlockCooGetRowOffset(BfMatBlockCoo const *matBlockCoo, BfSize i) {
  return ROW_OFFSET(matBlockCoo, i);
}

BfSize bfMatBlockCooGetColOffset(BfMatBlockCoo const *matBlockCoo, BfSize j) {
  return COL_OFFSET(matBlockCoo, j);
}

BfMat const *bfMatBlockCooGetBlockConst(BfMatBlockCoo const *matBlockCoo, BfSize i, BfSize j) {
  BF_ERROR_BEGIN();

  BfMatBlock const *matBlock = bfMatBlockCooConstToMatBlockConst(matBlockCoo);

  BfSize numBlocks = bfMatBlockNumBlocks(matBlock);

  BfMat const *block = NULL;

  BfSize k = 0;
  for (; k < numBlocks; ++k) {
    if (ROW_IND(matBlockCoo, k) == i && COL_IND(matBlockCoo, k) == j) {
      block = bfMatGet(BLOCK(matBlockCoo, k), BF_POLICY_VIEW);
      break;
    }
  }

  BfSize i0 = ROW_OFFSET(matBlockCoo, i);
  BfSize m = ROW_OFFSET(matBlockCoo, i + 1) - i0;

  BfSize j0 = COL_OFFSET(matBlockCoo, j);
  BfSize n = COL_OFFSET(matBlockCoo, j + 1) - j0;

  BfMatZero *zeroBlock = NULL;

  if (block == NULL) {
    zeroBlock = bfMatZeroNew();
    HANDLE_ERROR();

    bfMatZeroInit(zeroBlock, m, n);
    HANDLE_ERROR();

    block = bfMatZeroConstToMatConst((BfMatZero const *)zeroBlock);
  } else {
    if (bfMatGetNumRows(block) != m || bfMatGetNumCols(block) != n)
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return block;
}

/** Upcasting: MatBlockCoo -> Mat */

BfMat *bfMatBlockCooToMat(BfMatBlockCoo *matBlockCoo) {
  return &matBlockCoo->super.super;
}

BfMat const *bfMatBlockCooConstToMatConst(BfMatBlockCoo const *matBlockCoo) {
  return &matBlockCoo->super.super;
}

/** Upcasting: MatBlockCoo -> MatBlock */

BfMatBlock *bfMatBlockCooToMatBlock(BfMatBlockCoo *matBlockCoo) {
  return &matBlockCoo->super;
}

BfMatBlock const *bfMatBlockCooConstToMatBlockConst(BfMatBlockCoo const *matBlockCoo) {
  return &matBlockCoo->super;
}

/** Downcasting: Mat -> MatBlockCoo */

BfMatBlockCoo *bfMatToMatBlockCoo(BfMat *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_BLOCK_COO)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatBlockCoo *)mat;
  }
}

BfMatBlockCoo const *bfMatConstToMatBlockCooConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_BLOCK_COO)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatBlockCoo const *)mat;
  }
}

/** Downcasting: MatBlockCoo -> MatBlock */

BfMatBlockCoo const *bfMatBlockConstToMatBlockCooConst(BfMatBlock const *matBlock) {
  return bfMatConstToMatBlockCooConst(&matBlock->super);
}

/** Implementation: MatBlockCoo */

void bfMatBlockCooInvalidate(BfMatBlockCoo *matBlockCoo) {
  bfMatBlockInvalidate(&matBlockCoo->super);

  matBlockCoo->numBlocks = BF_SIZE_BAD_VALUE;
  matBlockCoo->rowInd = NULL;
  matBlockCoo->colInd = NULL;
}

BfMatBlockCoo *bfMatBlockCooNew() {
  BF_ERROR_BEGIN();

  BfMatBlockCoo *matBlockCoo = bfMemAlloc(1, sizeof(BfMatBlockCoo));
  if (matBlockCoo == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  bfMatBlockCooInvalidate(matBlockCoo);

  BF_ERROR_END() {
    matBlockCoo = NULL;
  }

  return matBlockCoo;
}

BfMatBlockCoo *bfMatBlockCooNewFromArrays(BfSizeArray const *rowOffsets,
                                          BfSizeArray const *colOffsets,
                                          BfSizeArray const *rowInds,
                                          BfSizeArray const *colInds,
                                          BfPtrArray const *blocks) {
  BF_ERROR_BEGIN();

  BfMatBlockCoo *matBlockCoo = bfMatBlockCooNew();
  HANDLE_ERROR();

  BfMatBlock *matBlock = bfMatBlockCooToMatBlock(matBlockCoo);

  BfSize numBlockRows = bfSizeArrayGetSize(rowOffsets) - 1;
  BfSize numBlockCols = bfSizeArrayGetSize(colOffsets) - 1;
  BfSize numBlocks = bfPtrArraySize(blocks);

  bfMatBlockCooInit(matBlockCoo, numBlockRows, numBlockCols, numBlocks);
  HANDLE_ERROR();

  bfPtrArrayCopyData(blocks, (BfPtr *)matBlock->block);
  bfSizeArrayCopyData(rowOffsets, matBlock->rowOffset);
  bfSizeArrayCopyData(colOffsets, matBlock->colOffset);
  bfSizeArrayCopyData(rowInds, matBlockCoo->rowInd);
  bfSizeArrayCopyData(colInds, matBlockCoo->colInd);

  BF_ERROR_END() {
    bfMatBlockCooDeinitAndDealloc(&matBlockCoo);
  }

  return matBlockCoo;
}

BfMatBlockCoo *bfMatBlockCooNewColFromBlocks(BfPtrArray *blocks, BfPolicy policy) {
  BF_ERROR_BEGIN();

  if (bfPtrArrayIsEmpty(blocks))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatBlockCoo *matBlockCoo = bfMatBlockCooNew();
  HANDLE_ERROR();

  BfSize numBlocks = bfPtrArraySize(blocks);

  bfMatBlockInit(&matBlockCoo->super, &MAT_VTABLE, &MAT_BLOCK_VTABLE,
                 numBlocks, numBlocks, 1);
  HANDLE_ERROR();

  matBlockCoo->numBlocks = numBlocks;

  BfMat *block = bfPtrArrayGet(blocks, 0);

  /* Compute the row and column offsets */

  COL_OFFSET(matBlockCoo, 0) = 0;
  COL_OFFSET(matBlockCoo, 1) = bfMatGetNumCols(block);

  ROW_OFFSET(matBlockCoo, 0) = 0;
  for (BfSize k = 0; k < numBlocks; ++k) {
    block = bfPtrArrayGet(blocks, k);
    BLOCK(matBlockCoo, k) = bfMatGet(block, policy);
    ROW_OFFSET(matBlockCoo, k + 1) = bfMatGetNumRows(block);
  }
  bfSizeRunningSum(numBlocks + 1, &ROW_OFFSET(matBlockCoo, 0));

  /* Set up the row and column indices */

  matBlockCoo->rowInd = bfMemAlloc(numBlocks, sizeof(BfSize));
  if (matBlockCoo->rowInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  matBlockCoo->colInd = bfMemAlloc(numBlocks, sizeof(BfSize));
  if (matBlockCoo->colInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize k = 0; k < numBlocks; ++k) {
    ROW_IND(matBlockCoo, k) = k;
    COL_IND(matBlockCoo, k) = 0;
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return matBlockCoo;
}

BfMatBlockCoo *bfMatBlockCooNewRowFromBlocks(BfPtrArray *blocks, BfPolicy policy) {
  BF_ERROR_BEGIN();

  if (bfPtrArrayIsEmpty(blocks))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatBlockCoo *matBlockCoo = bfMatBlockCooNew();
  HANDLE_ERROR();

  BfSize numBlocks = bfPtrArraySize(blocks);

  bfMatBlockInit(&matBlockCoo->super, &MAT_VTABLE, &MAT_BLOCK_VTABLE,
                 numBlocks, 1, numBlocks);
  HANDLE_ERROR();

  matBlockCoo->numBlocks = numBlocks;

  BfMat *block = bfPtrArrayGet(blocks, 0);

  /* Compute the row and column offsets */

  ROW_OFFSET(matBlockCoo, 0) = 0;
  ROW_OFFSET(matBlockCoo, 1) = bfMatGetNumRows(block);

  COL_OFFSET(matBlockCoo, 0) = 0;
  for (BfSize k = 0; k < numBlocks; ++k) {
    block = bfPtrArrayGet(blocks, k);
    BLOCK(matBlockCoo, k) = bfMatGet(block, policy);
    COL_OFFSET(matBlockCoo, k + 1) = bfMatGetNumCols(block);
  }
  bfSizeRunningSum(numBlocks + 1, &COL_OFFSET(matBlockCoo, 0));

  /* Set up the row and column indices */

  matBlockCoo->rowInd = bfMemAlloc(numBlocks, sizeof(BfSize));
  if (matBlockCoo->rowInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  matBlockCoo->colInd = bfMemAlloc(numBlocks, sizeof(BfSize));
  if (matBlockCoo->colInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize k = 0; k < numBlocks; ++k) {
    ROW_IND(matBlockCoo, k) = 0;
    COL_IND(matBlockCoo, k) = k;
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return matBlockCoo;
}

BfMatBlockCoo *bfMatBlockCooNewFromIndexedBlocks(BfSize numRows, BfSize numCols, BfPtrArray *indexedBlocks, BfPolicy policy) {
  BF_ERROR_BEGIN();

  BfMatBlockCoo *matBlockCoo = NULL;

  BfSize numBlocks = bfPtrArraySize(indexedBlocks);

  /* Extract the blocks from the array of indexed blocks */

  BfPtrArray blocks;
  bfInitPtrArrayWithDefaultCapacity(&blocks);
  HANDLE_ERROR();

  for (BfSize k = 0; k < numBlocks; ++k) {
    BfIndexedMat *indexedBlock = bfPtrArrayGet(indexedBlocks, k);

    BfMat *block = bfMatGet(indexedBlock->mat, policy);
    HANDLE_ERROR();

    bfPtrArrayAppend(&blocks, block);
    HANDLE_ERROR();
  }

  /* Build and accumulate arrays containing the row and column offsets
   * determined by `indexedBlocks`. */

  BfSizeArray rowOffsets;
  bfSizeArrayInitWithDefaultCapacity(&rowOffsets);
  HANDLE_ERROR();

  BfSizeArray colOffsets;
  bfSizeArrayInitWithDefaultCapacity(&colOffsets);
  HANDLE_ERROR();

  bfSizeArrayAppend(&rowOffsets, 0);
  bfSizeArrayAppend(&rowOffsets, numRows);

  bfSizeArrayAppend(&colOffsets, 0);
  bfSizeArrayAppend(&colOffsets, numCols);

  for (BfSize k = 0; k < numBlocks; ++k) {
    BfIndexedMat const *indexedBlock = bfPtrArrayGet(indexedBlocks, k);
    BfMat const *block = indexedBlock->mat;

    BfSize i0 = indexedBlock->i0;
    BfSize i1 = i0 + bfMatGetNumRows(block);

    BfSize j0 = indexedBlock->j0;
    BfSize j1 = j0 + bfMatGetNumCols(block);

    if (!bfSizeArrayContains(&rowOffsets, i0)) {
      bfSizeArrayInsertSorted(&rowOffsets, i0);
      HANDLE_ERROR();
    }

    if (!bfSizeArrayContains(&rowOffsets, i1)) {
      bfSizeArrayInsertSorted(&rowOffsets, i1);
      HANDLE_ERROR();
    }

    if (!bfSizeArrayContains(&colOffsets, j0)) {
      bfSizeArrayInsertSorted(&colOffsets, j0);
      HANDLE_ERROR();
    }

    if (!bfSizeArrayContains(&colOffsets, j1)) {
      bfSizeArrayInsertSorted(&colOffsets, j1);
      HANDLE_ERROR();
    }
  }

  /* Find rowInd and colInd for each block. */

  BfSizeArray rowInds;
  bfSizeArrayInitWithDefaultCapacity(&rowInds);
  HANDLE_ERROR();

  BfSizeArray colInds;
  bfSizeArrayInitWithDefaultCapacity(&colInds);
  HANDLE_ERROR();

  for (BfSize k = 0; k < numBlocks; ++k) {
    BfIndexedMat const *indexedBlock = bfPtrArrayGet(indexedBlocks, k);
    BfMat const *block = indexedBlock->mat;

    BfSize i0 = indexedBlock->i0;
    BfSize i1 = i0 + bfMatGetNumRows(block);

    BfSize j0 = indexedBlock->j0;
    BfSize j1 = j0 + bfMatGetNumCols(block);

    BfSize rowInd = bfSizeArrayFindFirst(&rowOffsets, i0);
    if (rowInd == BF_SIZE_BAD_VALUE)
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

    BF_ASSERT(bfSizeArrayGet(&rowOffsets, rowInd + 1) == i1);

    bfSizeArrayAppend(&rowInds, rowInd);
    HANDLE_ERROR();

    BfSize colInd = bfSizeArrayFindFirst(&colOffsets, j0);
    if (colInd == BF_SIZE_BAD_VALUE)
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

    BF_ASSERT(bfSizeArrayGet(&colOffsets, colInd + 1) == j1);

    bfSizeArrayAppend(&colInds, colInd);
    HANDLE_ERROR();
  }

  /* Create and initialize MatBlockCoo instance using inds and
   * offsets: */

  matBlockCoo = bfMatBlockCooNewFromArrays(
    &rowOffsets, &colOffsets, &rowInds, &colInds, &blocks);
  HANDLE_ERROR();

  BF_ERROR_END() {
    bfMatBlockCooDeinitAndDealloc(&matBlockCoo);
  }

  bfPtrArrayDeinit(&blocks);

  bfSizeArrayDeinit(&rowOffsets);
  bfSizeArrayDeinit(&colOffsets);

  bfSizeArrayDeinit(&rowInds);
  bfSizeArrayDeinit(&colInds);

  return matBlockCoo;
}

void bfMatBlockCooInit(BfMatBlockCoo *mat, BfSize numBlockRows,
                       BfSize numBlockCols, BfSize numBlocks)
{
  BF_ERROR_BEGIN();

  bfMatBlockInit(&mat->super,
                 &MAT_VTABLE, &MAT_BLOCK_VTABLE,
                 numBlocks, numBlockRows, numBlockCols);
  HANDLE_ERROR();

  mat->numBlocks = numBlocks;

  /* allocate and initialize row indices to dummy values */
  mat->rowInd = bfMemAlloc(numBlocks, sizeof(BfSize));
  if (mat->rowInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  for (BfSize i = 0; i < numBlocks; ++i)
    mat->rowInd[i] = BF_SIZE_BAD_VALUE;

  /* allocate and initialize column indices to dummy values */
  mat->colInd = bfMemAlloc(numBlocks, sizeof(BfSize));
  if (mat->colInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  for (BfSize i = 0; i < numBlocks; ++i)
    mat->colInd[i] = BF_SIZE_BAD_VALUE;

  BF_ERROR_END() {
    bfMatBlockDeinit(&mat->super);
    bfMemFree(mat->rowInd);
    bfMemFree(mat->colInd);
  }
}

void bfMatBlockCooDeinit(BfMatBlockCoo *matBlockCoo) {
  BfMat *mat = bfMatBlockCooToMat(matBlockCoo);
  if (!bfMatIsView(mat)) {
    bfMemFree(matBlockCoo->rowInd);
    bfMemFree(matBlockCoo->colInd);
  }

  matBlockCoo->rowInd = NULL;
  matBlockCoo->colInd = NULL;

  bfMatBlockDeinit(&matBlockCoo->super);

  matBlockCoo->numBlocks = BF_SIZE_BAD_VALUE;
}

void bfMatBlockCooDealloc(BfMatBlockCoo **mat) {
  bfMemFree(*mat);
  *mat = NULL;
}

void bfMatBlockCooDeinitAndDealloc(BfMatBlockCoo **mat) {
  bfMatBlockCooDeinit(*mat);
  bfMatBlockCooDealloc(mat);
}
