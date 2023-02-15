#include <bf/mat_block_coo.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/ptr_array.h>

/** Helper macros: */

#define ROW_IND(mat, k) mat->rowInd[k]
#define COL_IND(mat, k) mat->colInd[k]
#define ROW_OFFSET(mat, i) mat->super.rowOffset[i]
#define COL_OFFSET(mat, j) mat->super.colOffset[j]
#define BLOCK_ROW_OFFSET(mat, k) ROW_OFFSET(mat, mat->rowInd[k])
#define BLOCK_COL_OFFSET(mat, k) COL_OFFSET(mat, mat->colInd[k])
#define BLOCK(mat, k) mat->super.block[k]

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .Copy = (__typeof__(&bfMatBlockCooCopy))bfMatBlockCooCopy,
  .GetRowCopy = (__typeof__(&bfMatBlockCooGetRowCopy))bfMatBlockCooGetRowCopy,
  .Delete = (__typeof__(&bfMatBlockCooDelete))bfMatBlockCooDelete,
  .GetType = (__typeof__(&bfMatBlockCooGetType))bfMatBlockCooGetType,
  .NumBytes = (__typeof__(&bfMatBlockCooNumBytes))bfMatBlockCooNumBytes,
  .GetNumRows = (__typeof__(&bfMatBlockCooGetNumRows))bfMatBlockCooGetNumRows,
  .GetNumCols = (__typeof__(&bfMatBlockCooGetNumCols))bfMatBlockCooGetNumCols,
  .GetRowRangeCopy = (__typeof__(&bfMatGetRowRangeCopy))bfMatBlockCooGetRowRangeCopy,
  .Mul = (__typeof__(&bfMatBlockCooMul))bfMatBlockCooMul,
  .Negate = (__typeof__(&bfMatBlockCooNegate))bfMatBlockCooNegate,
};

BfMat *bfMatBlockCooCopy(BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

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

  END_ERROR_HANDLING() {
    bfMatDelete(&blockCopy);
    bfMatBlockCooDeinitAndDealloc(&matBlockCooCopy);
  }

  return bfMatBlockCooToMat(matBlockCooCopy);
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
  BEGIN_ERROR_HANDLING();

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
    BfMatBlockCooEntry *entry = malloc(sizeof(BfMatBlockCooEntry));
    entry->i0 = BLOCK_ROW_OFFSET(matBlockCoo, k);
    entry->j0 = BLOCK_COL_OFFSET(matBlockCoo, k);
    entry->block = BLOCK(matBlockCoo, k);
    bfPtrArrayAppend(&blocks, entry);
  }

  /* Sort the array of blocks into dictionary order */
  bfPtrArraySort(&blocks, (BfPtrCmp)cmpColumnOrder);

  /* Verify that none of the collected blocks are overlapping */
#if BF_DEBUG
  BfSize j1 = 0;
  for (BfSize k = 0; k < bfPtrArraySize(&blocks); ++k) {
    BfMatBlockCooEntry const *entry = bfPtrArrayGet(&blocks, k);
    assert(j1 <= entry->j0);
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

  assert(rowCopy->size == bfMatGetNumCols(mat));

  END_ERROR_HANDLING() {}

  for (BfSize k = 0; k < bfPtrArraySize(&blocks); ++k)
    free(bfPtrArrayGet(&blocks, k));
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
  BEGIN_ERROR_HANDLING();

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

    /* The starting column index of `block` with respect to
     * `matBlockCoo`'s index system. */
    BfSize j0_ = BLOCK_COL_OFFSET(matBlockCoo, k);

    /* Skip blocks which don't overlap with the row range [i0, i1). */
    assert(i0 <= i1 && i0_ <= i1_);
    if (i1_ <= i0 || i1 <= i0_)
      continue;

    /* The range of rows which should be copied from `block` with
     * respect to `block`'s index system. */
    BfSize i0__ = i0_ < i0 ? i0 - i0_ : 0;
    BfSize i1__ = m_ - (i1 < i1_ ? i1_ - i1 : 0);

    /* Copy the row range---the whole point here is that this can
     * easily be a subset of the rows of `block`! */
    BfMat *blockRowRange = bfMatGetRowRangeCopy(block, i0__, i1__);
    HANDLE_ERROR();

    BfIndexedMat *indexedRowBlock = malloc(sizeof(BfIndexedMat));
    if (indexedRowBlock == NULL)
      RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

    indexedRowBlock->i0 = i0_ < i0 ? 0 : i0_ - i0;
    indexedRowBlock->j0 = j0_;
    indexedRowBlock->mat = blockRowRange;

    assert(indexedRowBlock->i0 <= m_);
    assert(indexedRowBlock->mat != NULL);

    bfPtrArrayAppend(&indexedRowBlocks, indexedRowBlock);
    HANDLE_ERROR();
  }

  BfSize m = i1 - i0;
  BfSize n = bfMatGetNumCols(mat);

  blockRow = bfMatBlockCooNewFromIndexedBlocks(m, n, &indexedRowBlocks);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    bfMatBlockCooDeinitAndDealloc(&blockRow);
  }

  /* Free the BfIndexedMat wrappers (but not the wrapped mats!) */
  for (BfSize k = 0; k < bfPtrArraySize(&indexedRowBlocks); ++k)
    free(bfPtrArrayGet(&indexedRowBlocks, k));

  bfPtrArrayDeinit(&indexedRowBlocks);

  return bfMatBlockCooToMat(blockRow);
}

BfMat *bfMatBlockCooMul(BfMat const *op1, BfMat const *op2) {
  BEGIN_ERROR_HANDLING();

  BfMatBlock const *matBlock1 = bfMatConstToMatBlockConst(op1);
  BfMatBlockCoo const *matBlockCoo1 = bfMatConstToMatBlockCooConst(op1);

  BfSize numRows = bfMatGetNumRows(op1);
  BfSize numCols = bfMatGetNumCols(op2);
  BfSize numBlocks = bfMatBlockCooNumBlocks(matBlock1);

  BfMat *result = NULL;
  BfMat *block = NULL;
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

    block = matBlock1->block[k];
    op2Rows = bfMatGetRowRange((BfMat *)op2, j0, j1);
    tmp = bfMatMul(block, op2Rows);
    resultRows = bfMatGetRowRange(result, i0, i1);
    bfMatAddInplace(resultRows, tmp);

    bfMatDelete(&tmp);
  }

  END_ERROR_HANDLING()
    bfMatDelete(&result);

  return result;

}

void bfMatBlockCooNegate(BfMat *mat) {
  BEGIN_ERROR_HANDLING();

  BfMatBlock *matBlock = NULL;

  matBlock = bfMatToMatBlock(mat);
  HANDLE_ERROR();

  BfSize numBlocks = bfMatBlockNumBlocks(matBlock);

  for (BfSize i = 0; i < numBlocks; ++i)
    bfMatNegate(matBlock->block[i]);

  END_ERROR_HANDLING() {}
}

/** Interface: MatBlock */

static BfMatBlockVtable MAT_BLOCK_VTABLE = {
  .NumBlocks = (__typeof__(&bfMatBlockCooNumBlocks))bfMatBlockCooNumBlocks
};

BfSize bfMatBlockCooNumBlocks(BfMatBlock const *mat) {
  return bfMatBlockConstToMatBlockCooConst(mat)->numBlocks;
}

/** Upcasting: MatBlockCoo -> Mat */

BfMat *bfMatBlockCooToMat(BfMatBlockCoo *matBlockCoo) {
  return &matBlockCoo->super.super;
}

BfMat const *bfMatBlockCooConstToMatConst(BfMatBlockCoo const *matBlockCoo) {
  return &matBlockCoo->super.super;
}

/** Upcasting: MatBlockCoo -> MatBlock */

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
  return malloc(sizeof(BfMatBlockCoo));
}

void bfMatBlockCooInit(BfMatBlockCoo *mat, BfSize numBlockRows,
                       BfSize numBlockCols, BfSize numBlocks)
{
  BEGIN_ERROR_HANDLING();

  bfMatBlockInit(&mat->super,
                 &MAT_VTABLE, &MAT_BLOCK_VTABLE,
                 numBlocks, numBlockRows, numBlockCols);
  HANDLE_ERROR();

  mat->numBlocks = numBlocks;

  /* allocate and initialize row indices to dummy values */
  mat->rowInd = malloc(numBlocks*sizeof(BfSize));
  if (mat->rowInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  for (BfSize i = 0; i < numBlocks; ++i)
    mat->rowInd[i] = BF_SIZE_BAD_VALUE;

  /* allocate and initialize column indices to dummy values */
  mat->colInd = malloc(numBlocks*sizeof(BfSize));
  if (mat->colInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  for (BfSize i = 0; i < numBlocks; ++i)
    mat->colInd[i] = BF_SIZE_BAD_VALUE;

  END_ERROR_HANDLING() {
    bfMatBlockDeinit(&mat->super);
    free(mat->rowInd);
    free(mat->colInd);
  }
}

void bfMatBlockCooDeinit(BfMatBlockCoo *mat) {
  (void)mat;
  assert(false);
}

void bfMatBlockCooDealloc(BfMatBlockCoo **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatBlockCooDeinitAndDealloc(BfMatBlockCoo **mat) {
  bfMatBlockCooDeinit(*mat);
  bfMatBlockCooDealloc(mat);
}
