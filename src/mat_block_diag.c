#include <bf/mat_block_diag.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_zero.h>
#include <bf/ptr_array.h>
#include <bf/util.h>
#include <bf/vec_zero.h>

/** Helper macros: */

#define NUM_ROW_BLOCKS(mat) mat->super.super.numRows
#define NUM_COL_BLOCKS(mat) mat->super.super.numCols
#define NUM_BLOCKS(mat) NUM_ROW_BLOCKS(mat)*NUM_COL_BLOCKS(mat)
#define ROW_OFFSET(mat, k) mat->super.rowOffset[k]
#define COL_OFFSET(mat, k) mat->super.colOffset[k]
#define NUM_BLOCK_ROWS(mat, k) ROW_OFFSET(mat, k + 1) - ROW_OFFSET(mat, k)
#define NUM_BLOCK_COLS(mat, k) COL_OFFSET(mat, k + 1) - COL_OFFSET(mat, k)
#define BLOCK(mat, k) mat->super.block[k]

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .Copy = (__typeof__(&bfMatBlockDiagCopy))bfMatBlockDiagCopy,
  .GetRowCopy = (__typeof__(&bfMatBlockDiagGetRowCopy))bfMatBlockDiagGetRowCopy,
  .Delete = (__typeof__(&bfMatBlockDiagDelete))bfMatBlockDiagDelete,
  .GetType = (__typeof__(&bfMatBlockDiagGetType))bfMatBlockDiagGetType,
  .GetNumRows = (__typeof__(&bfMatBlockDiagGetNumRows))bfMatBlockDiagGetNumRows,
  .GetNumCols = (__typeof__(&bfMatBlockDiagGetNumCols))bfMatBlockDiagGetNumCols,
  .ScaleCols = (__typeof__(&bfMatBlockDiagScaleCols))bfMatBlockDiagScaleCols,
  .Mul = (__typeof__(&bfMatBlockDiagMul))bfMatBlockDiagMul,
  .Negate = (__typeof__(&bfMatBlockDiagNegate))bfMatBlockDiagNegate,
};

BfMat *bfMatBlockDiagCopy(BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfMatBlock const *matBlock = NULL;
  BfMatBlockDiag const *matBlockDiag = NULL;
  BfMatBlockDiag *matBlockDiagCopy = NULL;
  BfMat const *block = NULL;
  BfMat *blockCopy = NULL;

  matBlock = bfMatConstToMatBlockConst(mat);
  HANDLE_ERROR();

  matBlockDiag = bfMatConstToMatBlockDiagConst(mat);
  HANDLE_ERROR();

  BfSize numBlocks = bfMatBlockDiagNumBlocks(matBlock);

  matBlockDiagCopy = bfMatBlockDiagNew();
  HANDLE_ERROR();

  bfMatBlockDiagInit(matBlockDiagCopy, numBlocks, numBlocks);
  HANDLE_ERROR();

  for (BfSize i = 0; i <= numBlocks; ++i)
    matBlockDiagCopy->super.rowOffset[i] = matBlock->rowOffset[i];

  for (BfSize j = 0; j <= numBlocks; ++j)
    matBlockDiagCopy->super.colOffset[j] = matBlock->colOffset[j];

  for (BfSize i = 0; i < numBlocks; ++i) {
    block = bfMatBlockDiagGetBlockConst(matBlockDiag, i);
    blockCopy = bfMatCopy(block);
    HANDLE_ERROR();
    bfMatBlockDiagSetBlock(matBlockDiagCopy, i, blockCopy);
  }

  END_ERROR_HANDLING() {
    bfMatDelete(&blockCopy);
    bfMatBlockDiagDeinitAndDealloc(&matBlockDiagCopy);
  }

  return bfMatBlockDiagToMat(matBlockDiagCopy);
}

BfVec *bfMatBlockDiagGetRowCopy(BfMat const *mat, BfSize i) {
  BEGIN_ERROR_HANDLING();

  BfVec *rowCopy = NULL;
  BfMatBlock const *matBlock = NULL;
  BfMatBlockDiag const *matBlockDiag = NULL;
  BfMat const *block = NULL;

  matBlock = bfMatConstToMatBlockConst(mat);
  HANDLE_ERROR();

  matBlockDiag = bfMatConstToMatBlockDiagConst(mat);
  HANDLE_ERROR();

  BfSize numRows = bfMatGetNumRows(mat);
  if (i >= numRows)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numRowBlocks = bfMatBlockGetNumRowBlocks(matBlock);

  BfSize k = 0, i0, i1;
  for (; k < numRowBlocks; ++k) {
    i0 = matBlockDiag->super.rowOffset[k];
    i1 = matBlockDiag->super.rowOffset[k + 1];
    if (i0 <= i && i < i1)
      break;
  }

  block = BLOCK(matBlockDiag, k);

  rowCopy = bfMatGetRowCopy(block, i - i0);
  HANDLE_ERROR();

  /* Pad what we have so far on either side to make sure it matches
   * the width of `mat` */

  BfSize j0 = matBlockDiag->super.colOffset[k];

  BfVec *cat = NULL;
  BfVecZero *zeroPad = NULL;

  if (j0 > 0) {
    zeroPad = bfVecZeroNew();
    HANDLE_ERROR();

    bfVecZeroInit(zeroPad, j0);
    HANDLE_ERROR();

    cat = bfVecConcat(bfVecZeroToVec(zeroPad), rowCopy);
    HANDLE_ERROR();

    bfVecZeroDeinitAndDealloc(&zeroPad);
    bfVecDelete(&rowCopy);

    rowCopy = cat;
  }

  BfSize numCols = bfMatGetNumCols(mat);

  if (rowCopy->size < numCols) {
    zeroPad = bfVecZeroNew();
    HANDLE_ERROR();

    bfVecZeroInit(zeroPad, numCols - rowCopy->size);
    HANDLE_ERROR();

    cat = bfVecConcat(rowCopy, bfVecZeroToVec(zeroPad));
    HANDLE_ERROR();

    bfVecZeroDeinitAndDealloc(&zeroPad);
    bfVecDelete(&rowCopy);

    rowCopy = cat;
  }

  assert(rowCopy->size == numCols);

  END_ERROR_HANDLING() {}

  return rowCopy;
}

void bfMatBlockDiagDelete(BfMat **mat) {
  bfMatBlockDiagDeinitAndDealloc((BfMatBlockDiag **)mat);
}

BfType bfMatBlockDiagGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_BLOCK_DIAG;
}

BfSize bfMatBlockDiagGetNumRows(BfMat const *mat) {
  BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
  return bfMatIsTransposed(mat) ?
    matBlock->colOffset[bfMatBlockGetNumColBlocks(matBlock)] :
    matBlock->rowOffset[bfMatBlockGetNumRowBlocks(matBlock)];
}

BfSize bfMatBlockDiagGetNumCols(BfMat const *mat) {
  BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
  return bfMatIsTransposed(mat) ?
    matBlock->rowOffset[bfMatBlockGetNumRowBlocks(matBlock)] :
    matBlock->colOffset[bfMatBlockGetNumColBlocks(matBlock)];
}

void bfMatBlockDiagScaleCols(BfMat *mat, BfVec const *vec) {
  BEGIN_ERROR_HANDLING();

  BfMatBlock *matBlock = bfMatToMatBlock(mat);
  HANDLE_ERROR();

  BfSize numBlocks = bfMatBlockDiagNumBlocks(matBlock);

  for (BfSize k = 0; k < numBlocks; ++k) {
    BfSize j0 = matBlock->colOffset[k];
    BfSize j1 = matBlock->colOffset[k + 1];
    BfMat *block = matBlock->block[k];
    BfVec *subvec = bfVecGetSubvecCopy(vec, j0, j1);
    bfMatScaleCols(block, subvec);
    bfVecDelete(&subvec);
  }

  END_ERROR_HANDLING() {}
}

BfMat *bfMatBlockDiagMul(BfMat const *mat, BfMat const *other) {
  BEGIN_ERROR_HANDLING();

  BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);

  BfSize numRows = bfMatGetNumRows(mat);
  BfSize numCols = bfMatGetNumCols(other);
  BfSize numBlocks = bfMatBlockDiagNumBlocks(matBlock);

  BfMat *result = NULL;
  BfMat *block = NULL;
  BfMat *op2Rows = NULL;
  BfMat *resultRows = NULL;

  result = bfMatEmptyLike(other, numRows, numCols);
  HANDLE_ERROR();

  for (BfSize i = 0, i0, i1, j0, j1; i < numBlocks; ++i) {
    i0 = matBlock->rowOffset[i];
    i1 = matBlock->rowOffset[i + 1];
    j0 = matBlock->colOffset[i];
    j1 = matBlock->colOffset[i + 1];

    block = matBlock->block[i];
    op2Rows = bfMatGetRowRange((BfMat *)other, j0, j1);
    resultRows = bfMatMul(block, op2Rows);
    bfMatSetRowRange(result, i0, i1, resultRows);

    bfMatDelete(&resultRows);
  }

  END_ERROR_HANDLING()
    bfMatDelete(&result);

  return result;
}

void bfMatBlockDiagNegate(BfMat *mat) {
  BEGIN_ERROR_HANDLING();

  BfMatBlock *matBlock = NULL;
  BfMatBlockDiag *matBlockDiag = NULL;

  matBlock = bfMatToMatBlock(mat);
  HANDLE_ERROR();

  BfSize numBlocks = bfMatBlockNumBlocks(matBlock);

  matBlockDiag = bfMatToMatBlockDiag(mat);
  HANDLE_ERROR();

  for (BfSize i = 0; i < numBlocks; ++i)
    bfMatNegate(bfMatBlockDiagGetBlock(matBlockDiag, i));

  END_ERROR_HANDLING() {}
}

/** Interface: MatBlock */

static BfMatBlockVtable MAT_BLOCK_VTABLE = {
  .NumBlocks = (__typeof__(&bfMatBlockDiagNumBlocks))bfMatBlockDiagNumBlocks,
  .GetBlockCopy = (__typeof__(&bfMatBlockGetBlockCopy))bfMatBlockDiagGetBlockCopy,
};

BfSize bfMatBlockDiagNumBlocks(BfMatBlock const *matBlock) {
  BfMat const *mat = bfMatBlockConstToMatConst(matBlock);
  return mat->numRows < mat->numCols ? mat->numRows : mat->numCols;
}

BfMat *bfMatBlockDiagGetBlockCopy(BfMatBlockDiag const *matBlockDiag, BfSize i, BfSize j) {
  BEGIN_ERROR_HANDLING();

  BfMatBlock const *matBlock = bfMatBlockDiagConstToMatBlockConst(matBlockDiag);
  HANDLE_ERROR();

  BfMat *block = NULL;

  BfSize m = bfMatBlockGetNumBlockRows(matBlock, i);
  BfSize n = bfMatBlockGetNumBlockCols(matBlock, j);

  /* Copy the ith block if this is a diagonal block. */
  if (i == j) {
    if (matBlock->block[i] == NULL)
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

    block = bfMatCopy(matBlock->block[i]);
    HANDLE_ERROR();
  }

  /* Otherwise, return a new zero block. */
  else {
    BfMatZero *matZero = bfMatZeroNew();
    HANDLE_ERROR();

    bfMatZeroInit(matZero, m, n);
    HANDLE_ERROR();
  }

  END_ERROR_HANDLING() {
    bfMatDelete(&block);
  }

  return block;
}

/** Upcasting: */

BfMat *bfMatBlockDiagToMat(BfMatBlockDiag *matBlockDiag) {
  return &matBlockDiag->super.super;
}

BfMatBlock *bfMatBlockDiagToMatBlock(BfMatBlockDiag *matBlockDiag) {
  return &matBlockDiag->super;
}

BfMatBlock const *bfMatBlockDiagConstToMatBlockConst(BfMatBlockDiag const *matBlockDiag) {
  return &matBlockDiag->super;
}

/** Downcasting: */

BfMatBlockDiag *bfMatToMatBlockDiag(BfMat *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_BLOCK_DIAG)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatBlockDiag *)mat;
  }
}

BfMatBlockDiag const *bfMatConstToMatBlockDiagConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_BLOCK_DIAG)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatBlockDiag const *)mat;
  }
}

/** Implementation: MatBlockDiag */

BfMatBlockDiag *bfMatBlockDiagNew() {
  BEGIN_ERROR_HANDLING();

  BfMatBlockDiag *mat = malloc(sizeof(BfMatBlockDiag));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return mat;
}

BfMatBlockDiag *bfMatBlockDiagNewFromBlocks(BfPtrArray *blocks) {
  BEGIN_ERROR_HANDLING();

  BfMatBlockDiag *matBlockDiag = bfMatBlockDiagNew();
  HANDLE_ERROR();

  BfSize numBlocks = bfPtrArraySize(blocks);

  bfMatBlockInit(&matBlockDiag->super, &MAT_VTABLE, &MAT_BLOCK_VTABLE,
                 numBlocks, numBlocks, numBlocks);
  HANDLE_ERROR();

  for (BfSize i = 0; i < numBlocks; ++i) {
    BLOCK(matBlockDiag, i) = bfPtrArrayGet(blocks, i);
  }

  ROW_OFFSET(matBlockDiag, 0) = 0;
  COL_OFFSET(matBlockDiag, 0) = 0;
  for (BfSize i = 0; i < numBlocks; ++i) {
    BfMat const *block = bfPtrArrayGet(blocks, i);
    ROW_OFFSET(matBlockDiag, i + 1) = bfMatGetNumRows(block);
    COL_OFFSET(matBlockDiag, i + 1) = bfMatGetNumCols(block);
  }
  bfSizeRunningSum(numBlocks + 1, &ROW_OFFSET(matBlockDiag, 0));
  bfSizeRunningSum(numBlocks + 1, &COL_OFFSET(matBlockDiag, 0));

  END_ERROR_HANDLING() {
    assert(false);
  }

  return matBlockDiag;
}

void bfMatBlockDiagInit(BfMatBlockDiag *mat, BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  BfSize numBlocks = numRows < numCols ? numRows : numCols;

  bfMatBlockInit(&mat->super, &MAT_VTABLE, &MAT_BLOCK_VTABLE, numBlocks, numRows, numCols);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {}
}

void bfMatBlockDiagDeinit(BfMatBlockDiag *mat) {
  bfMatBlockDeinit(&mat->super);
}

void bfMatBlockDiagDealloc(BfMatBlockDiag **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatBlockDiagDeinitAndDealloc(BfMatBlockDiag **mat) {
  bfMatBlockDiagDeinit(*mat);
  bfMatBlockDiagDealloc(mat);
}

BfMat *bfMatBlockDiagGetBlock(BfMatBlockDiag *matBlockDiag, BfSize i) {
  BEGIN_ERROR_HANDLING();

  BfMatBlock *matBlock = bfMatBlockDiagToMatBlock(matBlockDiag);
  BfMat *block = NULL;

  BfSize numRowBlocks = bfMatBlockGetNumRowBlocks(matBlock);
  if (i >= numRowBlocks)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numColBlocks = bfMatBlockGetNumColBlocks(matBlock);
  if (i >= numColBlocks)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  block = matBlockDiag->super.block[i];

  END_ERROR_HANDLING()
    block = NULL;

  return block;
}

BfMat const *bfMatBlockDiagGetBlockConst(BfMatBlockDiag const *matBlockDiag, BfSize i) {
  BEGIN_ERROR_HANDLING();

  BfMatBlock const *matBlock = bfMatBlockDiagConstToMatBlockConst(matBlockDiag);
  BfMat const *block = NULL;

  BfSize numRowBlocks = bfMatBlockGetNumRowBlocks(matBlock);
  if (i >= numRowBlocks)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numColBlocks = bfMatBlockGetNumColBlocks(matBlock);
  if (i >= numColBlocks)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  block = matBlockDiag->super.block[i];

  END_ERROR_HANDLING()
    block = NULL;

  return block;
}

void bfMatBlockDiagSetBlock(BfMatBlockDiag *matBlockDiag, BfSize i, BfMat *mat) {
  BEGIN_ERROR_HANDLING();

  BfMatBlock const *matBlock = bfMatBlockDiagConstToMatBlockConst(matBlockDiag);

  BfSize numRowBlocks = bfMatBlockGetNumRowBlocks(matBlock);
  if (i >= numRowBlocks)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numColBlocks = bfMatBlockGetNumColBlocks(matBlock);
  if (i >= numColBlocks)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  matBlockDiag->super.block[i] = mat;

  END_ERROR_HANDLING() {}
}
