#include <bf/mat_block_diag.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/vec_zero.h>

#define INTERFACE BF_INTERFACE_Mat
BF_DEFINE_VTABLE(Mat, MatBlockDiag)
#undef INTERFACE

#define INTERFACE BF_INTERFACE_MatBlock
BF_DEFINE_VTABLE(MatBlock, MatBlockDiag)
#undef INTERFACE

/** Helper macros: */

#define NUM_ROW_BLOCKS(mat) mat->super.super.numRows
#define NUM_COL_BLOCKS(mat) mat->super.super.numCols
#define NUM_BLOCKS(mat) NUM_ROW_BLOCKS(mat)*NUM_COL_BLOCKS(mat)
#define ROW_OFFSET(mat, k) mat->super.rowOffset[k]
#define COL_OFFSET(mat, k) mat->super.colOffset[k]
#define NUM_BLOCK_ROWS(mat, k) ROW_OFFSET(mat, k + 1) - ROW_OFFSET(mat, k)
#define NUM_BLOCK_COLS(mat, k) COL_OFFSET(mat, k + 1) - COL_OFFSET(mat, k)
#define BLOCK(mat, k) mat->super.block[k]

/** BfMat interface: */

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

BF_STUB(BfMat *, MatBlockDiagGetView, BfMat *)

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

BF_STUB(BfVec *, MatBlockDiagGetRowView, BfMat *, BfSize)
BF_STUB(BfVec *, MatBlockDiagGetColView, BfMat *, BfSize)
BF_STUB(BfVec *, MatBlockDiagGetColRangeView, BfMat *, BfSize, BfSize, BfSize)

void bfMatBlockDiagDelete(BfMat **mat) {
  bfMatBlockDiagDeinitAndDealloc((BfMatBlockDiag **)mat);
}

BF_STUB(BfMat *, MatBlockDiagEmptyLike, BfMat const *, BfSize, BfSize)
BF_STUB(BfMat *, MatBlockDiagZerosLike, BfMat const *, BfSize, BfSize)

BfType bfMatBlockDiagGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_BLOCK_DIAG;
}

bool bfMatBlockDiagInstanceOf(BfMat const *mat, BfType type) {
  return bfTypeDerivedFrom(bfMatGetType(mat), type);
}

BF_STUB(BfSize, MatBlockDiagNumBytes, BfMat const *)
BF_STUB(void, MatBlockDiagSave, BfMat const *, char const *)
BF_STUB(void, MatBlockDiagPrint, BfMat const *, FILE *)

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

BF_STUB(void, MatBlockDiagSetRow, BfMat *, BfSize, BfVec const *)
BF_STUB(void, MatBlockDiagSetCol, BfMat *, BfSize, BfVec const *)
BF_STUB(void, MatBlockDiagSetColRange, BfMat *, BfSize, BfSize, BfSize, BfVec const *)
BF_STUB(BfMat *, MatBlockDiagGetRowRange, BfMat *, BfSize, BfSize)
BF_STUB(BfMat *, MatBlockDiagGetColRange, BfMat *, BfSize, BfSize)
BF_STUB(BfMat *, MatBlockDiagGetRowRangeCopy, BfMat const *, BfSize, BfSize)
BF_STUB(BfMat *, MatBlockDiagGetColRangeCopy, BfMat const *, BfSize, BfSize)
BF_STUB(void, MatBlockDiagSetRowRange, BfMat *, BfSize, BfSize, BfMat const *)
BF_STUB(void, MatBlockDiagPermuteRows, BfMat *, BfPerm const *)
BF_STUB(void, MatBlockDiagPermuteCols, BfMat *, BfPerm const *)
BF_STUB(BfVec *, MatBlockDiagRowDists, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatBlockDiagColDists, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatBlockDiagColDots, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatBlockDiagColNorms, BfMat const *)
BF_STUB(void, MatBlockDiagScaleRows, BfMat *, BfVec const *)

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

BF_STUB(BfVec *, MatBlockDiagSumCols, BfMat const *)
BF_STUB(void, MatBlockDiagAddInplace, BfMat *, BfMat const *)
BF_STUB(void, MatBlockDiagAddDiag, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatBlockDiagSub, BfMat const *, BfMat const *)
BF_STUB(void, MatBlockDiagSubInplace, BfMat *, BfMat const *)

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

BF_STUB(BfVec *, MatBlockDiagMulVec, BfMat const *, BfVec const *)
BF_STUB(void, MatBlockDiagMulInplace, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatBlockDiagSolveLU, BfMat const *, BfMat const *)
BF_STUB(BfMat *, MatBlockDiagLstSq, BfMat const *, BfMat const *)
BF_STUB(BfMat *, MatBlockDiagGetGivensRotation, BfVec const *, BfSize, BfSize)
BF_STUB(bool, MatBlockDiagIsUpperTri, BfMat const *)
BF_STUB(BfVec *, MatBlockDiagForwardSolveVec, BfMat const *, BfVec const *)
BF_STUB(BfVec *, MatBlockDiagBackwardSolveVec, BfMat const *, BfVec const *)
BF_STUB(bool, MatBlockDiagIsZero, BfMat const *)

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

BF_STUB(BfMat *, MatBlockDiagToType, BfMat const *, BfType)
BF_STUB(BfMat *, MatBlockDiagCholesky, BfMat const *)

/** Interface: MatBlock */

BfSize bfMatBlockDiagNumBlocks(BfMatBlock const *matBlock) {
  BfMat const *mat = bfMatBlockConstToMatConst(matBlock);
  return mat->numRows < mat->numCols ? mat->numRows : mat->numCols;
}

BF_STUB(BfSize, MatBlockDiagGetNumRowBlocks, BfMatBlock const *)
BF_STUB(BfSize, MatBlockDiagGetNumColBlocks, BfMatBlock const *)
BF_STUB(BfSize, MatBlockDiagGetNumBlockRows, BfMatBlock const *, BfSize)
BF_STUB(BfSize, MatBlockDiagGetNumBlockCols, BfMatBlock const *, BfSize)

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

/** Implementation: */

BfMatBlockDiag *bfMatBlockDiagNew() {
  BEGIN_ERROR_HANDLING();

  BfMatBlockDiag *mat = malloc(sizeof(BfMatBlockDiag));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return mat;
}

void bfMatBlockDiagInit(BfMatBlockDiag *mat, BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  BfSize numBlocks = numRows < numCols ? numRows : numCols;

  bfMatBlockInit(&mat->super, &MatVtbl, &MatBlockVtbl, numBlocks, numRows, numCols);
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
