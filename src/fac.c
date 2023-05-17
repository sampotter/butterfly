#include <bf/fac.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/indexed_mat.h>
#include <bf/linalg.h>
#include <bf/mat_block_coo.h>
#include <bf/mat_block_dense.h>
#include <bf/mat_block_diag.h>
#include <bf/mat_identity.h>
#include <bf/mat_product.h>
#include <bf/mem.h>
#include <bf/size_array.h>
#include <bf/tree_iter.h>

void bfFacDeinit(BfFac *fac) {
  fac->colNode = NULL;

  bfConstNodeArrayDeinit(&fac->rowNodes);

  bfMatDelete(&fac->Psi);
  for (BfSize i = 0; i < fac->numW; ++i)
    bfMatDelete(&fac->W[i]);
}

void bfFacDealloc(BfFac **facHandle) {
  bfMemFree(*facHandle);
  *facHandle = NULL;
}

void bfFacDelete(BfFac **facHandle) {
  bfFacDeinit(*facHandle);
  bfFacDealloc(facHandle);
}

BfMat *bfFacGetMat(BfFac const *fac) {
  BF_ERROR_BEGIN();

  BfMatProduct *matProduct = bfFacGetMatProduct(fac);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    BF_DIE();
  }

  return bfMatProductToMat(matProduct);
}

BfMatProduct *bfFacGetMatProduct(BfFac const *fac) {
  BF_ERROR_BEGIN();

  BfMatProduct *matProduct = bfMatProductNew();
  HANDLE_ERROR();

  bfMatProductInit(matProduct);
  HANDLE_ERROR();

  bfMatProductPostMultiply(matProduct, fac->Psi);
  HANDLE_ERROR();

  for (BfSize k = 0; k < fac->numW; ++k) {
    bfMatProductPostMultiply(matProduct, fac->W[k]);
    HANDLE_ERROR();
  }

  END_ERROR_HANDLING() {
    BF_DIE();
  }

  return matProduct;
}

BfSize bfFacGetNumBytes(BfFac const *fac) {
  BfSize numBytes = bfMatNumBytes(fac->Psi);
  for (BfSize k = 0; k < fac->numW; ++k)
    numBytes += bfMatNumBytes(fac->W[k]);
  return numBytes;
}

BfFac *makeLeafNodePartialFac(BfTreeNode const *colNode,
                              BfPtrArray *PsiBlocks,
                              BfPtrArray *WBlocks) {
  BF_ERROR_BEGIN();

  BfFac *fac = bfMemAlloc(1, sizeof(BfFac));
  if (fac == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  fac->colNode = colNode;

  bfConstNodeArrayInitWithDefaultCapacity(&fac->rowNodes);
  HANDLE_ERROR();

  BfMatBlockDiag *Psi = bfMatBlockDiagNewFromBlocks(PsiBlocks, BF_POLICY_COPY);
  HANDLE_ERROR();

  fac->Psi = bfMatBlockDiagToMat(Psi);

  fac->numW = 1;

  fac->W = bfMemAlloc(1, sizeof(BfMat *));
  if (fac->W == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfMatBlockDense *W0 = bfMatBlockDenseNewColFromBlocks(WBlocks, BF_POLICY_COPY);
  HANDLE_ERROR();

  fac->W[0] = bfMatBlockDenseToMat(W0);

  END_ERROR_HANDLING() {
    fac = NULL;
  }

  return fac;
}

void setPartialFacWBlock(BfFac *fac, BfSize k, BfMat *W) {
  BF_ASSERT(k < fac->numW);
  BF_ASSERT(W != NULL);
  fac->W[k] = W;
}

BfSize partialFacGetNumRows(BfFac const *fac) {
  return bfMatGetNumRows(fac->Psi);
}

BfVec *partialFacMulVec(BfFac const *fac, BfVec const *x) {
  BF_ASSERT(fac->numW > 0);
  BfSize k = fac->numW - 1;
  BfVec *y = NULL;
  BfVec *yPrev = bfMatMulVec(fac->W[k], x);
  while (k > 0) {
    y = bfMatMulVec(fac->W[--k], yPrev);
    bfVecDelete(&yPrev);
    yPrev = y;
  }
  y = bfMatMulVec(fac->Psi, yPrev);
  bfVecDelete(&yPrev);
  return y;
}

static void appendIndexedPsiSubblock(BfPtrArray *indexedPsiSubblocks,
                                     BfSize i0, BfSize j0, BfMat *mat) {
  BF_ERROR_BEGIN();

  BfIndexedMat *indexedMat = bfMemAlloc(1, sizeof(BfIndexedMat));
  if (indexedMat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  indexedMat->i0 = i0;
  indexedMat->j0 = j0;
  indexedMat->mat = mat;

  bfPtrArrayAppend(indexedPsiSubblocks, indexedMat);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    BF_ASSERT(false);
  }
}

static void getIndexedPsiSubblocksInRowRangeRec(BfMat *mat,
                                                BfSize i0Sel, BfSize i1Sel,
                                                BfSize i0Parent, BfSize j0Parent,
                                                BfPtrArray *indexedPsiSubblocks) {
  BF_ERROR_BEGIN();

  BfType type = bfMatGetType(mat);

  if (type == BF_TYPE_MAT_BLOCK_DIAG) {
    BfMatBlockDiag const *matBlockDiag = bfMatConstToMatBlockDiagConst(mat);
    BfSize numBlocks = bfMatBlockDiagNumBlocks(matBlockDiag);
    for (BfSize k = 0; k < numBlocks; ++k) {
      BfMat *block = matBlockDiag->super.block[k];

      /* Compute row range for current block */
      BfSize i0 = i0Parent + matBlockDiag->super.rowOffset[k];
      BfSize i1 = i0 + bfMatGetNumRows(block);

      /* Get column offset for current block */
      BfSize j0 = j0Parent + matBlockDiag->super.colOffset[k];

      /* If there's some overlap, continue recursively for this block: */
      if (!(i1Sel <= i0 || i1 <= i0Sel)) {
        getIndexedPsiSubblocksInRowRangeRec(
          block, i0Sel, i1Sel, i0, j0, indexedPsiSubblocks);
        HANDLE_ERROR();
      }
    }
  }

  else if (type == BF_TYPE_MAT_DENSE_REAL ||
           type == BF_TYPE_MAT_BLOCK_COO ||
           type == BF_TYPE_MAT_BLOCK_DENSE ||
           type == BF_TYPE_MAT_IDENTITY)
    appendIndexedPsiSubblock(indexedPsiSubblocks, i0Parent, j0Parent, mat);

  else RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  END_ERROR_HANDLING() {
    BF_ASSERT(false);
  }
}

BfPtrArray *getIndexedPsiSubblocksInRowRange(BfFac const *fac, BfSize i0Sel, BfSize i1Sel) {
  BF_ERROR_BEGIN();

  BfPtrArray *indexedPsiSubblocks = bfPtrArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  getIndexedPsiSubblocksInRowRangeRec(fac->Psi, i0Sel, i1Sel, 0, 0, indexedPsiSubblocks);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    BF_ASSERT(false);
  }

  return indexedPsiSubblocks;
}

void getPsiAndW0BlocksByRowNodeForPartialFac(BfFac const *fac,
                                             BfTreeNode const *rowNode,
                                             BfMat **PsiBlockPtr,
                                             BfMat **W0BlockPtr) {
  BF_ERROR_BEGIN();

  if (PsiBlockPtr == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (W0BlockPtr == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMat const *W0 = fac->W[0];

  BfMat *PsiBlock = NULL;
  BfMat *W0Block = NULL;

  BfPtrArray PsiSubblocks;
  bfInitPtrArrayWithDefaultCapacity(&PsiSubblocks);
  HANDLE_ERROR();

  BfPtrArray W0Subblocks;
  bfInitPtrArrayWithDefaultCapacity(&W0Subblocks);
  HANDLE_ERROR();

  BfSize i0 = bfTreeNodeGetFirstIndex(rowNode);
  BfSize i1 = bfTreeNodeGetLastIndex(rowNode);

  /* Verify that `rowNode` is the common parent for the subspan of
   * nodes in `fac->rowNodes` whose index ranges are a subset
   * of [i0, i1). */
#if BF_DEBUG
  for (BfSize k = 0; k < bfConstNodeArraySize(&fac->rowNodes); ++k) {
    BfTreeNode const *otherRowNode = bfConstNodeArrayGet(&fac->rowNodes, k);
    BfSize i0_ = bfTreeNodeGetFirstIndex(otherRowNode);
    BfSize i1_ = bfTreeNodeGetLastIndex(otherRowNode);
    if (i0 <= i0_ && i1_ <= i1 && rowNode != otherRowNode)
      BF_ASSERT(bfTreeNodeIsDescendant(otherRowNode, rowNode));
  }
#endif

  /* Find all of the leaf subblocks in Psi and their offsets. */
  BfPtrArray *indexedPsiSubblocks =
    getIndexedPsiSubblocksInRowRange(fac, i0, i1);

  BfSize numSubblocks = bfPtrArraySize(indexedPsiSubblocks);
  BF_ASSERT(numSubblocks > 0);

  BfSize i1Prev = BF_SIZE_BAD_VALUE;
  BfSize j1Prev = BF_SIZE_BAD_VALUE;

  /* Iterate over each leaf Psi subblock, grab a copy of it an the
   * corresponding W0 row range to our running lists of blocks */
  for (BfSize k = 0; k < numSubblocks; ++k) {
    BfIndexedMat *indexedMat = bfPtrArrayGet(indexedPsiSubblocks, k);
    BfMat *PsiSubblock = indexedMat->mat;

    /* Get row span for block */
    BfSize i0_ = indexedMat->i0;
    BfSize i1_ = i0_ + bfMatGetNumRows(PsiSubblock);
    BF_ASSERT(i0 <= i0_ && i1_ <= i1);

    /* Row spans should be adjacent */
    BF_ASSERT(i1Prev == BF_SIZE_BAD_VALUE || i0_ == i1Prev);
    i1Prev = i1_;

    /* Get column range for block */
    BfSize j0_ = indexedMat->j0;
    BfSize j1_ = j0_ + bfMatGetNumCols(PsiSubblock);

    /* Column spaces should be adjacent */
    BF_ASSERT(j1Prev == BF_SIZE_BAD_VALUE || j0_ == j1Prev);
    j1Prev = j1_;

    /* Get a copy of the current diagonal Psi block */
    BfMat *PsiSubblockCopy = bfMatCopy(PsiSubblock);
    HANDLE_ERROR();

    /* ... and append it to the running array of Psi subblocks */
    bfPtrArrayAppend(&PsiSubblocks, PsiSubblockCopy);
    HANDLE_ERROR();

    /* Get a copy of the current W row subblock */
    BfMat *W0Subblock = bfMatGetRowRangeCopy(W0, j0_, j1_);
    HANDLE_ERROR();

    /* .. and append it to the running array of W subblocks */
    bfPtrArrayAppend(&W0Subblocks, W0Subblock);
    HANDLE_ERROR();
  }

  BF_ASSERT(bfPtrArraySize(&PsiSubblocks) == numSubblocks);
  BF_ASSERT(bfPtrArraySize(&W0Subblocks) == numSubblocks);

  /* If we only found one Psi and W subblock each, we should just
   * return these directly instead of wrapping them in a MatBlockDiag
   * and MatBlockCoo, respectively (what is done otherwise). */
  if (numSubblocks == 1) {
    PsiBlock = bfMatGet(bfPtrArrayGet(&PsiSubblocks, 0), BF_POLICY_STEAL);
    W0Block = bfMatGet(bfPtrArrayGet(&W0Subblocks, 0), BF_POLICY_STEAL);
  } else {
    /* Diagonally concatenate the Psi subblocks we found */
    PsiBlock = bfMatBlockDiagToMat(
      bfMatBlockDiagNewFromBlocks(&PsiSubblocks, BF_POLICY_STEAL));
    HANDLE_ERROR();

    /* ... and vertically concatenate the W row subblocks */
    W0Block = bfMatBlockDenseToMat(
      bfMatBlockDenseNewColFromBlocks(&W0Subblocks, BF_POLICY_STEAL));
    HANDLE_ERROR();
  }

  END_ERROR_HANDLING() {
    BF_ASSERT(false); // T_T
  }

  for (BfSize i = 0; i < numSubblocks; ++i) {
    BfMat *PsiSubblock = bfPtrArrayGet(&PsiSubblocks, i);
    bfMatDelete(&PsiSubblock);
  }
  bfPtrArrayDeinit(&PsiSubblocks);

  for (BfSize i = 0; i < numSubblocks; ++i) {
    BfMat *W0Subblock = bfPtrArrayGet(&W0Subblocks, i);
    bfMatDelete(&W0Subblock);
  }
  bfPtrArrayDeinit(&W0Subblocks);

  for (BfSize k = 0; k < numSubblocks; ++k) {
    BfIndexedMat *indexedMat = bfPtrArrayGet(indexedPsiSubblocks, k);
    bfMemFree(indexedMat);
  }
  bfPtrArrayDelete(&indexedPsiSubblocks);

  BF_ASSERT(PsiBlock != NULL);
  BF_ASSERT(W0Block != NULL);
  BF_ASSERT(bfTreeNodeGetNumPoints(rowNode) == bfMatGetNumRows(PsiBlock));
  BF_ASSERT(bfMatGetNumCols(PsiBlock) == bfMatGetNumRows(W0Block));
  BF_ASSERT(bfMatGetNumCols(W0Block) == bfMatGetNumCols(fac->W[0]));

  *PsiBlockPtr = PsiBlock;
  *W0BlockPtr = W0Block;
}

static bool facsHaveContiguousRowSpans(BfPtrArray const *facs) {
  for (BfSize k = 0; k < bfPtrArraySize(facs); ++k) {
    BfFac const *fac = bfPtrArrayGet(facs, k);
    if (!bfConstNodeArrayIsContiguous(&fac->rowNodes))
      return false;
  }
  return true;
}

// TODO: this should be replaced with a function on the NodeSpan
// type. We should only ever have to call it in mergeAndSplit.
bool partialFacsHaveSameRowSpan(BfPtrArray const *facs) {
  BF_ASSERT(bfPtrArraySize(facs) >= 2);

  if (!facsHaveContiguousRowSpans(facs)) {
    // TODO: this case is a little annoying and complicated and
    // doesn't matter for us right now. Handle it later if we need to
    // for some reason. (A later comment: unnecessary to handle? When
    // would we ever have a Fac which doesn't have a contiguous row
    // span...?)
    BF_ASSERT(false);
  }

  BfFac *fac;
  BfTreeNode const *rowNode;

  BfSize i0, i1;

  bfPtrArrayGetFirst(facs, (BfPtr *)&fac);

  rowNode = bfConstNodeArrayGetFirst(&fac->rowNodes);
  i0 = bfTreeNodeGetFirstIndex(rowNode);

  rowNode = bfConstNodeArrayGetLast(&fac->rowNodes);
  i1 = bfTreeNodeGetLastIndex(rowNode);

  for (BfSize k = 1; k < bfPtrArraySize(facs); ++k) {
    BfSize i0_, i1_;

    bfPtrArrayGetFirst(facs, (BfPtr *)&fac);

    rowNode = bfConstNodeArrayGetFirst(&fac->rowNodes);
    i0_ = bfTreeNodeGetFirstIndex(rowNode);

    rowNode = bfConstNodeArrayGetLast(&fac->rowNodes);
    i1_ = bfTreeNodeGetLastIndex(rowNode);

    if (i0 != i0_ || i1 != i1_)
      return false;
  }

  return true;
}

/* Given a `PtrArray` filled with `PartialFac`s, get the first row
 * tree node in each `PartialFac`, collect them together into a
 * `ConstPtrArray`, and return it. */
BfConstNodeArray getFirstRowNodes(BfPtrArray const *facs) {
  BF_ERROR_BEGIN();

  BfConstNodeArray rowNodes;
  bfConstNodeArrayInitWithDefaultCapacity(&rowNodes);
  HANDLE_ERROR();

  for (BfSize k = 0; k < bfPtrArraySize(facs); ++k) {
    BfFac const *fac = bfPtrArrayGet(facs, k);

    if (bfConstNodeArrayIsEmpty(&fac->rowNodes))
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

    BfTreeNode const *rowNode = bfConstNodeArrayGetFirst(&fac->rowNodes);

    bfConstNodeArrayAppend(&rowNodes, rowNode);
  }

  END_ERROR_HANDLING() {
    bfConstNodeArrayDeinit(&rowNodes);
  }

  return rowNodes;
}

/* Given a `PtrArray` filled with `PartialFac`s, get the last row
 * tree node in each `PartialFac`, collect them together into a
 * `ConstPtrArray`, and return it. */
// TODO: this function is unnecessary (see the merge cut function below)
BfConstNodeArray getLastRowNodes(BfPtrArray const *facs) {
  BF_ERROR_BEGIN();

  BfConstNodeArray rowNodes;
  bfConstNodeArrayInitWithDefaultCapacity(&rowNodes);
  HANDLE_ERROR();

  for (BfSize k = 0; k < bfPtrArraySize(facs); ++k) {
    BfFac const *fac = bfPtrArrayGet(facs, k);

    if (bfConstNodeArrayIsEmpty(&fac->rowNodes))
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

    BfTreeNode const *rowNode = bfConstNodeArrayGetLast(&fac->rowNodes);

    bfConstNodeArrayAppend(&rowNodes, rowNode);
  }

  END_ERROR_HANDLING() {
    bfConstNodeArrayDeinit(&rowNodes);
  }

  return rowNodes;
}

BfConstNodeArray getRowNodesByFirstIndex(BfPtrArray const *facs, BfSize i0) {
  BF_ERROR_BEGIN();

  BfConstNodeArray rowNodes;
  bfConstNodeArrayInitWithDefaultCapacity(&rowNodes);
  HANDLE_ERROR();

  for (BfSize k = 0; k < bfPtrArraySize(facs); ++k) {
    BfFac const *fac = bfPtrArrayGet(facs, k);

    BfTreeNode const *rowNode = getNodeByFirstIndex(&fac->rowNodes, i0);
    if (rowNode == NULL)
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

    bfConstNodeArrayAppend(&rowNodes, rowNode);
    HANDLE_ERROR();
  }

  END_ERROR_HANDLING() {
    bfConstNodeArrayDeinit(&rowNodes);
  }

  return rowNodes;
}

BfConstNodeArray getMergeCut(BfPtrArray const *partialFacs) {
  BF_ERROR_BEGIN();

  BfSize numPartialFacs = bfPtrArraySize(partialFacs);

  if (numPartialFacs == 0)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* For the "merge cut" operation to be well defined, each of the
   * factorizations must have the same row span. */
  if (!partialFacsHaveSameRowSpan(partialFacs))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfConstNodeArray mergeCut;
  bfConstNodeArrayInitWithDefaultCapacity(&mergeCut);
  HANDLE_ERROR();

  BfConstNodeArray rowNodes = getFirstRowNodes(partialFacs);
  HANDLE_ERROR();

  BF_ASSERT(bfConstNodeArraySize(&rowNodes) == numPartialFacs);
  BF_ASSERT(nodesHaveSameFirstIndex(&rowNodes));

  BfConstNodeArray lastRowNodes = getLastRowNodes(partialFacs);
  HANDLE_ERROR();

  BfTreeNode const *rowNode = NULL;
  BfSize i1 = getMaxLastIndexForRowNodes(&rowNodes, &rowNode);

  bfConstNodeArrayAppend(&mergeCut, rowNode);
  HANDLE_ERROR();

  BfSize i1Final = getMaxLastIndexForRowNodes(&lastRowNodes, NULL);

  while (i1 != i1Final) {
    bfConstNodeArrayDeinit(&rowNodes);

    rowNodes = getRowNodesByFirstIndex(partialFacs, i1);
    HANDLE_ERROR();

    i1 = getMaxLastIndexForRowNodes(&rowNodes, &rowNode);
    HANDLE_ERROR();

    bfConstNodeArrayAppend(&mergeCut, rowNode);
    HANDLE_ERROR();

    BF_ASSERT(bfConstNodeArraySize(&rowNodes) == numPartialFacs);
  }

  END_ERROR_HANDLING() {
    bfConstNodeArrayDeinit(&mergeCut);
  }

  return mergeCut;
}

void getPsiAndW0BlocksByRowNode(BfPtrArray const *currentPartialFacs,
                               BfTreeNode const *rowNode,
                               BfMat **PsiPtr, BfMat **W0Ptr) {
  BF_ERROR_BEGIN();

  BfMatBlockDense *Psi = NULL;
  BfMatBlockDiag *W0 = NULL;

  if (PsiPtr == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  if (W0Ptr == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfPtrArray PsiBlocks;
  bfInitPtrArrayWithDefaultCapacity(&PsiBlocks);
  HANDLE_ERROR();

  BfPtrArray W0Blocks;
  bfInitPtrArrayWithDefaultCapacity(&W0Blocks);
  HANDLE_ERROR();

  for (BfSize k = 0; k < bfPtrArraySize(currentPartialFacs); ++k) {
    BfFac const *partialFac = bfPtrArrayGet(currentPartialFacs, k);

    /* Get diagonal psi blocks for kth partial fac and corresponding
     * index range by indexing with rowNode */
    BfMat *PsiBlock = NULL;
    BfMat *W0Block = NULL;
    getPsiAndW0BlocksByRowNodeForPartialFac(partialFac, rowNode, &PsiBlock, &W0Block);

    /* Append the new Psi block */
    bfPtrArrayAppend(&PsiBlocks, PsiBlock);
    HANDLE_ERROR();

    /* Append the new W0 block */
    bfPtrArrayAppend(&W0Blocks, W0Block);
    HANDLE_ERROR();
  }

  /* Horizontally concatenate together the Psi blocks we found to get
   * the leading "Psi*" block. */
  Psi = bfMatBlockDenseNewRowFromBlocks(&PsiBlocks, BF_POLICY_STEAL);
  HANDLE_ERROR();

  /* Diagonally concatenate together the W0 row blocks we found to get
   * the first column block of the permuted W0 factor. */
  W0 = bfMatBlockDiagNewFromBlocks(&W0Blocks, BF_POLICY_STEAL);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    BF_DIE();
  }

  for (BfSize i = 0; i < bfPtrArraySize(&PsiBlocks); ++i) {
    BfMat *PsiBlock = bfPtrArrayGet(&PsiBlocks, i);
    bfMatDelete(&PsiBlock);
  }
  bfPtrArrayDeinit(&PsiBlocks);

  for (BfSize i = 0; i < bfPtrArraySize(&W0Blocks); ++i) {
    BfMat *W0Block = bfPtrArrayGet(&W0Blocks, i);
    bfMatDelete(&W0Block);
  }
  bfPtrArrayDeinit(&W0Blocks);

  *PsiPtr = bfMatBlockDenseToMat(Psi);
  *W0Ptr = bfMatBlockDiagToMat(W0);

  BF_ASSERT(bfMatGetNumCols(*PsiPtr) == bfMatGetNumRows(*W0Ptr));
}

static bool getPsiAndW_skinny(BfMat const *block, BfMat **PsiPtr, BfMat **WPtr) {
  BF_ERROR_BEGIN();

  bool success = true;

  BfSize n = bfMatGetNumCols(block);

  BfMat *Psi = bfMatCopy(block);
  HANDLE_ERROR();

  BfMatIdentity *W = bfMatIdentityNew();
  HANDLE_ERROR();

  bfMatIdentityInit(W, n);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    success = false;

    bfMatDelete(&Psi);
    bfMatIdentityDeinitAndDealloc(&W);
  }

  *PsiPtr = Psi;
  *WPtr = bfMatIdentityToMat(W);

  return success;
}

static bool getPsiAndW_normal(BfMat const *block, BfReal tol, BfMat **PsiPtr, BfMat **WPtr) {
  BF_ERROR_BEGIN();

  /* Compute truncated SVD of the block. */
  BfMat *Psi = NULL;
  BfMatDiagReal *S = NULL;
  BfMat *W = NULL;
  BfTruncSpec truncSpec = {.usingTol = true, .tol = tol};
  bool success = bfGetTruncatedSvd(block, &Psi, &S, &W, &truncSpec, BF_BACKEND_LAPACK);
  HANDLE_ERROR();

  /* Scale W if we successfully compressed this block. */
  if (success) {
    BfVec *s = bfMatDiagRealGetVecView(S);
    HANDLE_ERROR();

    bfMatScaleRows(W, s);
    HANDLE_ERROR();
  }

  END_ERROR_HANDLING() {
    success = false;
  }

  if (success) {
    *PsiPtr = Psi;
    *WPtr = W;
  } else {
    bfMatDelete(&Psi);
    bfMatDelete(&W);
  }

  bfMatDiagRealDeinitAndDealloc(&S);

  return success;
}

bool getPsiAndW(BfFacSpec const *facSpec,
                BfMat const *mat, BfTreeNode const *rowNode,
                BfMat **PsiPtr, BfMat **WPtr) {
  BF_ERROR_BEGIN();

  /** NOTE! This function is actually a bit redundant... The other
   * place in this file where we check the extreme cases of having too
   * few rows or columns is `findEpsilonRankCutAndGetNewBlocks`. We
   * should do what we need to do to replace `getPsiAndW` with that
   * function. (I think this is the right idea, anyway...) */

  BfSize i0 = bfTreeNodeGetFirstIndex(rowNode);
  BfSize i1 = bfTreeNodeGetLastIndex(rowNode);

  /* TODO: should be something like: `bfMatGetRowRangeViewConst` */
  BfMat *block = bfMatGetRowRange((BfMat *)mat, i0, i1);
  HANDLE_ERROR();

  bool success = false;

  /* If there are too few columns in the current block, we don't
   * bother trying to compress. Instead, we just pass the current
   * block through for the Psi block and emit an identity matrix for
   * the new W block. */
  if (bfMatGetNumCols(mat) < facSpec->minNumCols)
    success = getPsiAndW_skinny(block, PsiPtr, WPtr);

  /* If there are too few rows, we pass through the current block as W
   * and emit an identity matrix for Psi. */
  else if (i1 - i0 < facSpec->minNumRows) {
    BfMatIdentity *Psi = bfMatIdentityNew();
    HANDLE_ERROR();

    bfMatIdentityInit(Psi, i1 - i0);
    HANDLE_ERROR();

    BfMat *W = bfMatCopy(block);
    HANDLE_ERROR();

    *PsiPtr = bfMatIdentityToMat(Psi);
    *WPtr = W;

    success = true;
  }

  /* In the normal mode of operation, we compute a truncated SVD of
   * the current block. We emit Psi := U and W := S*V'. */
  else
    success = getPsiAndW_normal(block, facSpec->tol, PsiPtr, WPtr);

  END_ERROR_HANDLING() {
    success = false;

    /* TODO: make sure to free Psi and W if I allocated them... */
    BF_ASSERT(false);
  }

  bfMatDelete(&block);

  return success;
}

bool getLowRankApproximation(BfFacSpec const *facSpec, BfMat const *PsiStarSubblock,
                             BfMat **PsiSubblockPtr, BfMat **W0SubblockPtr) {
  BF_ERROR_BEGIN();

  BfMat *U = NULL;
  BfMatDiagReal *S = NULL;
  BfMat *VT = NULL;

  BfTruncSpec truncSpec = {
    .usingTol = true,
    .tol = facSpec->tol
  };

  BfBackend backend = BF_BACKEND_LAPACK;

  bool truncated = bfGetTruncatedSvd(PsiStarSubblock, &U, &S, &VT, &truncSpec, backend);
  HANDLE_ERROR();

  BfMat *W0Subblock = VT;

  BfVec *s = bfMatDiagRealGetVecView(S);
  HANDLE_ERROR();

  bfMatScaleRows(W0Subblock, s);
  HANDLE_ERROR();

  BfSizeArray *nonzeroColumnRanges = bfMatGetNonzeroColumnRanges(PsiStarSubblock);
  HANDLE_ERROR();

  BfSize numNonzeroColumnRanges = bfSizeArrayGetSize(nonzeroColumnRanges)/2;

  bool shouldFixSparsity = numNonzeroColumnRanges > 1 ||
    bfSizeArrayGet(nonzeroColumnRanges, 0) > 0 ||
    bfSizeArrayGet(nonzeroColumnRanges, 1) < bfMatGetNumCols(PsiStarSubblock);

  if (shouldFixSparsity) {
    // puts("WARN: make sure to check that remaining columns are actually zero!");

    BfPtrArray indexedBlocks;
    bfInitPtrArrayWithDefaultCapacity(&indexedBlocks);
    HANDLE_ERROR();

    for (BfSize k = 0; k < numNonzeroColumnRanges; ++k) {
      BfSize j0 = bfSizeArrayGet(nonzeroColumnRanges, 2*k);
      BfSize j1 = bfSizeArrayGet(nonzeroColumnRanges, 2*k + 1);

      BfMat *block = bfMatGetColRangeCopy(W0Subblock, j0, j1);
      HANDLE_ERROR();

      BfIndexedMat *indexedBlock = bfIndexedMatNewFromMat(
        0, j0, block, BF_POLICY_STEAL);
      HANDLE_ERROR();

      bfPtrArrayAppend(&indexedBlocks, indexedBlock);
    }

    BfMatBlockCoo *W0SubblockSparsified = bfMatBlockCooNewFromIndexedBlocks(
      bfMatGetNumRows(W0Subblock),
      bfMatGetNumCols(W0Subblock),
      &indexedBlocks,
      BF_POLICY_STEAL);
    HANDLE_ERROR();

    bfMatDelete(&W0Subblock);

    W0Subblock = bfMatBlockCooToMat(W0SubblockSparsified);

    for (BfSize k = 0; k < bfPtrArraySize(&indexedBlocks); ++k) {
      BfIndexedMat *indexedBlock = bfPtrArrayGet(&indexedBlocks, k);
      bfIndexedMatDelete(&indexedBlock);
    }
    bfPtrArrayDeinit(&indexedBlocks);
  }

  END_ERROR_HANDLING() {
    BF_ASSERT(false);
  }

  *PsiSubblockPtr = U;
  *W0SubblockPtr = W0Subblock;

  bfVecDelete(&s);
  bfMatDiagRealDeinitAndDealloc(&S);

  return truncated;
}

void findEpsilonRankCutAndGetNewBlocks(BfFacSpec const *facSpec,
                                       BfTreeNode const *rootRowNode,
                                       BfMat const *PsiStarBlock,
                                       BfConstNodeArray **epsRankCutPtr,
                                       BfMat **PsiBlockPtr,
                                       BfMat **W0BlockPtr) {
  BF_ERROR_BEGIN();

  BfMatBlockDiag *PsiBlock = NULL;
  BfMatBlockDense *W0Block = NULL;

  if (epsRankCutPtr == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (PsiBlockPtr == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (W0BlockPtr == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (bfMatGetNumRows(PsiStarBlock) != bfTreeNodeGetNumPoints(rootRowNode))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* The row nodes in the epsilon-rank cut. */
  BfConstNodeArray *epsRankCut = bfConstNodeArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  /* Array used to accumulate the diagonal subblocks of `PsiBlock`. */
  BfPtrArray PsiSubblocks;
  bfInitPtrArrayWithDefaultCapacity(&PsiSubblocks);
  HANDLE_ERROR();

  /* Array used to accumulate column subblocks of `W0Block`. */
  BfPtrArray W0Subblocks;
  bfInitPtrArrayWithDefaultCapacity(&W0Subblocks);
  HANDLE_ERROR();

  /* Stack used to find the epsilon-rank cut. */
  BfConstPtrArray stack;
  bfConstPtrArrayInitWithDefaultCapacity(&stack);
  HANDLE_ERROR();

  /* Compute index range of the root row node for use below. */
  BfSize i0 = bfTreeNodeGetFirstIndex(rootRowNode);
  BfSize i1 = bfTreeNodeGetLastIndex(rootRowNode);

  /* Push the root row tree node onto the stack to start. */
  bfConstPtrArrayAppend(&stack, rootRowNode);
  HANDLE_ERROR();

  while (!bfConstPtrArrayIsEmpty(&stack)) {
    bool success = true;

    BfTreeNode const *rowNode = bfConstPtrArrayPopLast(&stack);

    /* Get index range of current row node. */
    BfSize i0_ = bfTreeNodeGetFirstIndex(rowNode);
    BfSize i1_ = bfTreeNodeGetLastIndex(rowNode);
    BF_ASSERT(i0 <= i0_ && i0_ < i1_ && i1_ <= i1);

    /* Index range for the rows to copy: */
    BfSize i0__ = i0_ - i0;
    BfSize i1__ = i1_ - i0;
    BF_ASSERT(i0__ < i1__ && i1__ <= bfMatGetNumRows(PsiStarBlock));

    BfMat *PsiStarSubblock = bfMatGetRowRangeCopy(PsiStarBlock, i0__, i1__);
    HANDLE_ERROR();

    BfSize numRowsPsiStarSubblock = bfMatGetNumRows(PsiStarSubblock);
    BfSize numColsPsiStarSubblock = bfMatGetNumCols(PsiStarSubblock);

    BfMat *PsiSubblock = NULL;
    BfMat *W0Subblock = NULL;

    /* If the new Psi subblock has too few rows, we omit an identity
     * matrix for the Psi subblock and pass Psi through as the W0
     * subblock. */
    if (numRowsPsiStarSubblock < facSpec->minNumRows) {
      BfMatIdentity *identity = bfMatIdentityNew();
      HANDLE_ERROR();

      bfMatIdentityInit(identity, numRowsPsiStarSubblock);
      HANDLE_ERROR();

      PsiSubblock = bfMatIdentityToMat(identity);

      W0Subblock = PsiStarSubblock;

      goto next;
    }

    /* If the new Psi subblock has too few columns, we omit an
     * identity matrix for the W0 subblock, and pass the diagonal Psi
     * subblock through as-is.
     *
     * TODO: can we do this before slicing above to save a bit of time? */
    if (numColsPsiStarSubblock < facSpec->minNumCols) {
      PsiSubblock = PsiStarSubblock;

      BfMatIdentity *identity = bfMatIdentityNew();
      HANDLE_ERROR();

      bfMatIdentityInit(identity, numColsPsiStarSubblock);
      HANDLE_ERROR();

      W0Subblock = bfMatIdentityToMat(identity);

      goto next;
    }

    bool truncated = getLowRankApproximation(
      facSpec, PsiStarSubblock, &PsiSubblock, &W0Subblock);
    HANDLE_ERROR();

    bool compressed = bfMatNumBytes(W0Subblock) < bfMatNumBytes(PsiStarSubblock);

    success = truncated && compressed;

    /* If we failed to compute a truncated SVD (i.e., we didn't drop
     * any terms), we push `rowNode`'s children onto the stack, moving
     * deeper into the row tree. */
    if (!success) {
      /* Push nodes onto the stack in the correct order to make sure
       * we continue traversing rows from top to bottom. */
      for (BfSize k = 0; k < rowNode->maxNumChildren; ++k) {
        BfTreeNode const *child = rowNode->child[rowNode->maxNumChildren - k - 1];
        if (child != NULL)
          bfConstPtrArrayAppend(&stack, child);
      }
      continue;
    }

  next:
    /* Append the current row node to the epsilon-rank cut. */
    bfConstNodeArrayAppend(epsRankCut, rowNode);
    HANDLE_ERROR();

    /* Append the Psi subblock. */
    bfPtrArrayAppend(&PsiSubblocks, PsiSubblock);
    HANDLE_ERROR();

    /* Append the W0 subblock. */
    bfPtrArrayAppend(&W0Subblocks, W0Subblock);
    HANDLE_ERROR();
  }

  /* Diagonally concatenate the Psi subblocks. */
  PsiBlock = bfMatBlockDiagNewFromBlocks(&PsiSubblocks, BF_POLICY_STEAL);
  HANDLE_ERROR();

  /* Vertically concatenate the W0 subblocks. */
  W0Block = bfMatBlockDenseNewColFromBlocks(&W0Subblocks, BF_POLICY_STEAL);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    bfConstNodeArrayDeinitAndDealloc(&epsRankCut);

    // TODO: free blocks
    BF_ASSERT(false);
  }

  bfPtrArrayDeinit(&PsiSubblocks);
  bfPtrArrayDeinit(&W0Subblocks);
  bfConstPtrArrayDeinit(&stack);

  *epsRankCutPtr = epsRankCut;
  *PsiBlockPtr = bfMatBlockDiagToMat(PsiBlock);
  *W0BlockPtr = bfMatBlockDenseToMat(W0Block);
}

/* Used below in `mergeAndSplit`.
 *
 * The column node of a merged BF is the parent of each of the merged
 * BFs. This function gets that node from an array of BFs, and
 * additionally checks that every BF has the same parent. */
static BfTreeNode const *getCommonParent(BfPtrArray const *facs) {
  BF_ERROR_BEGIN();

  BfTreeNode const *commonParent = NULL;

  if (bfPtrArrayIsEmpty(facs))
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  for (BfSize k = 0; k < bfPtrArraySize(facs); ++k) {
    BfFac const *fac = bfPtrArrayGet(facs, k);
    BfTreeNode const *parent = bfTreeNodeGetParentConst(fac->colNode);
    if (commonParent == NULL)
      commonParent = parent;
    if (commonParent != parent)
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);
  }

  END_ERROR_HANDLING() {
    BF_DIE();
  }

  return commonParent;
}

BfFac *mergeAndSplit(BfPtrArray const *facs, BfFacSpec const *facSpec) {
  static int CALL_NUMBER = 0;

  BF_ERROR_BEGIN();

  BfFac *mergedFac = NULL;

  BfMatBlockDiag *Psi = NULL;
  BfMatBlockDiag *W0 = NULL;
  BfMatBlockDense *W1 = NULL;

  BF_ASSERT(!bfPtrArrayIsEmpty(facs));

  BfTreeNode const *colNode = getCommonParent(facs);
  HANDLE_ERROR();

  /* TODO: For now, we assume that all current partial factorizations
   * have the same number of W factors. This can be relaxed by filling
   * in with identity matrices below when we concatenate together the
   * trailing W factors. */
  BfSize maxNumW;
  { BfFac *fac = bfPtrArrayGet(facs, 0);
    maxNumW = fac->numW;
    for (BfSize k = 1; k < bfPtrArraySize(facs); ++k) {
      BfFac *otherFac = bfPtrArrayGet(facs, k);
      if (fac->numW != otherFac->numW)
        RAISE_ERROR(BF_ERROR_RUNTIME_ERROR); } }

  // compute merge cut...
  BfConstNodeArray mergeCut = getMergeCut(facs);
  HANDLE_ERROR();

  /* Array used to accumulate row nodes from each of the epsilon rank
   * cuts found when merging the current partial factorizations
   * starting from the merge cut. The resulting span of nodes should
   * be beneath the merge cut and compatible with it. */
  BfConstNodeArray rowNodes;
  bfConstNodeArrayInitWithDefaultCapacity(&rowNodes);
  HANDLE_ERROR();

  /* Array used to accumulate diagonal Psi blocks obtained by merging
   * and splitting all of the current partial factorizations using the
   * merge cut. */
  BfPtrArray PsiBlocks;
  bfInitPtrArrayWithDefaultCapacity(&PsiBlocks);
  HANDLE_ERROR();

  BfPtrArray W0Blocks;
  bfInitPtrArrayWithDefaultCapacity(&W0Blocks);
  HANDLE_ERROR();

  BfPtrArray W1Blocks;
  bfInitPtrArrayWithDefaultCapacity(&W1Blocks);
  HANDLE_ERROR();

  /* Iterate over each of the row tree nodes in the merge cut, and do
   * the following:
   * - extract merged Psi blocks
   * - correctly sift the blocks from leading W factors
   * - find the epsilon rank cuts of the Psi blocks
   * - emit the new Psi and W blocks
   *
   * Note that we take care to ensure that we take full advantage of
   * the sparsity of the new W blocks here to compress the
   * factorization as much as possible. */
  for (BfSize k = 0; k < bfConstNodeArraySize(&mergeCut); ++k) {
    BfTreeNode const *rowNode = bfConstNodeArrayGet(&mergeCut, k);

    BfConstNodeArray *epsRankCut = NULL;
    BfMat *PsiStarBlock = NULL;
    BfMat *PsiBlock = NULL;
    BfMat *W0Block = NULL;
    BfMat *W1Block = NULL;

    /* Get current row block of Psi. We simultaneously get the column
     * index ranges corresponding to each column subblock. */
    getPsiAndW0BlocksByRowNode(facs, rowNode, &PsiStarBlock, &W1Block);

    bfPtrArrayAppend(&W1Blocks, W1Block);
    HANDLE_ERROR();

    /* Find the epsilon-rank cut and compute the correspond blocks of
     * the new Psi and W0 factors while simultaneously sifting the W0
     * blocks from the current partial factorizations together to get
     * the next block of the W1 factor. */
    findEpsilonRankCutAndGetNewBlocks(
      facSpec,
      rowNode,
      PsiStarBlock,
      &epsRankCut,
      &PsiBlock,
      &W0Block);
    HANDLE_ERROR();

    bfConstNodeArrayExtend(&rowNodes, epsRankCut);
    HANDLE_ERROR();

    bfPtrArrayAppend(&PsiBlocks, PsiBlock);
    HANDLE_ERROR();

    bfPtrArrayAppend(&W0Blocks, W0Block);
    HANDLE_ERROR();

    bfConstNodeArrayDeinitAndDealloc(&epsRankCut);
  }

  /* Diagonally concatenate the new Psi blocks to get the leading
   * factor of the merged butterfly factorization. */
  Psi = bfMatBlockDiagNewFromBlocks(&PsiBlocks, BF_POLICY_STEAL);
  HANDLE_ERROR();

  /* Diagonally concatenate together the new W0 blocks. */
  W0 = bfMatBlockDiagNewFromBlocks(&W0Blocks, BF_POLICY_STEAL);
  HANDLE_ERROR();

  /* Vertically concatenate the collected W1 blocks. */
  W1 = bfMatBlockDenseNewColFromBlocks(&W1Blocks, BF_POLICY_COPY);
  HANDLE_ERROR();

  /* Create and initialize the merged factorization. */
  // TODO: should really break this out into a separate function once
  // we refactor PartialFac into an honest ButterflyFactorization
  {
    mergedFac = bfMemAlloc(1, sizeof(BfFac));
    if (mergedFac == NULL)
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

    mergedFac->colNode = colNode;

    mergedFac->rowNodes = rowNodes;

    mergedFac->Psi = bfMatBlockDiagToMat(Psi);

    mergedFac->numW = maxNumW + 1;
    BF_ASSERT(mergedFac->numW >= 2);

    mergedFac->W = bfMemAllocAndZero(mergedFac->numW, sizeof(BfMat *));
    if (mergedFac->W == NULL)
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

    setPartialFacWBlock(mergedFac, 0, bfMatBlockDiagToMat(W0));
    setPartialFacWBlock(mergedFac, 1, bfMatBlockDenseToMat(W1));

    /* Block together trailing W factors from the current partial
     * butterfly factorizations to get the trailing W factors for the
     * newly merged butterfly factorization */
    for (BfSize k = 1; k < maxNumW; ++k) {
      BfPtrArray WkBlocks;
      bfInitPtrArrayWithDefaultCapacity(&WkBlocks);
      HANDLE_ERROR();

      for (BfSize l = 0; l < bfPtrArraySize(facs); ++l) {
        BfFac const *fac = bfPtrArrayGet(facs, l);

        BfMat *WkBlock = bfMatGet(fac->W[k], BF_POLICY_COPY);
        BF_ASSERT(WkBlock != NULL);

        bfPtrArrayAppend(&WkBlocks, WkBlock);
        HANDLE_ERROR();
      }

      BfMatBlockDiag *Wk = bfMatBlockDiagNewFromBlocks(&WkBlocks, BF_POLICY_STEAL);
      HANDLE_ERROR();

      setPartialFacWBlock(mergedFac, k + 1, bfMatBlockDiagToMat(Wk));

      for (BfSize l = 0; l < bfPtrArraySize(&WkBlocks); ++l) {
        BfMat *WkBlock = bfPtrArrayGet(&WkBlocks, l);
        bfMatDelete(&WkBlock);
      }
      bfPtrArrayDeinit(&WkBlocks);
    }

    /* Make sure everything was set successfully: */
#if BF_DEBUG
    BF_ASSERT(mergedFac->Psi != NULL);
    for (BfSize i = 0; i < bfConstNodeArraySize(&mergedFac->rowNodes); ++i)
      BF_ASSERT(bfConstNodeArrayGet(&mergedFac->rowNodes, i) != NULL);
    BF_ASSERT(mergedFac->colNode != NULL);
    for (BfSize k = 0; k < mergedFac->numW; ++k)
      BF_ASSERT(mergedFac->W[k] != NULL);
#endif
  }

  END_ERROR_HANDLING() {
    BF_ASSERT(false);
  }

  bfPtrArrayDeinit(&PsiBlocks);
  bfPtrArrayDeinit(&W0Blocks);
  bfPtrArrayDeinit(&W1Blocks);

  ++CALL_NUMBER;

  return mergedFac;
}
