#include "bf_streaming.h"

#include <assert.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_block_diag.h>
#include <bf/ptr_array.h>
#include <bf/tree.h>

#undef BF_ARRAY_DEFAULT_CAPACITY
#define BF_ARRAY_DEFAULT_CAPACITY 32

struct BfFacStreamer {
  BfReal tol;
  bool adaptive;
  BfSize rowTreeInitDepth;

  BfTree const *rowTree;
  BfTree const *colTree;

  BfTreeIter colTreeIter;
  BfTreeNode const *currentColNode;

}

static bool getPsiAndW(BfMat const *mat, BfTreeNode const *rowNode,
                       BfTreeNode const *colNode, BfReal tol,
                       BfMat **Psi, BfMat **W) {
  BEGIN_ERROR_HANDLING();

  BfSize i0 = bfTreeNodeGetFirstIndex(rowNode);
  BfSize i1 = bfTreeNodeGetLastIndex(rowNode);

  BfSize m = i1 - i0;
  BfSize n = bfMatGetNumCols(mat);
  BfSize kmax = m < n ? m : n;

  /* TODO: should be something like: `bfMatGetRowRangeViewConst` */
  BfMat const *block = bfMatGetRowRange((BfMat *)mat, i0, i1);
  HANDLE_ERROR();

  /* Compute truncated SVD of the block. */
  BfMatDiagReal *S;
  bfMatGetTruncatedSvd(block, *Psi, S, *W, tol);
  HANDLE_ERROR();

  /* Bail if we didn't compress this block */
  BfSize k = bfMatGetNumCols(U);
  if (k == kmax) {
    bfMatDelete(Psi);
    bfMatDelete(W);
    bfMatDiagRealDeallocAndDelete(&S);
    return false;
  }

  BfVec *Svec = bfVecRealToVec(bfMatDiagRealGetVecView(S));
  bfMatScaleRows(*W, Svec);

  END_ERROR_HANDLING() {}

  bfVecDelete(&Svec);
  bfMatDiagRealDeallocAndDelete(&S);
  bfMatDelete(&block);

  return true;
}

static void continueFactorizing(BfFacStreamer *facStreamer) {
  BEGIN_ERROR_HANDLING();

  BfTreeNode const *currentNode = NULL;

  while (!bfTreeLevelIterIsDone(facStreamer->colTreeIter)
         && !bfTreeNodeIsLeaf(facStreamer->colTreeIter->currentNode)) {
    currentNode = facStreamer->colTreeIter->currentNode;

    BfSize currentDepth = bfTreeNodeGetDepth(currentNode);
    BfSize numChildren = bfTreeNodeGetNumChildren(currentNode);

    // get Psi blocks for children


    // merge

    // split

    // emit

  }

  facStreamer->currentColNode = facStreamer->colTreeIter->currentNode;

  END_ERROR_HANDLING() {}
}

/* Notes:
 * - The rows and columns of `mat` should already be permuted into the
 *   orders defined by `facStreamer->rowTree` and
 *   `facStreamer->colTree`. */
void bfFacStreamerFeed(BfFacStreamer *facStreamer, BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  /* Get the current column tree leaf node */
  BfTreeNode const *colNode = facStreamer->currentColNode;
  assert(bfTreeNodeIsLeaf(colNode));

  /* Make sure `mat` has the right number of rows */
  if (bfMatGetNumRows(mat) != bfTreeGetNumPoints(facStreamer->rowTree));
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* Make sure `mat` has the right number of columns */
  if (bfMatGetNumCols(mat) != bfTreeNodeGetNumPoints(colNode))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* Array used to accumulate blocks for block diagonal Psi factor */
  BfPtrArray Psi_blocks;
  bfInitPtrArray(&Psi_blocks, BF_ARRAY_DEFAULT_CAPACITY);
  HANDLE_ERROR();

  /* Array for blocks making up the W factor */
  BfPtrArray W_blocks;
  bfInitPtrArray(&W_blocks, BF_ARRAY_DEFAULT_CAPACITY);
  HANDLE_ERROR();

  /* Create the stack used to adaptively find the cut through the row
   * tree to start the adaptive butterfly factorization */
  BfPtrArray stack;
  bfInitPtrArray(&stack, BF_ARRAY_DEFAULT_CAPACITY);
  HANDLE_ERROR();

  /* Fill it with nodes on the starting level of the row tree
   *
   * TODO: simplify this by adding `bfTreeGetLevelNodes` function */
  { BfTreeLevelIter *iter = bfTreeGetLevelIter(
      facStreamer->rowTree, BF_TREE_TRAVERSAL_LR_LEVEL_ORDER,
      facStreamer->rowTreeInitDepth);
    for (BfSize i = 0; i < bfPtrArraySize(&iter->levelNodes); ++i)
      bfPtrArrayAppend(&stack, bfPtrArrayGet(&iter->levelNodes, i));
    bfTreeLevelIterFree(iter); }

  while (!bfPtrArrayIsEmpty(&stack)) {
    BfTreeNode const *rowNode = bfPtrArrayPopLast(&stack);

    BfMat *Psi = NULL, *W = NULL;
    bool metTol = getPsiAndW(mat, rowNode, colNode, facStreamer->tol, &Psi, &W);
    HANDLE_ERROR();

    /* Accumulate the Psi and W blocks and continue if we successfully
     * compressed the current block */
    if (metTol) {
      bfPtrArrayAppend(&Psi_blocks, Psi);
      bfPtrArrayAppend(&W_blocks, W);
      continue;
    }

    /* Push the children of the current row node onto the stack in
     * reverse order so that we traverse `mat` top to bottom */
    for (BfSize i = bfTreeNodeGetNumChildren(rowNode) - 1; i >= 0; --i)
      bfPtrArrayAppend(bfTreeNodeGetChild(rowNode, i));
  }

  continueFactorizing(facStreamer);

  END_ERROR_HANDLING() {}

  bfPtrArrayDeinit(&W_blocks);
  bfPtrArrayDeinit(&Psi_blocks);
  bfPtrArrayDeinit(&stack);
}

bool bfFacStreamerDone(BfFacStreamer const *facStreamer) {
  return bfTreeLevelIterIsDone(facStreamer->colTreeIter);
}
