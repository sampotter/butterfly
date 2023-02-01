#include <bf/fac_streamer.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/linalg.h>
#include <bf/mat_block_diag.h>
#include <bf/mat_diag_real.h>
#include <bf/ptr_array.h>
#include <bf/tree.h>
#include <bf/tree_node.h>
#include <bf/tree_iter.h>
#include <bf/vec_real.h>

#undef BF_ARRAY_DEFAULT_CAPACITY
#define BF_ARRAY_DEFAULT_CAPACITY 32

struct BfFacStreamer {
  BfReal tol;
  BfSize rowTreeInitDepth;
  BfSize colTreeInitDepth;

  BfTree const *rowTree;
  BfTree const *colTree;

  BfTreeIter *colTreeIter;
  BfTreeNode const *currentColNode;
};

BfFacStreamer *bfFacStreamerNew() {
  BEGIN_ERROR_HANDLING();

  BfFacStreamer *facStreamer = malloc(sizeof(BfFacStreamer));
  if (facStreamer == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING()
    facStreamer = NULL;

  return facStreamer;
}

void bfFacStreamerInit(BfFacStreamer *facStreamer, BfTree const *rowTree,
                       BfTree const *colTree, BfSize rowTreeInitDepth,
                       BfSize colTreeInitDepth, BfReal tol) {
  facStreamer->tol = tol;

  facStreamer->rowTreeInitDepth = rowTreeInitDepth;
  facStreamer->colTreeInitDepth = colTreeInitDepth;

  facStreamer->rowTree = rowTree;
  facStreamer->colTree = colTree;

  facStreamer->colTreeIter = NULL;
  facStreamer->currentColNode = NULL;
}

static bool getPsiAndW(BfMat const *mat, BfTreeNode const *rowNode,
                       BfTreeNode const *colNode, BfReal tol,
                       BfMat **Psi, BfMat **W) {
  (void)colNode;

  BEGIN_ERROR_HANDLING();

  BfSize i0 = bfTreeNodeGetFirstIndex(rowNode);
  BfSize i1 = bfTreeNodeGetLastIndex(rowNode);

  BfSize m = i1 - i0;
  BfSize n = bfMatGetNumCols(mat);
  BfSize kmax = m < n ? m : n;

  /* TODO: should be something like: `bfMatGetRowRangeViewConst` */
  BfMat *block = bfMatGetRowRange((BfMat *)mat, i0, i1);
  HANDLE_ERROR();

  /* Compute truncated SVD of the block. */
  BfMatDiagReal *S;
  BfTruncSpec truncSpec = {.usingTol = true, .tol = tol};
  bfGetTruncatedSvd(block, Psi, &S, W, truncSpec, BF_BACKEND_LAPACK);
  HANDLE_ERROR();

  /* Bail if we didn't compress this block */
  BfSize k = bfMatGetNumCols(*Psi);
  if (k == kmax) {
    bfMatDelete(Psi);
    bfMatDelete(W);
    bfMatDiagRealDeinitAndDealloc(&S);
    return false;
  }

  BfVec *Svec = bfMatDiagRealGetVecView(S);
  bfMatScaleRows(*W, Svec);

  END_ERROR_HANDLING() {}

  bfVecDelete(&Svec);
  bfMatDiagRealDeinitAndDealloc(&S);
  bfMatDelete(&block);

  return true;
}

static void continueFactorizing(BfFacStreamer *facStreamer) {
  BfTreeNode const *currentNode = NULL;

  while (!bfTreeIterIsDone(facStreamer->colTreeIter)) {
    currentNode = facStreamer->colTreeIter->currentNode;
    if (bfTreeNodeIsLeaf(currentNode))
      break;

    BfSize currentDepth = currentNode->depth;
    BfSize numChildren = currentNode->maxNumChildren;


    (void)currentDepth;
    (void)numChildren;

    // get Psi blocks for children


    // merge

    // split

    // emit

  }

  facStreamer->currentColNode = facStreamer->colTreeIter->currentNode;
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
  if (bfMatGetNumRows(mat) != bfTreeGetNumPoints(facStreamer->rowTree))
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
  BfPtrArray stack = bfTreeGetLevelPtrArray(
    facStreamer->rowTree, facStreamer->rowTreeInitDepth);
  HANDLE_ERROR();

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
    for (BfSize i = rowNode->maxNumChildren; i > 0; --i)
      bfPtrArrayAppend(&stack, rowNode->child[i - 1]);
  }

  continueFactorizing(facStreamer);

  END_ERROR_HANDLING() {}

  bfPtrArrayDeinit(&W_blocks);
  bfPtrArrayDeinit(&Psi_blocks);
  bfPtrArrayDeinit(&stack);
}

bool bfFacStreamerDone(BfFacStreamer const *facStreamer) {
  return bfTreeIterIsDone(facStreamer->colTreeIter);
}

BfMat *bfFacStreamerGetFac(BfFacStreamer const *facStreamer) {
  (void)facStreamer;

  assert(false);
}
