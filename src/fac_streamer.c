#include <bf/fac_streamer.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/linalg.h>
#include <bf/mat_block_diag.h>
#include <bf/mat_diag_real.h>
#include <bf/mat_identity.h>
#include <bf/ptr_array.h>
#include <bf/tree.h>
#include <bf/tree_node.h>
#include <bf/tree_iter_post_order.h>
#include <bf/vec_real.h>

#undef BF_ARRAY_DEFAULT_CAPACITY
#define BF_ARRAY_DEFAULT_CAPACITY 32

struct BfFacStreamer {
  BfSize rowTreeInitDepth;
  BfSize colTreeInitDepth;

  BfTree *rowTree;
  BfTree *colTree;

  BfTreeIter *colTreeIter;

  BfReal tol;
  BfSize minNumCols;

  BfTreeNode *currentColNode;
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

void bfFacStreamerInit(BfFacStreamer *facStreamer, BfFacStreamerSpec const *spec) {
  BEGIN_ERROR_HANDLING();

  facStreamer->rowTree = spec->rowTree;
  facStreamer->colTree = spec->colTree;

  facStreamer->rowTreeInitDepth = spec->rowTreeInitDepth;
  facStreamer->colTreeInitDepth = spec->colTreeInitDepth;

  facStreamer->tol = spec->tol;
  facStreamer->minNumCols = spec->minNumCols;

  BfTreeIterPostOrder *iter = bfTreeIterPostOrderNew();
  HANDLE_ERROR();

  bfTreeIterPostOrderInit(iter, facStreamer->colTree);
  HANDLE_ERROR();

  facStreamer->colTreeIter = bfTreeIterPostOrderToTreeIter(iter);

  facStreamer->currentColNode = bfTreeIterGetCurrentNode(facStreamer->colTreeIter);
  if (facStreamer->currentColNode == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  END_ERROR_HANDLING() {
    bfFacStreamerDeinit(facStreamer);
  }
}

void bfFacStreamerDeinit(BfFacStreamer *facStreamer) {
  (void)facStreamer;
  assert(false);
}

static bool getPsiAndW_skinny(BfMat const *block, BfMat **PsiPtr, BfMat **WPtr) {
  BEGIN_ERROR_HANDLING();

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
  BEGIN_ERROR_HANDLING();

  bool success = true;

  BfSize m = bfMatGetNumRows(block);
  BfSize n = bfMatGetNumCols(block);

  BfSize kmax = m < n ? m : n;

  /* Compute truncated SVD of the block. */
  BfMat *Psi = NULL;
  BfMatDiagReal *S = NULL;
  BfMat *W = NULL;
  BfTruncSpec truncSpec = {.usingTol = true, .tol = tol};
  bfGetTruncatedSvd(block, &Psi, &S, &W, truncSpec, BF_BACKEND_LAPACK);
  HANDLE_ERROR();

  BfSize k = bfMatGetNumCols(Psi);

  /* Scale W if we successfully compressed this block. */
  success = k < kmax;
  if (success)
    bfMatMulInplace(W, bfMatDiagRealToMat(S));

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

static bool getPsiAndW(BfMat const *mat, BfTreeNode const *rowNode,
                       BfTreeNode const *colNode,
                       BfReal tol, BfSize minNumCols,
                       BfMat **PsiPtr, BfMat **WPtr) {
  BEGIN_ERROR_HANDLING();

  BfSize i0 = bfTreeNodeGetFirstIndex(rowNode);
  BfSize i1 = bfTreeNodeGetLastIndex(rowNode);

  /* TODO: should be something like: `bfMatGetRowRangeViewConst` */
  BfMat *block = bfMatGetRowRange((BfMat *)mat, i0, i1);
  HANDLE_ERROR();

  /* If there are too few columns in the current block, we don't
   * bother trying to compress. Instead, we just pass the current
   * block through for the Psi block and emit an identity matrix for
   * the new W block. Otherwise, we compute a truncated SVD of the
   * current block. We emit Psi := U and W := S*V'. */
  bool success = bfMatGetNumCols(mat) < minNumCols ?
    getPsiAndW_skinny(block, PsiPtr, WPtr) :
    getPsiAndW_normal(block, tol, PsiPtr, WPtr);

  END_ERROR_HANDLING() {
    success = false;

    /* TODO: make sure to free Psi and W if I allocated them... */
    assert(false);
  }

  bfMatDelete(&block);

  return success;
}

static void continueFactorizing(BfFacStreamer *facStreamer) {
  BfTreeNode const *currentNode = NULL;

  while (!bfTreeIterIsDone(facStreamer->colTreeIter)) {
    bfTreeIterNext(facStreamer->colTreeIter);

    currentNode = bfTreeIterGetCurrentNode(facStreamer->colTreeIter);
    if (bfTreeNodeIsLeaf(currentNode))
      break;

    BfSize currentDepth = currentNode->depth;
    BfSize numChildren = currentNode->maxNumChildren;

    (void)currentDepth;
    (void)numChildren;

    assert(false);

    // get Psi blocks for children

    // merge

    // split

    // emit

  }

  facStreamer->currentColNode = bfTreeIterGetCurrentNode(facStreamer->colTreeIter);
}

/* Notes:
 * - The rows and columns of `mat` should already be permuted into the
 *   orders defined by `facStreamer->rowTree` and
 *   `facStreamer->colTree`. */
void bfFacStreamerFeed(BfFacStreamer *facStreamer, BfMat const *Phi, BfVecReal const *Lam) {
  BEGIN_ERROR_HANDLING();

  /* Get the current column tree leaf node */
  BfTreeNode const *colNode = facStreamer->currentColNode;
  assert(bfTreeNodeIsLeaf(colNode));

  /* Make sure `mat` has the right number of rows */
  if (bfMatGetNumRows(Phi) != bfTreeGetNumPoints(facStreamer->rowTree))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* Make sure `mat` has the right number of columns */
  if (bfMatGetNumCols(Phi) != bfTreeNodeGetNumPoints(colNode))
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

  /* Reverse the contents of the stack to make sure we traverse it
   * from the top of the column block down. */
  bfPtrArrayReverse(&stack);

  while (!bfPtrArrayIsEmpty(&stack)) {
    BfTreeNode const *rowNode = bfPtrArrayPopLast(&stack);

    BfMat *Psi = NULL, *W = NULL;
    bool metTol = getPsiAndW(
      Phi, rowNode, colNode, facStreamer->tol, facStreamer->minNumCols, &Psi, &W);
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

bool bfFacStreamerIsDone(BfFacStreamer const *facStreamer) {
  return bfTreeIterIsDone(facStreamer->colTreeIter);
}

BfMat *bfFacStreamerGetFac(BfFacStreamer const *facStreamer) {
  (void)facStreamer;

  assert(false);
}

BfTreeNode *bfFacStreamerGetCurrentColumnNode(BfFacStreamer const *facStreamer) {
  return facStreamer->currentColNode;
}
