#include <bf/fac_streamer.h>

#if BF_DEBUG
#include <stdio.h>
#endif

#include <bf/assert.h>
#include <bf/const.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/fac_span.h>
#include <bf/indexed_mat.h>
#include <bf/linalg.h>
#include <bf/logging.h>
#include <bf/mat_block_coo.h>
#include <bf/mat_block_dense.h>
#include <bf/mat_block_diag.h>
#include <bf/mat_dense_real.h>
#include <bf/mat_diag_real.h>
#include <bf/mat_identity.h>
#include <bf/mat_product.h>
#include <bf/mem.h>
#include <bf/ptr_array.h>
#include <bf/size_array.h>
#include <bf/tree.h>
#include <bf/tree_node.h>
#include <bf/tree_iter_post_order.h>
#include <bf/vec_real.h>

typedef struct {
  BfTreeNode const *treeNode;
  BfMat *mat;
} MatWithTreeNodeKey;

struct BfFacStreamer {
  BfFacSpec const *facSpec;

  BfTreeIter *colTreeIter;

  BfPerm *rowTreeReversePerm;

  /* This is basically an association list which maps from column tree
   * nodes to lists of Psi and W blocks. We use this to track blocks
   * as we traverse the column tree, merging and splitting blocks. */
  // TODO: change this to `BfConstPtrArray *partialFacs;`
  BfPtrArray partialFacs;

  /* An association list mapping from column tree nodes to contiguous
   * blocks of Phi. We use this to check how accurate each partial
   * factorization is as we stream them. */
  BfPtrArray *prevPhis;
};

BfFacStreamer *bfFacStreamerNew() {
  BF_ERROR_BEGIN();

  BfFacStreamer *facStreamer = bfMemAlloc(1, sizeof(BfFacStreamer));
  if (facStreamer == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BF_ERROR_END()
    facStreamer = NULL;

  return facStreamer;
}

void bfFacStreamerInit(BfFacStreamer *facStreamer, BfFacSpec const *facSpec) {
  BF_ERROR_BEGIN();

  facStreamer->facSpec = facSpec;

  facStreamer->rowTreeReversePerm = bfPermGetReversePerm(&facSpec->rowTree->perm);
  HANDLE_ERROR();

  BfTreeIterPostOrder *iter = bfTreeIterPostOrderNew();
  HANDLE_ERROR();

  bfTreeIterPostOrderInit(iter, facSpec->colTree);
  HANDLE_ERROR();

  facStreamer->colTreeIter = bfTreeIterPostOrderToTreeIter(iter);

  facStreamer->partialFacs = bfGetUninitializedPtrArray();
  HANDLE_ERROR();

  bfInitPtrArrayWithDefaultCapacity(&facStreamer->partialFacs);
  HANDLE_ERROR();

  if (facSpec->compareRelativeErrors) {
    facStreamer->prevPhis = bfPtrArrayNewWithDefaultCapacity();
    HANDLE_ERROR();
  } else {
    facStreamer->prevPhis = NULL;
  }

  BF_ERROR_END() {
    bfFacStreamerDeinit(facStreamer);
  }
}

void bfFacStreamerDeinit(BfFacStreamer *facStreamer) {
  bfTreeIterDelete(&facStreamer->colTreeIter);

  bfPermDeinitAndDealloc(&facStreamer->rowTreeReversePerm);

  for (BfSize i = 0; i < bfPtrArraySize(&facStreamer->partialFacs); ++i) {
    BfFac *fac = bfPtrArrayGet(&facStreamer->partialFacs, i);
    bfFacDelete(&fac);
  }
  bfPtrArrayDeinit(&facStreamer->partialFacs);

  if (facStreamer->facSpec->compareRelativeErrors) {
    for (BfSize i = 0; i < bfPtrArraySize(facStreamer->prevPhis); ++i) {
      MatWithTreeNodeKey *entry = bfPtrArrayGet(facStreamer->prevPhis, i);
      BF_ASSERT(!bfMatIsView(entry->mat));
      bfMatDelete(&entry->mat);
      bfMemFree(entry);
    }
    bfPtrArrayDelete(&facStreamer->prevPhis);
  }

  facStreamer->facSpec = NULL;
}

void bfFacStreamerDealloc(BfFacStreamer **facStreamer) {
  bfMemFree(*facStreamer);
  *facStreamer = NULL;
}

void bfFacStreamerDelete(BfFacStreamer **facStreamer) {
  bfFacStreamerDeinit(*facStreamer);
  bfFacStreamerDealloc(facStreamer);
}

BfSize bfFacStreamerGetNumRows(BfFacStreamer const *facStreamer) {
  return bfTreeGetNumPoints(facStreamer->facSpec->rowTree);
}

/* Find the current partial factorizations associated with each of the
 * children of the current column node and return them in PtrArray. */
static BfPtrArray getCurrentPartialFacs(BfFacStreamer const *facStreamer) {
  BF_ERROR_BEGIN();

  BfTreeNode const *currentColNode = bfFacStreamerGetCurrentColumnNode(facStreamer);
  BF_ASSERT(!bfTreeNodeIsLeaf(currentColNode));

  BfPtrArray currentPartialFacs;
  bfInitPtrArrayWithDefaultCapacity(&currentPartialFacs);
  HANDLE_ERROR();

  for (BfSize k = 0; k < currentColNode->maxNumChildren; ++k) {
    BfTreeNode const *childColNode = currentColNode->child[k];
    if (childColNode == NULL)
      continue;

    BfFac *partialFac = NULL;
    for (BfSize l = 0; l < bfPtrArraySize(&facStreamer->partialFacs); ++l) {
      partialFac = bfPtrArrayGet(&facStreamer->partialFacs, l);
      if (partialFac->colNode == childColNode)
        break;
      partialFac = NULL;
    }
    BF_ASSERT(partialFac != NULL);

    bfPtrArrayAppend(&currentPartialFacs, partialFac);
    HANDLE_ERROR();
  }

  BF_ASSERT(!bfPtrArrayIsEmpty(&currentPartialFacs));

  BF_ERROR_END() {
    BF_DIE();
  }

  return currentPartialFacs;
}

static BfMat const *getPrevPhiViewByColNode(BfFacStreamer const *facStreamer,
                                            BfTreeNode const *treeNode) {
  BfMat const *mat = NULL;
  for (BfSize k = 0; k < bfPtrArraySize(facStreamer->prevPhis); ++k) {
    MatWithTreeNodeKey const *entry = bfPtrArrayGet(facStreamer->prevPhis, k);
    if (entry->treeNode == treeNode) {
      mat = bfMatGetView(entry->mat);
      break;
    }
  }
  return mat;
}

/* Have `facStreamer` store `Phi` so that we can check relative errors
 * without having to recompute bands of eigenvectors later. After
 * calling this function, `facStreamer` owns `Phi`. */
static void addPrevPhi(BfFacStreamer *facStreamer, BfTreeNode const *treeNode, BfMat *Phi) {
  MatWithTreeNodeKey *entry = bfMemAlloc(1, sizeof(MatWithTreeNodeKey));
  entry->treeNode = treeNode;
  entry->mat = bfMatGet(Phi, BF_POLICY_STEAL);
  bfPtrArrayAppend(facStreamer->prevPhis, entry);
}

static void deletePrevPhiByColNode(BfFacStreamer *facStreamer, BfTreeNode const *colNode) {
  for (BfSize k = bfPtrArraySize(facStreamer->prevPhis); k > 0; --k) {
    MatWithTreeNodeKey *entry = bfPtrArrayGet(facStreamer->prevPhis, k - 1);
    if (entry->treeNode == colNode) {
      BF_ASSERT(!bfMatIsView(entry->mat));
      bfMatDelete(&entry->mat);
      bfMemFree(entry);
      bfPtrArrayRemove(facStreamer->prevPhis, k - 1);
    }
  }
}

static void addPrevPhiForChildNodesAndDeletePrev(BfFacStreamer *facStreamer, BfTreeNode const *currentColNode) {
  BfSize m = bfFacStreamerGetNumRows(facStreamer);
  BfSize n = bfTreeNodeGetNumPoints(currentColNode);
  BfMatDenseReal *prevPhi = bfMatDenseRealNewWithValue(m, n, BF_NAN);

  /** Concatenate together children of the current colum node: */
  for (BfSize k = 0; k < currentColNode->maxNumChildren; ++k) {
    BfTreeNode const *childColNode = currentColNode->child[k];
    if (childColNode == NULL)
      continue;

    /* The column node offsets give us the offset into the array of
     * points backing the column tree. To get column indices into
     * `Phi`, we need to offset them appropriately. */
    BfSize jOffset = currentColNode->offset[0];
    BfSize j0 = currentColNode->offset[k] - jOffset;
    BfSize j1 = currentColNode->offset[k + 1] - jOffset;

    BfMat const *prevPhiBlock = getPrevPhiViewByColNode(facStreamer, childColNode);
    bfMatDenseRealSetBlock(prevPhi, 0, m, j0, j1, bfMatConstToMatDenseRealConst(prevPhiBlock));
    bfMatDelete((BfMat **)&prevPhiBlock);
  }

  addPrevPhi(facStreamer, currentColNode, bfMatDenseRealToMat(prevPhi));

  /** Delete the dense Phi blocks we just merged together: */
  for (BfSize k = 0; k < currentColNode->maxNumChildren; ++k) {
    BfTreeNode const *childColNode = currentColNode->child[k];
    if (childColNode == NULL)
      continue;
    deletePrevPhiByColNode(facStreamer, childColNode);
  }

  BF_ASSERT(bfMatIsView(bfMatDenseRealToMat(prevPhi)));
  bfMatDenseRealDelete(&prevPhi);
}

static void deletePrevFacs(BfFacStreamer *facStreamer, BfTreeNode const *currentColNode) {
  /* Iterate over each of the children of the current column node: */
  for (BfSize k = 0; k < currentColNode->maxNumChildren; ++k) {
    BfTreeNode const *childColNode = currentColNode->child[k];
    if (childColNode == NULL)
      continue;

    /* Find the partial factorization corresponding to the current
     * child column node: */
    for (BfSize j = bfPtrArraySize(&facStreamer->partialFacs); j > 0; --j) {
      BfFac *fac = bfPtrArrayGet(&facStreamer->partialFacs, j - 1);
      if (fac->colNode != childColNode)
        continue;

      /* Remove it from the array of partial factorizations: */
      bfPtrArrayRemove(&facStreamer->partialFacs, j - 1);

      /* Free the partial factorization for the current column node.
       *
       * NOTE: there are a lot of assumptions resting on the following
       * line being correct. The most recently merged factorization
       * has stolen some of the guts of this factorization. We want to
       * make sure we free everything we need to free here, without
       * accidentally freeing parts of the new factorization! */
      bfFacDelete(&fac);
    }
  }
}

/* Check the current relative error. */
static void checkRelError(BfMat const *Phi, BfFac const *fac) {
  BfSize n = bfMatGetNumCols(Phi);
  if (n == 0) {
    bfLogInfo("- merged fac has no columns---skipping\n");
    return;
  }

  BfVec *x = bfVecRealToVec(bfVecRealNewRandn(n));
  BfVec *y_gt = bfMatMulVec(Phi, x);
  BfVec *y = partialFacMulVec(fac, x);
  BfReal rel_error = bfVecDistMax(y, y_gt)/bfVecNormMax(y_gt);
  bfLogInfo("- rel max error for random MVP: %g\n", rel_error);
  bfVecDelete(&y);
  bfVecDelete(&y_gt);
  bfVecDelete(&x);
}

static void continueFactorizing(BfFacStreamer *facStreamer) {
  BF_ERROR_BEGIN();

  /* Continue the post-order traversal until the next leaf node */
  while (!bfTreeIterIsDone(facStreamer->colTreeIter)) {
    BfTreeNode const *currentColNode = bfTreeIterGetCurrentNode(facStreamer->colTreeIter);
    if (bfTreeNodeIsLeaf(currentColNode))
      break;

    bfLogInfo("merging [%lu, %lu)\n",
              bfTreeNodeGetFirstIndex(currentColNode),
              bfTreeNodeGetLastIndex(currentColNode));

    if (facStreamer->facSpec->compareRelativeErrors)
      addPrevPhiForChildNodesAndDeletePrev(facStreamer, currentColNode);

    BfPtrArray currentPartialFacs = getCurrentPartialFacs(facStreamer);
    HANDLE_ERROR();

    BfFac *mergedFac = mergeAndSplit(
      &currentPartialFacs, facStreamer->facSpec, BF_POLICY_STEAL);
    HANDLE_ERROR();

#if BF_DEBUG
    /* Dump the merged factorization's blocks for plotting: */
    FILE *fp = fopen("Psi.txt", "w");
    bfMatPrintBlocksDeep(mergedFac->Psi, fp, 0, 0, 0);
    fclose(fp);
    for (BfSize k = 0; k < mergedFac->numW; ++k) {
      char path[256];
      sprintf(path, "W%lu.txt", k);
      fp = fopen(path, "w");
      bfMatPrintBlocksDeep(mergedFac->W[k], fp, 0, 0, 0);
      fclose(fp);
    }
#endif

    /* Check the relative error so far: */
    if (facStreamer->facSpec->compareRelativeErrors) {
      BfMat const *Phi = getPrevPhiViewByColNode(facStreamer, currentColNode);
      BF_ASSERT(Phi != NULL);
      checkRelError(Phi, mergedFac);
      bfMatDelete((BfMat **)&Phi);
    }

    bfPtrArrayAppend(&facStreamer->partialFacs, mergedFac);
    HANDLE_ERROR();

    /* Reclaim memory by freeing unused things at this point: */
    deletePrevFacs(facStreamer, currentColNode);

    bfTreeIterNext(facStreamer->colTreeIter);
    HANDLE_ERROR();

    bfPtrArrayDeinit(&currentPartialFacs);
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

static void addPartialFac(BfFacStreamer *facStreamer, BfFac *partialFac) {
  BF_ERROR_BEGIN();

  BF_ASSERT(bfFacStreamerGetNumRows(facStreamer) == partialFacGetNumRows(partialFac));

  bfLogTodo("reminder: do a sorted insert in addPartialFac\n");

  bfPtrArrayAppend(&facStreamer->partialFacs, partialFac);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }
}

/* Notes:
 * - The rows and columns of `mat` should already be permuted into the
 *   orders defined by `facStreamer->rowTree` and
 *   `facStreamer->colTree`.
 *
 * TODO: this function is a mess! should be cleaned up! */
void bfFacStreamerFeed(BfFacStreamer *facStreamer, BfMat *Phi) {
  BF_ERROR_BEGIN();

  /* Get the current column tree leaf node */
  BfTreeNode const *colNode = bfFacStreamerGetCurrentColumnNode(facStreamer);
  BF_ASSERT(bfTreeNodeIsLeaf(colNode));

  /* Make sure `mat` has the right number of rows */
  if (bfMatGetNumRows(Phi) != bfFacStreamerGetNumRows(facStreamer))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* Make sure `mat` has the right number of columns */
  if (bfMatGetNumCols(Phi) != bfTreeNodeGetNumPoints(colNode))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfPtrArray PsiBlocks;
  bfInitPtrArrayWithDefaultCapacity(&PsiBlocks);
  HANDLE_ERROR();

  BfPtrArray WBlocks;
  bfInitPtrArrayWithDefaultCapacity(&WBlocks);
  HANDLE_ERROR();

  BfConstPtrArray rowNodes;
  bfConstPtrArrayInitWithDefaultCapacity(&rowNodes);
  HANDLE_ERROR();

  /* Create the stack used to adaptively find the cut through the row
   * tree to start the adaptive butterfly factorization */
  BfPtrArray stack = bfTreeGetLevelPtrArray(
    facStreamer->facSpec->rowTree, facStreamer->facSpec->rowTreeInitDepth);
  HANDLE_ERROR();

  /* Reverse the contents of the stack to make sure we traverse it
   * from the top of the column block down. */
  bfPtrArrayReverse(&stack);

  while (!bfPtrArrayIsEmpty(&stack)) {
    BfTreeNode const *rowNode = bfPtrArrayPopLast(&stack);

    BfMat *Psi = NULL, *W = NULL;
    bool metTol = getPsiAndW(facStreamer->facSpec, Phi, rowNode, &Psi, &W);
    HANDLE_ERROR();

    /* Accumulate the Psi and W blocks and continue if we successfully
     * compressed the current block */
    if (metTol) {
      bfPtrArrayAppend(&PsiBlocks, Psi);
      HANDLE_ERROR();

      bfPtrArrayAppend(&WBlocks, W);
      HANDLE_ERROR();

      bfConstPtrArrayAppend(&rowNodes, rowNode);
      HANDLE_ERROR();

      continue;
    }

    BF_ASSERT(Psi == NULL);
    BF_ASSERT(W == NULL);
    BF_ASSERT(!bfTreeNodeIsLeaf(rowNode));

    /* Push the children of the current row node onto the stack in
     * reverse order so that we traverse `mat` top to bottom */
    for (BfSize i = rowNode->maxNumChildren; i > 0; --i) {
      BfTreeNode const *childRowNode = rowNode->child[i - 1];
      if (childRowNode != NULL)
        bfPtrArrayAppend(&stack, (BfPtr)childRowNode);
    }
  }

  /* Create a new leaf node BF. When we do this, we steal the blocks
   * from PsiBlocks and WBlocks, so that we don't need to worry about
   * freeing them from this scope. */
  BfFac *partialFac = makeLeafNodePartialFac(
    colNode, &PsiBlocks, &WBlocks, BF_POLICY_STEAL, BF_POLICY_STEAL);
  HANDLE_ERROR();

  // TODO: this should go into makeLeafNodePartialFac
  for (BfSize i = 0; i < bfConstPtrArraySize(&rowNodes); ++i) {
    bfConstNodeArrayAppend(&partialFac->rowNodes, bfConstPtrArrayGet(&rowNodes, i));
    HANDLE_ERROR();
  }

  addPartialFac(facStreamer, partialFac);
  HANDLE_ERROR();

  if (facStreamer->facSpec->compareRelativeErrors)
    addPrevPhi(facStreamer, colNode, Phi);

  bfLogInfo("streamed [%lu, %lu)\n",
            bfTreeNodeGetFirstIndex(colNode),
            bfTreeNodeGetLastIndex(colNode));

#if BF_DEBUG
  FILE *fp = fopen("Psi.txt", "w");
  bfMatPrintBlocksDeep(partialFac->Psi, fp, 0, 0, 0);
  fclose(fp);
  fp = fopen("W0.txt", "w");
  bfMatPrintBlocksDeep(partialFac->W[0], fp, 0, 0, 0);
  fclose(fp);
#endif

  if (facStreamer->facSpec->compareRelativeErrors)
    checkRelError(Phi, partialFac);

  /* We're done with this node---move to the next one before
   * continuing to factorize. */
  bfTreeIterNext(facStreamer->colTreeIter);
  HANDLE_ERROR();

  continueFactorizing(facStreamer);

  BF_ERROR_END() {
    BF_DIE();
  }

  for (BfSize i = 0; i < bfPtrArraySize(&PsiBlocks); ++i) {
    BfMat *PsiBlock = bfPtrArrayGet(&PsiBlocks, i);
    bfMatDelete(&PsiBlock);
  }
  bfPtrArrayDeinit(&PsiBlocks);

  for (BfSize i = 0; i < bfPtrArraySize(&WBlocks); ++i) {
    BfMat *WBlock = bfPtrArrayGet(&WBlocks, i);
    bfMatDelete(&WBlock);
  }
  bfPtrArrayDeinit(&WBlocks);

  bfConstPtrArrayDeinit(&rowNodes);
  bfPtrArrayDeinit(&stack);
}

bool bfFacStreamerIsDone(BfFacStreamer const *facStreamer) {
  return bfTreeIterIsDone(facStreamer->colTreeIter);
}

BfFac *bfFacStreamerGetFac(BfFacStreamer const *facStreamer) {
  BF_ERROR_BEGIN();

  BfFac *fac = NULL;

  BfSize numFacs = bfPtrArraySize(&facStreamer->partialFacs);
  if (numFacs != 1)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  fac = bfPtrArrayGet(&facStreamer->partialFacs, 0);

  BF_ERROR_END() {}

  return fac;
}

BfFacSpan *bfFacStreamerGetFacSpan(BfFacStreamer const *facStreamer) {
  return bfFacSpanNewFromPtrArray(&facStreamer->partialFacs);
}

BfTreeNode *bfFacStreamerGetCurrentColumnNode(BfFacStreamer const *facStreamer) {
  return bfTreeIterGetCurrentNode(facStreamer->colTreeIter);
}

BfPerm const *bfFacStreamerGetRowTreeReversePerm(BfFacStreamer const *facStreamer) {
  return facStreamer->rowTreeReversePerm;
}
