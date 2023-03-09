#include <bf/fac_streamer.h>

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include <bf/const.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/indexed_mat.h>
#include <bf/linalg.h>
#include <bf/mat_block_coo.h>
#include <bf/mat_block_dense.h>
#include <bf/mat_block_diag.h>
#include <bf/mat_dense_real.h>
#include <bf/mat_diag_real.h>
#include <bf/mat_identity.h>
#include <bf/ptr_array.h>
#include <bf/size_array.h>
#include <bf/tree.h>
#include <bf/tree_node.h>
#include <bf/tree_iter_post_order.h>
#include <bf/vec_real.h>

#undef BF_ARRAY_DEFAULT_CAPACITY
#define BF_ARRAY_DEFAULT_CAPACITY 32

////////////////////////////////////////////////////////////////////////
// BEGIN PARTIAL FAC STUFF...

typedef struct {
  BfTreeNode const *colNode;

  BfConstPtrArray rowNodes;

  BfMat *Psi;

  BfSize numW;
  BfMat **W;
} PartialFac;

static PartialFac *makeLeafNodePartialFac(BfTreeNode const *colNode,
                                          BfPtrArray *PsiBlocks,
                                          BfPtrArray *WBlocks) {
  BEGIN_ERROR_HANDLING();

  PartialFac *partialFac = malloc(sizeof(PartialFac));
  if (partialFac == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  partialFac->colNode = colNode;

  bfConstPtrArrayInitWithDefaultCapacity(&partialFac->rowNodes);
  HANDLE_ERROR();

  BfMatBlockDiag *Psi = bfMatBlockDiagNewFromBlocks(PsiBlocks);
  HANDLE_ERROR();

  partialFac->Psi = bfMatBlockDiagToMat(Psi);

  partialFac->numW = 1;

  partialFac->W = malloc(sizeof(BfMatBlockCoo *));
  if (partialFac->W == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfMatBlockCoo *W0 = bfMatBlockCooNewColFromBlocks(WBlocks);
  HANDLE_ERROR();

  partialFac->W[0] = bfMatBlockCooToMat(W0);

  END_ERROR_HANDLING() {
    partialFac = NULL;
  }

  return partialFac;
}

void setPartialFacWBlock(PartialFac *partialFac, BfSize k, BfMat *W) {
  assert(k < partialFac->numW);
  assert(W != NULL);
  partialFac->W[k] = W;
}

typedef struct {
  BfTreeNode const *treeNode;
  BfMat const *mat;
} MatWithTreeNodeKey;

struct BfFacStreamer {
  BfSize rowTreeInitDepth;
  BfSize colTreeInitDepth;

  BfTree *rowTree;
  BfTree *colTree;

  BfTreeIter *colTreeIter;

  /* The tolerance used to compute truncated SVDs of blocks when
   * streaming the butterfly factorization. */
  BfReal tol;

  /* The minimum number of rows in a block needed before a truncated
   * SVD is attempted. This is used to bottom out when finding
   * epsilon-rank cuts in the row tree. */
  BfSize minNumRows;

  /* The minimum number of columns in a block needed before trying a
   * truncated SVD. This is used to prevent us from starting to low in
   * the column tree. */
  BfSize minNumCols;

  /* This is basically an association list which maps from column tree
   * nodes to lists of Psi and W blocks. We use this to track blocks
   * as we traverse the column tree, merging and splitting blocks. */
  BfPtrArray partialFacs;

#if BF_DEBUG
  /* An association list mapping from column tree nodes to contiguous
   * blocks of Phi. We use this to check how accurate each partial
   * factorization is as we stream them. */
  BfPtrArray *prevPhis;
#endif
};

static BfSize partialFacGetNumRows(PartialFac const *partialFac) {
  return bfMatGetNumRows(partialFac->Psi);
}

static void addPartialFac(BfFacStreamer *facStreamer, PartialFac *partialFac) {
  BEGIN_ERROR_HANDLING();

  assert(bfFacStreamerGetNumRows(facStreamer) == partialFacGetNumRows(partialFac));

  // printf("reminder: do a sorted insert in addPartialFac\n");

  bfPtrArrayAppend(&facStreamer->partialFacs, partialFac);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    assert(false);
  }
}

static BfVec *partialFacMulVec(PartialFac *partialFac, BfVec const *x) {
  assert(partialFac->numW > 0);
  BfSize k = partialFac->numW - 1;
  BfVec *y = NULL;
  BfVec *yPrev = bfMatMulVec(partialFac->W[k], x);
  while (k > 0) {
    y = bfMatMulVec(partialFac->W[--k], yPrev);
    bfVecDelete(&yPrev);
    yPrev = y;
  }
  y = bfMatMulVec(partialFac->Psi, yPrev);
  bfVecDelete(&yPrev);
  return y;
}

// END PARTIAL FAC STUFF...
////////////////////////////////////////////////////////////////////////

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
  facStreamer->minNumRows = spec->minNumRows;
  facStreamer->minNumCols = spec->minNumCols;

  BfTreeIterPostOrder *iter = bfTreeIterPostOrderNew();
  HANDLE_ERROR();

  bfTreeIterPostOrderInit(iter, facStreamer->colTree);
  HANDLE_ERROR();

  facStreamer->colTreeIter = bfTreeIterPostOrderToTreeIter(iter);

  facStreamer->partialFacs = bfGetUninitializedPtrArray();
  HANDLE_ERROR();

  facStreamer->prevPhis = bfPtrArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  bfInitPtrArrayWithDefaultCapacity(&facStreamer->partialFacs);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    bfFacStreamerDeinit(facStreamer);
  }
}

void bfFacStreamerDeinit(BfFacStreamer *facStreamer) {
  (void)facStreamer;
  assert(false);
}

BfSize bfFacStreamerGetNumRows(BfFacStreamer const *facStreamer) {
  return bfTreeGetNumPoints(facStreamer->rowTree);
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

static bool getPsiAndW(BfFacStreamer const *facStreamer,
                       BfMat const *mat, BfTreeNode const *rowNode,
                       BfMat **PsiPtr, BfMat **WPtr) {
  BEGIN_ERROR_HANDLING();

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
  if (bfMatGetNumCols(mat) < facStreamer->minNumCols)
    success = getPsiAndW_skinny(block, PsiPtr, WPtr);

  /* If there are too few rows, we pass through the current block as W
   * and emit an identity matrix for Psi. */
  else if (i1 - i0 < facStreamer->minNumRows) {
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
    success = getPsiAndW_normal(block, facStreamer->tol, PsiPtr, WPtr);

  END_ERROR_HANDLING() {
    success = false;

    /* TODO: make sure to free Psi and W if I allocated them... */
    assert(false);
  }

  bfMatDelete(&block);

  return success;
}

/* Find the current partial factorizations associated with each of the
 * children of the current column node and return them in PtrArray. */
static BfPtrArray getCurrentPartialFacs(BfFacStreamer const *facStreamer) {
  BEGIN_ERROR_HANDLING();

  BfTreeNode const *currentColNode = bfFacStreamerGetCurrentColumnNode(facStreamer);
  assert(!bfTreeNodeIsLeaf(currentColNode));

  BfPtrArray currentPartialFacs;
  bfInitPtrArrayWithDefaultCapacity(&currentPartialFacs);
  HANDLE_ERROR();

  for (BfSize k = 0; k < currentColNode->maxNumChildren; ++k) {
    BfTreeNode const *childColNode = currentColNode->child[k];
    if (childColNode == NULL)
      continue;

    PartialFac *partialFac = NULL;
    for (BfSize l = 0; l < bfPtrArraySize(&facStreamer->partialFacs); ++l) {
      partialFac = bfPtrArrayGet(&facStreamer->partialFacs, l);
      if (partialFac->colNode == childColNode)
        break;
      partialFac = NULL;
    }
    assert(partialFac != NULL);

    bfPtrArrayAppend(&currentPartialFacs, partialFac);
    HANDLE_ERROR();
  }

  assert(!bfPtrArrayIsEmpty(&currentPartialFacs));

  END_ERROR_HANDLING() {
    assert(false);
  }

  return currentPartialFacs;
}

static bool partialFacHasContiguousRowSpan(PartialFac const *partialFac) {
  BfSize numRowNodes = bfConstPtrArraySize(&partialFac->rowNodes);
  if (numRowNodes <= 1)
    return true;
  BfTreeNode const *rowNode = bfConstPtrArrayGet(&partialFac->rowNodes, 0);
  BfSize i1Prev = bfTreeNodeGetLastIndex(rowNode);
  for (BfSize k = 1; k < numRowNodes; ++k) {
    rowNode = bfConstPtrArrayGet(&partialFac->rowNodes, k);
    BfSize i0 = bfTreeNodeGetFirstIndex(rowNode);
    if (i0 != i1Prev)
      return false;
    i1Prev = bfTreeNodeGetLastIndex(rowNode);
  }
  return true;
}

static bool partialFacsHaveContinguousRowSpans(BfPtrArray const *partialFacs) {
  for (BfSize k = 0; k < bfPtrArraySize(partialFacs); ++k) {
    PartialFac const *partialFac = bfPtrArrayGet(partialFacs, k);
    if (!partialFacHasContiguousRowSpan(partialFac))
      return false;
  }
  return true;
}

static bool partialFacsHaveSameRowSpan(BfPtrArray const *partialFacs) {
  assert(bfPtrArraySize(partialFacs) >= 2);

  if (!partialFacsHaveContinguousRowSpans(partialFacs)) {
    // TODO: this case is a little annoying and complicated and
    // doesn't matter for us right now. Handle it later if we need to
    // for some reason.
    assert(false);
  }

  PartialFac *partialFac;
  BfTreeNode const *rowNode;

  BfSize i0, i1;

  bfPtrArrayGetFirst(partialFacs, (BfPtr *)&partialFac);

  rowNode = bfConstPtrArrayGetFirst(&partialFac->rowNodes);
  i0 = bfTreeNodeGetFirstIndex(rowNode);

  rowNode = bfConstPtrArrayGetLast(&partialFac->rowNodes);
  i1 = bfTreeNodeGetLastIndex(rowNode);

  for (BfSize k = 1; k < bfPtrArraySize(partialFacs); ++k) {
    BfSize i0_, i1_;

    bfPtrArrayGetFirst(partialFacs, (BfPtr *)&partialFac);

    rowNode = bfConstPtrArrayGetFirst(&partialFac->rowNodes);
    i0_ = bfTreeNodeGetFirstIndex(rowNode);

    rowNode = bfConstPtrArrayGetLast(&partialFac->rowNodes);
    i1_ = bfTreeNodeGetLastIndex(rowNode);

    if (i0 != i0_ || i1 != i1_)
      return false;
  }

  return true;
}

/* Given a `PtrArray` filled with `PartialFac`s, get the first row
 * tree node in each `PartialFac`, collect them together into a
 * `ConstPtrArray`, and return it. */
static BfConstPtrArray getFirstRowNodes(BfPtrArray const *partialFacs) {
  BEGIN_ERROR_HANDLING();

  BfConstPtrArray rowNodes;
  bfConstPtrArrayInitWithDefaultCapacity(&rowNodes);
  HANDLE_ERROR();

  for (BfSize k = 0; k < bfPtrArraySize(partialFacs); ++k) {
    PartialFac const *partialFac = bfPtrArrayGet(partialFacs, k);

    if (bfConstPtrArrayIsEmpty(&partialFac->rowNodes))
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

    BfTreeNode const *rowNode = bfConstPtrArrayGetFirst(&partialFac->rowNodes);

    bfConstPtrArrayAppend(&rowNodes, rowNode);
  }

  END_ERROR_HANDLING() {
    bfConstPtrArrayDeinit(&rowNodes);
  }

  return rowNodes;
}

/* Given a `PtrArray` filled with `PartialFac`s, get the last row
 * tree node in each `PartialFac`, collect them together into a
 * `ConstPtrArray`, and return it. */
static BfConstPtrArray getLastRowNodes(BfPtrArray const *partialFacs) {
  BEGIN_ERROR_HANDLING();

  BfConstPtrArray rowNodes;
  bfConstPtrArrayInitWithDefaultCapacity(&rowNodes);
  HANDLE_ERROR();

  for (BfSize k = 0; k < bfPtrArraySize(partialFacs); ++k) {
    PartialFac const *partialFac = bfPtrArrayGet(partialFacs, k);

    if (bfConstPtrArrayIsEmpty(&partialFac->rowNodes))
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

    BfTreeNode const *rowNode = bfConstPtrArrayGetLast(&partialFac->rowNodes);

    bfConstPtrArrayAppend(&rowNodes, rowNode);
  }

  END_ERROR_HANDLING() {
    bfConstPtrArrayDeinit(&rowNodes);
  }

  return rowNodes;
}

static bool nodesHaveSameFirstIndex(BfConstPtrArray const *nodes) {
  BEGIN_ERROR_HANDLING();

  BfSize numNodes = bfConstPtrArraySize(nodes);

  if (numNodes == 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  bool sameFirstIndex = true;

  BfTreeNode const *node = bfConstPtrArrayGetFirst(nodes);
  BfSize i0 = bfTreeNodeGetFirstIndex(node);
  for (BfSize k = 1; k < numNodes; ++k) {
    node = bfConstPtrArrayGet(nodes, k);
    if (bfTreeNodeGetFirstIndex(node) != i0) {
      sameFirstIndex = false;
      break;
    }
  }

  END_ERROR_HANDLING() {}

  return sameFirstIndex;
}

static BfSize getMaxLastIndexForRowNodes(BfConstPtrArray const *nodes,
                                         BfTreeNode const **argmaxNodePtr) {
  BEGIN_ERROR_HANDLING();

  BfSize i1Max = BF_SIZE_BAD_VALUE;
  BfTreeNode const *argmaxNode = NULL;

  BfSize numNodes = bfConstPtrArraySize(nodes);
  if (numNodes == 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  argmaxNode = bfConstPtrArrayGetFirst(nodes);
  i1Max = bfTreeNodeGetLastIndex(argmaxNode);
  for (BfSize k = 1; k < numNodes; ++k) {
    BfTreeNode const *node = bfConstPtrArrayGet(nodes, k);
    BfSize i1 = bfTreeNodeGetLastIndex(node);
    if (i1 > i1Max) {
      i1Max = i1;
      argmaxNode = node;
    }
  }

  END_ERROR_HANDLING() {
    i1Max = BF_SIZE_BAD_VALUE;
    argmaxNode = NULL;
  }

  if (argmaxNodePtr != NULL)
    *argmaxNodePtr = argmaxNode;

  return i1Max;
}

static BfTreeNode const *getNodeByFirstIndex(BfConstPtrArray const *nodes, BfSize i0) {
  BfTreeNode const *node = NULL;
  for (BfSize k = 0; k < bfConstPtrArraySize(nodes); ++k) {
    node = bfConstPtrArrayGet(nodes, k);
    if (bfTreeNodeGetFirstIndex(node) == i0)
      break;
  }
  return node;
}

static BfConstPtrArray getRowNodesByFirstIndex(BfPtrArray const *partialFacs, BfSize i0) {
  BEGIN_ERROR_HANDLING();

  BfConstPtrArray rowNodes;
  bfConstPtrArrayInitWithDefaultCapacity(&rowNodes);
  HANDLE_ERROR();

  for (BfSize k = 0; k < bfPtrArraySize(partialFacs); ++k) {
    PartialFac const *partialFac = bfPtrArrayGet(partialFacs, k);

    BfTreeNode const *rowNode = getNodeByFirstIndex(&partialFac->rowNodes, i0);
    if (rowNode == NULL)
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

    bfConstPtrArrayAppend(&rowNodes, rowNode);
    HANDLE_ERROR();
  }

  END_ERROR_HANDLING() {
    bfConstPtrArrayDeinit(&rowNodes);
  }

  return rowNodes;
}

static BfConstPtrArray getMergeCut(BfPtrArray const *partialFacs) {
  BEGIN_ERROR_HANDLING();

  BfSize numPartialFacs = bfPtrArraySize(partialFacs);

  if (numPartialFacs == 0)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* For the "merge cut" operation to be well defined, each of the
   * factorizations must have the same row span. */
  if (!partialFacsHaveSameRowSpan(partialFacs))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfConstPtrArray mergeCut;
  bfConstPtrArrayInitWithDefaultCapacity(&mergeCut);
  HANDLE_ERROR();

  BfConstPtrArray rowNodes = getFirstRowNodes(partialFacs);
  HANDLE_ERROR();

  assert(bfConstPtrArraySize(&rowNodes) == numPartialFacs);
  assert(nodesHaveSameFirstIndex(&rowNodes));

  BfConstPtrArray lastRowNodes = getLastRowNodes(partialFacs);
  HANDLE_ERROR();

  BfTreeNode const *rowNode = NULL;
  BfSize i1 = getMaxLastIndexForRowNodes(&rowNodes, &rowNode);

  bfConstPtrArrayAppend(&mergeCut, rowNode);
  HANDLE_ERROR();

  BfSize i1Final = getMaxLastIndexForRowNodes(&lastRowNodes, NULL);

  while (i1 != i1Final) {
    bfConstPtrArrayDeinit(&rowNodes);

    rowNodes = getRowNodesByFirstIndex(partialFacs, i1);
    HANDLE_ERROR();

    i1 = getMaxLastIndexForRowNodes(&rowNodes, &rowNode);
    HANDLE_ERROR();

    bfConstPtrArrayAppend(&mergeCut, rowNode);
    HANDLE_ERROR();

    assert(bfConstPtrArraySize(&rowNodes) == numPartialFacs);
  }

  END_ERROR_HANDLING() {
    bfConstPtrArrayDeinit(&mergeCut);
  }

  return mergeCut;
}

static void
getPsiAndW0BlocksByRowNodeForPartialFac(PartialFac const *partialFac,
                                        BfTreeNode const *rowNode,
                                        BfMat **PsiBlockPtr,
                                        BfMat **W0BlockPtr) {
  BEGIN_ERROR_HANDLING();

  if (PsiBlockPtr == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (W0BlockPtr == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatBlockDiag const *Psi = bfMatToMatBlockDiag(partialFac->Psi);
  HANDLE_ERROR();

  BfMat const *W0 = partialFac->W[0];

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
   * nodes in `partialFac->rowNodes` whose index ranges are a subset
   * of [i0, i1). */
#if BF_DEBUG
  for (BfSize k = 0; k < bfConstPtrArraySize(&partialFac->rowNodes); ++k) {
    BfTreeNode const *otherRowNode = bfConstPtrArrayGet(&partialFac->rowNodes, k);
    BfSize i0_ = bfTreeNodeGetFirstIndex(otherRowNode);
    BfSize i1_ = bfTreeNodeGetLastIndex(otherRowNode);
    if (i0 <= i0_ && i1_ <= i1 && rowNode != otherRowNode)
      assert(bfTreeNodeIsDescendant(otherRowNode, rowNode));
  }
#endif

  BfSize numSubblocks = 0;

  for (BfSize k = 0; k < bfMatBlockDiagNumBlocks(Psi); ++k) {
    /* Get row offsets for current block of Psi */
    BfSize i0_ = Psi->super.rowOffset[k];
    BfSize i1_ = Psi->super.rowOffset[k + 1];

    /* Skip blocks which aren't indexed by [i0, i1) */
    if (!(i0 <= i0_ && i1_ <= i1))
      continue;

    ++numSubblocks;

    /* Get a copy of the current diagonal Psi block */
    BfMat *PsiSubblock = bfMatBlockDiagGetBlockCopy(Psi, k, k);
    HANDLE_ERROR();

    /* ... and append it to the running array of Psi subblocks */
    bfPtrArrayAppend(&PsiSubblocks, PsiSubblock);
    HANDLE_ERROR();

    /* Get column offsets for current block of Psi */
    BfSize j0_ = Psi->super.colOffset[k];
    BfSize j1_ = Psi->super.colOffset[k + 1];

    /* Get a copy of the current W row subblock */
    BfMat *W0Subblock = bfMatGetRowRangeCopy(W0, j0_, j1_);
    HANDLE_ERROR();

    /* .. and append it to the running array of W subblocks */
    bfPtrArrayAppend(&W0Subblocks, W0Subblock);
    HANDLE_ERROR();
  }

  assert(numSubblocks > 0);
  assert(bfPtrArraySize(&PsiSubblocks) == numSubblocks);
  assert(bfPtrArraySize(&W0Subblocks) == numSubblocks);

  /* If we only found one Psi and W subblock each, we should just
   * return these directly instead of wrapping them in a MatBlockDiag
   * and MatBlockCoo, respectively (what is done otherwise). */
  if (numSubblocks == 1) {
    bfPtrArrayGetFirst(&PsiSubblocks, (BfPtr *)&PsiBlock);
    bfPtrArrayGetFirst(&W0Subblocks, (BfPtr *)&W0Block);
  } else {
    /* Diagonally concatenate the Psi subblocks we found */
    PsiBlock = bfMatBlockDiagToMat(bfMatBlockDiagNewFromBlocks(&PsiSubblocks));
    HANDLE_ERROR();

    /* ... and vertically concatenate the W row subblocks */
    W0Block = bfMatBlockCooToMat(bfMatBlockCooNewColFromBlocks(&W0Subblocks));
    HANDLE_ERROR();
  }

  END_ERROR_HANDLING() {
    assert(false); // T_T
  }

  bfPtrArrayDeinit(&PsiSubblocks);
  bfPtrArrayDeinit(&W0Subblocks);

  assert(PsiBlock != NULL);
  assert(W0Block != NULL);
  assert(bfTreeNodeGetNumPoints(rowNode) == bfMatGetNumRows(PsiBlock));
  assert(bfMatGetNumCols(PsiBlock) == bfMatGetNumRows(W0Block));
  assert(bfMatGetNumCols(W0Block) == bfMatGetNumCols(partialFac->W[0]));

  *PsiBlockPtr = PsiBlock;
  *W0BlockPtr = W0Block;
}

static void getPsiAndWBlocksByRowNode(BfPtrArray const *currentPartialFacs,
                                      BfTreeNode const *rowNode,
                                      BfMat **PsiPtr,
                                      BfMat **WPtr) {
  BEGIN_ERROR_HANDLING();

  if (PsiPtr == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  if (WPtr == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfPtrArray PsiBlocks;
  bfInitPtrArrayWithDefaultCapacity(&PsiBlocks);
  HANDLE_ERROR();

  BfPtrArray W0Blocks;
  bfInitPtrArrayWithDefaultCapacity(&W0Blocks);
  HANDLE_ERROR();

  for (BfSize k = 0; k < bfPtrArraySize(currentPartialFacs); ++k) {
    PartialFac const *partialFac = bfPtrArrayGet(currentPartialFacs, k);

    /* Get diagonal psi blocks for kth partial fac and corresponding
     * index range by indexing with rowNode */
    BfMat *PsiBlock = NULL;
    BfMat *W0Block = NULL;
    getPsiAndW0BlocksByRowNodeForPartialFac(partialFac, rowNode, &PsiBlock, &W0Block);

    /* Append the new Psi block */
    bfPtrArrayAppend(&PsiBlocks, PsiBlock);
    HANDLE_ERROR();

    /* Append the new W block */
    bfPtrArrayAppend(&W0Blocks, W0Block);
    HANDLE_ERROR();
  }

  /* Horizontally concatenate together the Psi blocks we found to get
   * the leading "Psi*" block. */
  BfMatBlockCoo *Psi = bfMatBlockCooNewRowFromBlocks(&PsiBlocks);
  HANDLE_ERROR();

  /* Diagonally concatenate together the W row blocks we found to get
   * the first column block of the permuted W factor. */
  BfMatBlockDiag *W = bfMatBlockDiagNewFromBlocks(&W0Blocks);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {}

  *PsiPtr = bfMatBlockCooToMat(Psi);
  *WPtr = bfMatBlockDiagToMat(W);

  assert(bfMatGetNumCols(*PsiPtr) == bfMatGetNumRows(*WPtr));
}

static bool getLowRankApproximation(BfFacStreamer const *facStreamer,
                                    BfMat const *PsiStarSubblock,
                                    BfMat **PsiSubblockPtr,
                                    BfMat **W0SubblockPtr) {
  BEGIN_ERROR_HANDLING();

  BfMat *U = NULL;
  BfMatDiagReal *S = NULL;
  BfMat *VT = NULL;

  BfTruncSpec truncSpec = {
    .usingTol = true,
    .tol = facStreamer->tol
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

    BfPtrArray blocks;
    bfInitPtrArrayWithDefaultCapacity(&blocks);
    HANDLE_ERROR();

    for (BfSize k = 0; k < numNonzeroColumnRanges; ++k) {
      BfSize j0 = bfSizeArrayGet(nonzeroColumnRanges, 2*k);
      BfSize j1 = bfSizeArrayGet(nonzeroColumnRanges, 2*k + 1);
      BfMat *block = bfMatGetColRangeCopy(W0Subblock, j0, j1);
      BfIndexedMat *indexedBlock = malloc(sizeof(BfIndexedMat));
      indexedBlock->i0 = 0;
      indexedBlock->j0 = j0;
      indexedBlock->mat = block;
      bfPtrArrayAppend(&blocks, indexedBlock);
    }

    BfMatBlockCoo *W0SubblockSparsified = bfMatBlockCooNewFromIndexedBlocks(
      bfMatGetNumRows(W0Subblock), bfMatGetNumCols(W0Subblock), &blocks);
    HANDLE_ERROR();

    bfMatDelete(&W0Subblock);

    W0Subblock = bfMatBlockCooToMat(W0SubblockSparsified);
  }

  END_ERROR_HANDLING() {
    assert(false);
  }

  *PsiSubblockPtr = U;
  *W0SubblockPtr = W0Subblock;

  bfVecDelete(&s);
  bfMatDiagRealDeinitAndDealloc(&S);

  return truncated;
}

static void
findEpsilonRankCutAndGetNewBlocks(BfFacStreamer const *facStreamer,
                                  BfTreeNode const *rootRowNode,
                                  BfMat const *PsiStarBlock,
                                  BfConstPtrArray **epsRankCutPtr,
                                  BfMat **PsiBlockPtr,
                                  BfMat **W0BlockPtr) {
  BEGIN_ERROR_HANDLING();

  BfMatBlockDiag *PsiBlock = NULL;
  BfMatBlockCoo *W0Block = NULL;

  if (epsRankCutPtr == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (PsiBlockPtr == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (W0BlockPtr == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (bfMatGetNumRows(PsiStarBlock) != bfTreeNodeGetNumPoints(rootRowNode))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* The row nodes in the epsilon-rank cut. */
  BfConstPtrArray *epsRankCut = bfConstPtrArrayNewWithDefaultCapacity();
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
    assert(i0 <= i0_ && i0_ < i1_ && i1_ <= i1);

    /* Index range for the rows to copy: */
    BfSize i0__ = i0_ - i0;
    BfSize i1__ = i1_ - i0;
    assert(i0__ < i1__ && i1__ <= bfMatGetNumRows(PsiStarBlock));

    BfMat *PsiStarSubblock = bfMatGetRowRangeCopy(PsiStarBlock, i0__, i1__);
    HANDLE_ERROR();

    BfSize numRowsPsiStarSubblock = bfMatGetNumRows(PsiStarSubblock);
    BfSize numColsPsiStarSubblock = bfMatGetNumCols(PsiStarSubblock);

    BfMat *PsiSubblock = NULL;
    BfMat *W0Subblock = NULL;

    /* If the new Psi subblock has too few rows, we omit an identity
     * matrix for the Psi subblock and pass Psi through as the W0
     * subblock. */
    if (numRowsPsiStarSubblock < facStreamer->minNumRows) {
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
    if (numColsPsiStarSubblock < facStreamer->minNumCols) {
      PsiSubblock = PsiStarSubblock;

      BfMatIdentity *identity = bfMatIdentityNew();
      HANDLE_ERROR();

      bfMatIdentityInit(identity, numColsPsiStarSubblock);
      HANDLE_ERROR();

      W0Subblock = bfMatIdentityToMat(identity);

      goto next;
    }

    bool truncated = getLowRankApproximation(
      facStreamer, PsiStarSubblock, &PsiSubblock, &W0Subblock);
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
    bfConstPtrArrayAppend(epsRankCut, rowNode);
    HANDLE_ERROR();

    /* Append the Psi subblock. */
    bfPtrArrayAppend(&PsiSubblocks, PsiSubblock);
    HANDLE_ERROR();

    /* Append the W0 subblock. */
    bfPtrArrayAppend(&W0Subblocks, W0Subblock);
    HANDLE_ERROR();
  }

  /* Diagonally concatenate the Psi subblocks. */
  PsiBlock = bfMatBlockDiagNewFromBlocks(&PsiSubblocks);
  HANDLE_ERROR();

  /* Vertically concatenate the W0 subblocks. */
  W0Block = bfMatBlockCooNewColFromBlocks(&W0Subblocks);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    bfConstPtrArrayDeinitAndDealloc(&epsRankCut);

    // TODO: free blocks
    assert(false);
  }

  bfPtrArrayDeinit(&PsiSubblocks);
  bfPtrArrayDeinit(&W0Subblocks);
  bfConstPtrArrayDeinit(&stack);

  *epsRankCutPtr = epsRankCut;
  *PsiBlockPtr = bfMatBlockDiagToMat(PsiBlock);
  *W0BlockPtr = bfMatBlockCooToMat(W0Block);
}

static PartialFac *mergeAndSplit(BfFacStreamer *facStreamer) {
  static int CALL_NUMBER = 0;

  BEGIN_ERROR_HANDLING();

  PartialFac *mergedFac = NULL;

  BfMatBlockDiag *Psi = NULL;
  BfMatBlockDiag *W0 = NULL;
  BfMatBlockCoo *W1 = NULL; // TODO: change this to BfMatBlockDense

  // get PartialFacs for children of colNode
  BfPtrArray currentPartialFacs = getCurrentPartialFacs(facStreamer);
  HANDLE_ERROR();

  assert(!bfPtrArrayIsEmpty(&currentPartialFacs));

  /* TODO: For now, we assume that all current partial factorizations
   * have the same number of W factors. This can be relaxed by filling
   * in with identity matrices below when we concatenate together the
   * trailing W factors. */
  BfSize maxNumW;
  { PartialFac *fac = bfPtrArrayGet(&currentPartialFacs, 0);
    maxNumW = fac->numW;
    for (BfSize k = 1; k < bfPtrArraySize(&currentPartialFacs); ++k) {
      PartialFac *otherFac = bfPtrArrayGet(&currentPartialFacs, k);
      if (fac->numW != otherFac->numW)
        RAISE_ERROR(BF_ERROR_RUNTIME_ERROR); } }

  // compute merge cut...
  BfConstPtrArray mergeCut = getMergeCut(&currentPartialFacs);
  HANDLE_ERROR();

  /* Array used to accumulate row nodes from each of the epsilon rank
   * cuts found when merging the current partial factorizations
   * starting from the merge cut. The resulting span of nodes should
   * be beneath the merge cut and compatible with it. */
  BfConstPtrArray rowNodes;
  bfConstPtrArrayInitWithDefaultCapacity(&rowNodes);
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
  for (BfSize k = 0; k < bfConstPtrArraySize(&mergeCut); ++k) {
    BfTreeNode const *rowNode = bfConstPtrArrayGet(&mergeCut, k);

    BfConstPtrArray *epsRankCut = NULL;
    BfMat *PsiStarBlock = NULL;
    BfMat *PsiBlock = NULL;
    BfMat *W0Block = NULL;
    BfMat *W1Block = NULL;

    /* Get current row block of Psi. We simultaneously get the column
     * index ranges corresponding to each column subblock. */
    getPsiAndWBlocksByRowNode(&currentPartialFacs, rowNode, &PsiStarBlock, &W1Block);

    bfPtrArrayAppend(&W1Blocks, W1Block);
    HANDLE_ERROR();

    /* Find the epsilon-rank cut and compute the correspond blocks of
     * the new Psi and W0 factors while simultaneously sifting the W0
     * blocks from the current partial factorizations together to get
     * the next block of the W1 factor. */
    findEpsilonRankCutAndGetNewBlocks(
      facStreamer,
      rowNode,
      PsiStarBlock,
      &epsRankCut,
      &PsiBlock,
      &W0Block);
    HANDLE_ERROR();

    bfConstPtrArrayExtend(&rowNodes, epsRankCut);
    HANDLE_ERROR();

    bfPtrArrayAppend(&PsiBlocks, PsiBlock);
    HANDLE_ERROR();

    bfPtrArrayAppend(&W0Blocks, W0Block);
    HANDLE_ERROR();

    bfConstPtrArrayDeinitAndDealloc(&epsRankCut);
  }

  /* Diagonally concatenate the new Psi blocks to get the leading
   * factor of the merged butterfly factorization. */
  Psi = bfMatBlockDiagNewFromBlocks(&PsiBlocks);
  HANDLE_ERROR();

  /* Diagonally concatenate together the new W0 blocks. */
  W0 = bfMatBlockDiagNewFromBlocks(&W0Blocks);
  HANDLE_ERROR();

  /* Vertically concatenate the collected W1 blocks. */
  W1 = bfMatBlockCooNewColFromBlocks(&W1Blocks);
  HANDLE_ERROR();

  /* Create and initialize the merged factorization. */
  // TODO: should really break this out into a separate function once
  // we refactor PartialFac into an honest ButterflyFactorization
  {
    mergedFac = malloc(sizeof(PartialFac));
    if (mergedFac == NULL)
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

    mergedFac->colNode = bfTreeIterGetCurrentNode(facStreamer->colTreeIter);

    mergedFac->rowNodes = rowNodes;

    mergedFac->Psi = bfMatBlockDiagToMat(Psi);

    mergedFac->numW = maxNumW + 1;
    assert(mergedFac->numW >= 2);

    mergedFac->W = calloc(mergedFac->numW, sizeof(BfMat *));
    if (mergedFac->W == NULL)
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

    setPartialFacWBlock(mergedFac, 0, bfMatBlockDiagToMat(W0));
    setPartialFacWBlock(mergedFac, 1, bfMatBlockCooToMat(W1));

    /* Block together trailing W factors from the current partial
     * butterfly factorizations to get the trailing W factors for the
     * newly merged butterfly factorization */
    for (BfSize k = 1; k < maxNumW; ++k) {
      BfPtrArray WkBlocks;
      bfInitPtrArrayWithDefaultCapacity(&WkBlocks);
      HANDLE_ERROR();

      for (BfSize l = 0; l < bfPtrArraySize(&currentPartialFacs); ++l) {
        PartialFac const *partialFac = bfPtrArrayGet(&currentPartialFacs, l);

        BfMat *WkBlock = partialFac->W[k];
        assert(WkBlock != NULL);

        bfPtrArrayAppend(&WkBlocks, WkBlock);
        HANDLE_ERROR();
      }

      BfMatBlockDiag *Wk = bfMatBlockDiagNewFromBlocks(&WkBlocks);
      HANDLE_ERROR();

      setPartialFacWBlock(mergedFac, k + 1, bfMatBlockDiagToMat(Wk));

      bfPtrArrayDeinit(&WkBlocks);
    }

    /* Make sure everything was set successfully: */
#if BF_DEBUG
    assert(mergedFac->Psi != NULL);
    for (BfSize i = 0; i < bfConstPtrArraySize(&mergedFac->rowNodes); ++i)
      assert(bfConstPtrArrayGet(&mergedFac->rowNodes, i) != NULL);
    assert(mergedFac->colNode != NULL);
    for (BfSize k = 0; k < mergedFac->numW; ++k)
      assert(mergedFac->W[k] != NULL);
#endif
  }

  END_ERROR_HANDLING() {
    assert(false);
  }

  bfPtrArrayDeinit(&PsiBlocks);
  bfPtrArrayDeinit(&W0Blocks);
  bfPtrArrayDeinit(&W1Blocks);

  ++CALL_NUMBER;

  return mergedFac;
}

static BfMat const *getPhiByColNode(BfFacStreamer const *facStreamer,
                                    BfTreeNode const *treeNode) {
  BfMat const *mat = NULL;
  for (BfSize k = 0; k < bfPtrArraySize(facStreamer->prevPhis); ++k) {
    MatWithTreeNodeKey const *entry = bfPtrArrayGet(facStreamer->prevPhis, k);
    if (entry->treeNode == treeNode) {
      mat = entry->mat;
      break;
    }
  }
  return mat;
}

static void addPrevPhi(BfFacStreamer *facStreamer, BfTreeNode const *treeNode, BfMat const *Phi) {
  MatWithTreeNodeKey *entry = malloc(sizeof(MatWithTreeNodeKey));
  entry->treeNode = treeNode;
  entry->mat = bfMatCopy(Phi);
  bfPtrArrayAppend(facStreamer->prevPhis, entry);
}

static void addPrevPhiForChildNodes(BfFacStreamer *facStreamer, BfTreeNode const *currentColNode) {
  BfSize m = bfFacStreamerGetNumRows(facStreamer);
  BfSize n = bfTreeNodeGetNumPoints(currentColNode);
  BfMatDenseReal *Phi = bfMatDenseRealNewWithValue(m, n, BF_NAN);

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

    BfMatDenseReal const *PhiBlock = bfMatConstToMatDenseRealConst(
      getPhiByColNode(facStreamer, childColNode));

    bfMatDenseRealSetBlock(Phi, 0, m, j0, j1, PhiBlock);
  }

  addPrevPhi(facStreamer, currentColNode, bfMatDenseRealToMat(Phi));
}

static void continueFactorizing(BfFacStreamer *facStreamer) {
  BEGIN_ERROR_HANDLING();

  /* Continue the post-order traversal until the next leaf node */
  while (!bfTreeIterIsDone(facStreamer->colTreeIter)) {
    bfTreeIterNext(facStreamer->colTreeIter);
    HANDLE_ERROR();

    BfTreeNode const *currentColNode = bfTreeIterGetCurrentNode(facStreamer->colTreeIter);
    if (bfTreeNodeIsLeaf(currentColNode))
      break;

    printf("merging [%lu, %lu)",
           bfTreeNodeGetFirstIndex(currentColNode),
           bfTreeNodeGetLastIndex(currentColNode));

#if BF_DEBUG
    addPrevPhiForChildNodes(facStreamer, currentColNode);
#endif

    PartialFac *mergedFac = mergeAndSplit(facStreamer);
    HANDLE_ERROR();

    printf(" -> depths: Psi=%lu", bfMatBlockGetMaxDepth(bfMatToMatBlock(mergedFac->Psi)));
    for (BfSize k = 0; k < mergedFac->numW; ++k)
      printf(", W%lu=%lu", k, bfMatBlockGetMaxDepth(bfMatToMatBlock(mergedFac->W[k])));
    printf("\n");

    // printf("writing to Psi.txt");

    FILE *fp = fopen("Psi.txt", "w");
    bfMatPrintBlocksDeep(mergedFac->Psi, fp, 0, 0, 0);
    fclose(fp);

    for (BfSize k = 0; k < mergedFac->numW; ++k) {
      char path[256];
      sprintf(path, "W%lu.txt", k);

      // printf(", %s", path);

      fp = fopen(path, "w");
      bfMatPrintBlocksDeep(mergedFac->W[k], fp, 0, 0, 0);
      fclose(fp);
    }

    // printf("\n");

    bfPtrArrayAppend(&facStreamer->partialFacs, mergedFac);
    HANDLE_ERROR();

    // TODO: throw out old partial facs
#if BF_DEBUG
    // TODO: throw out old prevPhis
#endif
  }

  END_ERROR_HANDLING() {
    assert(false);
  }
}

/* Notes:
 * - The rows and columns of `mat` should already be permuted into the
 *   orders defined by `facStreamer->rowTree` and
 *   `facStreamer->colTree`.
 *
 * TODO: this function is a mess! should be cleaned up! */
void bfFacStreamerFeed(BfFacStreamer *facStreamer, BfMat const *Phi) {
  BEGIN_ERROR_HANDLING();

  /* Get the current column tree leaf node */
  BfTreeNode const *colNode = bfFacStreamerGetCurrentColumnNode(facStreamer);
  assert(bfTreeNodeIsLeaf(colNode));

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
    facStreamer->rowTree, facStreamer->rowTreeInitDepth);
  HANDLE_ERROR();

  /* Reverse the contents of the stack to make sure we traverse it
   * from the top of the column block down. */
  bfPtrArrayReverse(&stack);

  while (!bfPtrArrayIsEmpty(&stack)) {
    BfTreeNode const *rowNode = bfPtrArrayPopLast(&stack);

    BfMat *Psi = NULL, *W = NULL;
    bool metTol = getPsiAndW(facStreamer, Phi, rowNode, &Psi, &W);
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

    /* Push the children of the current row node onto the stack in
     * reverse order so that we traverse `mat` top to bottom */
    for (BfSize i = rowNode->maxNumChildren; i > 0; --i)
      bfPtrArrayAppend(&stack, rowNode->child[i - 1]);
  }

  PartialFac *partialFac = makeLeafNodePartialFac(colNode, &PsiBlocks, &WBlocks);
  HANDLE_ERROR();


  addPartialFac(facStreamer, partialFac);
  HANDLE_ERROR();

  printf("streamed [%lu, %lu)",
#if BF_DEBUG
  addPrevPhi(facStreamer, colNode, Phi);
#endif

         bfTreeNodeGetFirstIndex(colNode),
         bfTreeNodeGetLastIndex(colNode));
  printf(" -> depths: Psi=%lu", bfMatBlockGetMaxDepth(bfMatToMatBlock(partialFac->Psi)));
  for (BfSize k = 0; k < partialFac->numW; ++k)
    printf(", W%lu=%lu", k, bfMatBlockGetMaxDepth(bfMatToMatBlock(partialFac->W[k])));
  printf("\n");


  for (BfSize i = 0; i < bfConstPtrArraySize(&rowNodes); ++i) {
    bfConstPtrArrayAppend(&partialFac->rowNodes, bfConstPtrArrayGet(&rowNodes, i));
    HANDLE_ERROR();
  }

  continueFactorizing(facStreamer);

  END_ERROR_HANDLING() {
    assert(false);
  }

  bfPtrArrayDeinit(&PsiBlocks);
  bfPtrArrayDeinit(&WBlocks);
  bfConstPtrArrayDeinit(&rowNodes);
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
  return bfTreeIterGetCurrentNode(facStreamer->colTreeIter);
}
