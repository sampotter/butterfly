#include <bf/fac_helm2.h>

#include <math.h>

#include <bf/assert.h>
#include <bf/circle.h>
#include <bf/error_macros.h>
#include <bf/helm2.h>
#include <bf/mat_block.h>
#include <bf/mat_block_coo.h>
#include <bf/mat_block_dense.h>
#include <bf/mat_block_diag.h>
#include <bf/mat_product.h>
#include <bf/mem.h>
#include <bf/points.h>
#include <bf/quadtree.h>
#include <bf/quadtree_node.h>
#include <bf/util.h>

static BfSize const MAX_DENSE_MATRIX_SIZE = 128*128;

#if BF_DEBUG
static void facAuxDelete(BfFacAux *facAux) {
  bfPoints2DeinitAndDealloc(&facAux->srcPts[0]);
  bfPoints2DeinitAndDealloc(&facAux->srcPts[1]);
  bfPoints2DeinitAndDealloc(&facAux->tgtPts);
}
#endif

static BfSize
getChildren(BfQuadtreeNode const *node, BfQuadtreeNode const *child[4]) {
  for (BfSize i = 0; i < 4; ++i)
    child[i] = NULL;

  BfSize numChildren = 0;
  for (BfSize i = 0; i < 4; ++i)
    if (node->super.child[i] != NULL)
      child[numChildren++] = bfTreeNodeConstToQuadtreeNodeConst(node->super.child[i]);
  return numChildren;
}

static BfMat *makeFirstFactor(BfHelm2 const *helm,
                              BfQuadtree const *srcTree,
                              BfPtrArray const *srcLevelNodes,
                              BfPtrArray const *tgtLevelNodes) {
  BF_ERROR_BEGIN();

  /* there should only be one target node at the starting level of
   * the target tree when we begin the traversal */
  BF_ASSERT(bfPtrArraySize(tgtLevelNodes) == 1);

  /* the current level of the source node tree shouldn't be empty */
  BF_ASSERT(!bfPtrArrayIsEmpty(srcLevelNodes));

  BfQuadtreeNode const *srcNode = NULL;
  BfQuadtreeNode const *tgtNode = NULL;

  BfHelm2 helmProxy = *helm;
  helmProxy.layerPot = BF_PROXY_LAYER_POT[helm->layerPot];

  /* get the lone target node and its bounding circle */
  bfPtrArrayGetFirst(tgtLevelNodes, (BfPtr *)&tgtNode);
  BfCircle tgtCirc = bfQuadtreeNodeGetBoundingCircle(tgtNode);

  /* the number of diagonal blocks in the first factor equals the
   * number of nodes at the current depth in the source tree */
  BfSize numBlocks = bfPtrArraySize(srcLevelNodes);

  /* initialize the first block diagonal butterfly factor */
  BfMatBlockDiag *mat = bfMatBlockDiagNew();
  bfMatBlockDiagInit(mat, numBlocks, numBlocks);
  HANDLE_ERROR();

  /* iterate over each node at the starting level of the source tree
   * and compute the matrix which will maps the charges at the
   * original points in each source node to equivalent sources on a
   * circle surrounding the source box
   *
   * in the following loop, we initially set `rowOffset[i + 1]` and
   * `colOffset[i + 1]` to the number of rows in each block row and
   * column */
  BfPoints2 *srcPts = NULL, *srcCircPts = NULL, *tgtCircPts = NULL;
  BfVectors2 *srcNormals = NULL, *srcCircNormals = NULL;
  for (BfSize i = 0; i < numBlocks; ++i) {
    /* get the current source node and its bounding circle */
    srcNode = bfPtrArrayGet(srcLevelNodes, i);
    BfCircle srcCirc = bfQuadtreeNodeGetBoundingCircle(srcNode);

    /* get the original points in the current source node box */
    srcPts = bfQuadtreeNodeGetPoints(srcNode, srcTree);
    HANDLE_ERROR();

    /* Get the original normals in the current source node box: */
    srcNormals = bfQuadtreeNodeGetUnitNormals(srcNode, srcTree);
    HANDLE_ERROR();

    /* verify that the source bounding circle contains the points */
    BF_ASSERT(bfCircle2ContainsPoints(&srcCirc, srcPts));

    /* get the rank estimate for the current pair of source and
     * target bounding circles */
    BfSize p = bfHelm2RankEstForTwoCircles(helm, &srcCirc, &tgtCirc, 1, 1e-15);

    /* sample points on the source circle */
    srcCircPts = bfCircle2SamplePoints(&srcCirc, p);
    HANDLE_ERROR();

    /* sample normals on the source circle */
    srcCircNormals = bfCircle2SampleUnitNormals(&srcCirc, p);
    HANDLE_ERROR();

    /* sample points on the target circle */
    tgtCircPts = bfCircle2SamplePoints(&tgtCirc, p);
    HANDLE_ERROR();

    /* compute the shift matrix and store it in the current block */
    mat->super.block[i] = bfHelm2GetReexpansionMatrix(
      &helmProxy, srcPts, srcCircPts, srcNormals, srcCircNormals, tgtCircPts);
    HANDLE_ERROR();

    /* continue initializing the row and column offsets */
    mat->super.rowOffset[i + 1] = bfMatGetNumRows(mat->super.block[i]);
    mat->super.colOffset[i + 1] = bfMatGetNumCols(mat->super.block[i]);

    /* if we're debugging, store this block's points---free them
     * otherwise */
#if BF_DEBUG
    BfFacAux *facAux = bfMemAlloc(1, sizeof(BfFacAux));
    facAux->srcPts[0] = srcPts;
    facAux->srcPts[1] = srcCircPts;
    facAux->tgtPts = tgtCircPts;
    mat->super.block[i]->aux = facAux;
    mat->super.block[i]->auxDelete = (BfAuxDeleteFunc)facAuxDelete;
#else
    bfFreePoints2(&srcPts);
    bfFreePoints2(&tgtCircPts);
    bfFreePoints2(&srcCircPts);
#endif

    if (srcNormals != NULL)
      bfVectors2DeinitAndDealloc(&srcNormals);

    if (srcCircNormals != NULL)
      bfVectors2DeinitAndDealloc(&srcCircNormals);
  }

  BF_ERROR_END() {
    bfPoints2DeinitAndDealloc(&srcPts);
    bfPoints2DeinitAndDealloc(&tgtCircPts);
    bfPoints2DeinitAndDealloc(&srcCircPts);
    bfMatBlockDiagDeinitAndDealloc(&mat);
  }

  /* finish initializing row and column offsets */
  mat->super.rowOffset[0] = mat->super.colOffset[0] = 0;
  bfSizeRunningSum(numBlocks + 1, mat->super.rowOffset);
  bfSizeRunningSum(numBlocks + 1, mat->super.colOffset);

  return bfMatBlockDiagToMat(mat);
}

static BfSize getTotalNumChildren(BfPtrArray const *levelNodes) {
  BfSize totalNumChildren = 0;
  for (BfSize i = 0; i < bfPtrArraySize(levelNodes); ++i) {
    BfTreeNode const *node = bfPtrArrayGet(levelNodes, i);
    totalNumChildren += bfTreeNodeGetNumChildren(node);
  }
  return totalNumChildren;
}

/* An iterator which makes it easy to iterate over an array of
 * internal nodes and their children. The nodes contained in
 * `levelNodes` are the parents, and a precondition is that they are
 * all internal nodes (i.e., have children/aren't leaves). */
typedef struct {
  BfPtrArray const *levelNodes;
  BfSize nodeIndex;
  BfSize childIndex;
  BfSize numChildren;
  BfQuadtreeNode const *node;
  BfQuadtreeNode const *children[4];
  BfQuadtreeNode const *child;
  BfCircle circ;
  BfCircle childCirc;
} MakeFactorIter;

static void resetMakeFactorIter(MakeFactorIter *iter, BfPtrArray const *levelNodes) {
#if BF_DEBUG
  for (BfSize i = 0; i < bfPtrArraySize(levelNodes); ++i)
    BF_ASSERT(!bfTreeNodeIsLeaf(bfPtrArrayGet(levelNodes, i)));
#endif

  iter->levelNodes = levelNodes;
  iter->nodeIndex = iter->childIndex = 0;
  iter->node = bfPtrArrayGet(iter->levelNodes, iter->nodeIndex);
  iter->numChildren = getChildren(iter->node, iter->children);
  iter->child = iter->children[iter->childIndex];
  iter->circ = bfQuadtreeNodeGetBoundingCircle(iter->node);
  iter->childCirc = bfQuadtreeNodeGetBoundingCircle(iter->child);
}

static bool makeFactorIterNext(MakeFactorIter *iter) {
  if (++iter->childIndex == iter->numChildren) {
    if (++iter->nodeIndex == bfPtrArraySize(iter->levelNodes))
      return false;
    iter->childIndex = 0;
    iter->node = bfPtrArrayGet(iter->levelNodes, iter->nodeIndex);
    iter->numChildren = getChildren(iter->node, iter->children);
    iter->circ = bfQuadtreeNodeGetBoundingCircle(iter->node);
  }
  BF_ASSERT(iter->childIndex < iter->numChildren);
  iter->child = iter->children[iter->childIndex];
  iter->childCirc = bfQuadtreeNodeGetBoundingCircle(iter->child);
  return true;
}

/* This function computes an inner factor in a butterfly
 * factorization.
 *
 * TODO: explain how it works...
 */
static BfMat *makeFactor(BfHelm2 const *helm, BfMat const *prevMat, BfPtrArray const *srcLevelNodes, BfPtrArray const *tgtLevelNodes) {
  BF_ERROR_BEGIN();

  /* neither the source nor target levels should be empty */
  BF_ASSERT(!bfPtrArrayIsEmpty(srcLevelNodes));
  BF_ASSERT(!bfPtrArrayIsEmpty(tgtLevelNodes));

  /* get the layer potential we should use for reexpansion */
  BfHelm2 helmProxy = *helm;
  helmProxy.layerPot = BF_PROXY_LAYER_POT[helm->layerPot];

  /* count number of target nodes on this level */
  BfSize totalNumTgtNodes = bfPtrArraySize(tgtLevelNodes);

  /* count the total number of target children on this level */
  BfSize totalNumTgtChildren = getTotalNumChildren(tgtLevelNodes);

  /* count number of source nodes on this level */
  BfSize totalNumSrcNodes = bfPtrArraySize(srcLevelNodes);

  /* count the total number of source children on this level */
  BfSize totalNumSrcChildren = getTotalNumChildren(srcLevelNodes);

  BfMatBlock const *prevMatBlock = bfMatConstToMatBlockConst(prevMat);

  /* compute the number of block rows and columns as well as the total
   * number of blocks from this level's layout */
  BfSize numBlocks = totalNumSrcChildren*totalNumTgtChildren;
  BfSize numBlockRows = totalNumTgtChildren*totalNumSrcNodes;
  BfSize numBlockCols = totalNumSrcChildren*totalNumTgtNodes;
  BF_ASSERT(numBlockCols == bfMatBlockGetNumRowBlocks(prevMatBlock));

  BfMatBlockCoo *mat = bfMatBlockCooNew();
  bfMatBlockCooInit(mat, numBlockRows, numBlockCols, numBlocks);
  HANDLE_ERROR();

  /* set the number of columns in each block column to equal the
   * number of rows in each block row of the previous factor  */
  for (BfSize j = 0; j < numBlockCols; ++j)
    mat->super.colOffset[j + 1] = bfMatBlockGetNumBlockRows(prevMatBlock, j);

  /* ... and compute their running sum to get the column offsets */
  mat->super.colOffset[0] = 0;
  bfSizeRunningSum(numBlockCols + 1, mat->super.colOffset);

  /* next, we set the number of block rows and block columns
   *
   * we need to do this as an intermediate step since the individual
   * rank estimates for each pair of source and target circles could
   * vary slightly---we need the number of rows and columns of each
   * block to be compatible. hence, we use the corresponding maximum
   * rank estimate. */

  MakeFactorIter tgtIter, srcIter;

  BfSize blockIndex = 0, i_offset = 0, i, j;

  /* iterate over all target child nodes */
  resetMakeFactorIter(&tgtIter, tgtLevelNodes);
  do {
    /* iterate over all source child nodes */
    resetMakeFactorIter(&srcIter, srcLevelNodes);
    j = tgtIter.nodeIndex*totalNumSrcChildren;
    do {
      i = i_offset + srcIter.nodeIndex;
      BF_ASSERT(i < numBlockRows);
      BF_ASSERT(j < numBlockCols);

      /* a priori rank estimate for the original circles */
      BfSize rankOr = bfHelm2RankEstForTwoCircles(
        helm, &srcIter.childCirc, &tgtIter.circ, 1, 1e-15);

      /* a priori rank estimate for the new circles */
      BfSize rankEq = bfHelm2RankEstForTwoCircles(
        helm, &srcIter.circ, &tgtIter.childCirc, 1, 1e-15);

      /* use the larger of the two rank estimates... not sure if
       * this is totally necessary, probably being a little
       * paranoid... they should be nearly the same, but might
       * differ a little due to rounding */
      BfSize rank = rankOr > rankEq ? rankOr : rankEq;

      /* update number of rows for current block */
      BfSize rowOffset = mat->super.rowOffset[i + 1];
      if (rowOffset == BF_SIZE_BAD_VALUE || rank > rowOffset)
        mat->super.rowOffset[i + 1] = rank;

      /* set block row and column indices */
      BF_ASSERT(blockIndex < mat->numBlocks);
      mat->rowInd[blockIndex] = i;
      mat->colInd[blockIndex] = j;

      ++blockIndex;
      ++j;
    } while (makeFactorIterNext(&srcIter));
    i_offset += bfPtrArraySize(srcLevelNodes);
  } while (makeFactorIterNext(&tgtIter));

  /* compute running sum of row sizes to get the row offsets */
  mat->super.rowOffset[0] = 0;
  bfSizeRunningSum(numBlockRows + 1, mat->super.rowOffset);

  /* finally, we traverse the current source and target levels again,
   * sample proxy points, and compute shift matrices to assemble the
   * butterfly factor */

  BfPoints2 *srcChildPts, *srcPts, *tgtChildPts;
  BfVectors2 *srcChildNormals, *srcNormals;

  blockIndex = 0;

  /* iterate over all target child nodes */
  resetMakeFactorIter(&tgtIter, tgtLevelNodes);
  do {
    /* iterate over all source child nodes */
    resetMakeFactorIter(&srcIter, srcLevelNodes);
    do {
      BfSize numRows = bfMatBlockGetNumBlockRows(&mat->super, mat->rowInd[blockIndex]);
      BfSize numCols = bfMatBlockGetNumBlockCols(&mat->super, mat->colInd[blockIndex]);
      BF_ASSERT(numRows > 0 && numCols > 0);

      /* sample points on the source child circle */
      srcChildPts = bfCircle2SamplePoints(&srcIter.childCirc, numCols);
      HANDLE_ERROR();

      /* sample unit normals on the source child circle */
      srcChildNormals = bfCircle2SampleUnitNormals(&srcIter.childCirc, numCols);
      HANDLE_ERROR();

      /* sample points on the source circle */
      srcPts = bfCircle2SamplePoints(&srcIter.circ, numRows);
      HANDLE_ERROR();

      /* sample unit normals on the source circle */
      srcNormals = bfCircle2SampleUnitNormals(&srcIter.circ, numRows);
      HANDLE_ERROR();

      /* sample points on the target child circle */
      tgtChildPts = bfCircle2SamplePoints(&tgtIter.childCirc, numRows);
      HANDLE_ERROR();

      /* compute the shift matrix for this configuration of circles */
      mat->super.block[blockIndex] = bfHelm2GetReexpansionMatrix(
        &helmProxy, srcChildPts, srcPts, srcChildNormals, srcNormals, tgtChildPts);
      HANDLE_ERROR();

      /* if we're debugging, store this block's points---free them
       * otherwise */
#if BF_DEBUG
      BfFacAux *facAux = bfMemAlloc(1, sizeof(BfFacAux));
      facAux->srcPts[0] = srcChildPts;
      facAux->srcPts[1] = srcPts;
      facAux->tgtPts = tgtChildPts;
      mat->super.block[blockIndex]->aux = facAux;
      mat->super.block[blockIndex]->auxDelete = (BfAuxDeleteFunc)facAuxDelete;
#else
      bfFreePoints2(&srcChildPts);
      bfFreePoints2(&srcPts);
      bfFreePoints2(&tgtChildPts);
#endif

      if (srcChildNormals != NULL)
        bfVectors2DeinitAndDealloc(&srcChildNormals);

      if (srcNormals != NULL)
        bfVectors2DeinitAndDealloc(&srcNormals);

      ++blockIndex;
    } while (makeFactorIterNext(&srcIter));
  } while (makeFactorIterNext(&tgtIter));

  BF_ERROR_END() {
    bfPoints2DeinitAndDealloc(&srcChildPts);
    bfPoints2DeinitAndDealloc(&srcPts);
    bfPoints2DeinitAndDealloc(&tgtChildPts);
    bfMatBlockCooDeinitAndDealloc(&mat);
  }

  return bfMatBlockCooToMat(mat);
}

static BfMat *makeLastFactor(BfHelm2 const *helm, BfQuadtree const *tgtTree, BfMat const *prevMat, BfPtrArray const *srcLevelNodes, BfPtrArray const *tgtLevelNodes) {
  BF_ERROR_BEGIN();

  /* the current level of the target node tree shouldn't be empty */
  BF_ASSERT(!bfPtrArrayIsEmpty(tgtLevelNodes));

  /* the should only be one source node when we make the last
   * butterfly factor */
  BF_ASSERT(bfPtrArraySize(srcLevelNodes) == 1);

  /* the number of diagonal blocks in the last factor equals the
   * number of nodes at the current depth of the target tree ... */
  BfSize numBlocks = bfPtrArraySize(tgtLevelNodes);

  /* ... and this number of blocks should match the number of block
   *  rows of the previous factor */
  BF_ASSERT(numBlocks == prevMat->numRows);

  BfMatBlock const *prevMatBlock = bfMatConstToMatBlockConst(prevMat);

  BfQuadtreeNode const *srcNode = NULL;
  BfQuadtreeNode const *tgtNode = NULL;

  /* get the lone source node and its bounding circle */
  bfPtrArrayGetFirst(srcLevelNodes, (BfPtr *)&srcNode);
  BfCircle srcCirc = bfQuadtreeNodeGetBoundingCircle(srcNode);

  /* initialize the last butterfly factor */
  BfMatBlockDiag *mat = bfMatBlockDiagNew();
  bfMatBlockDiagInit(mat, numBlocks, numBlocks);
  HANDLE_ERROR();

  BfPoints2 *srcCircPts = NULL, *tgtPts = NULL;
  BfVectors2 *srcNormals = NULL, *tgtNormals = NULL;

  /* iterate over each node of the final level of the target node tree
   * and compute the matrix which will evaluate the potential at each
   * target point due to the charges on the source proxy points */
  for (BfSize i = 0; i < numBlocks; ++i) {
    BfSize prevNumRows = bfMatBlockGetNumBlockRows(prevMatBlock, i);

    /* get the current target node */
    tgtNode = bfPtrArrayGet(tgtLevelNodes, i);

    /* get the proxy points on the source circle for the current
     * source block */
    srcCircPts = bfCircle2SamplePoints(&srcCirc, prevNumRows);
    HANDLE_ERROR();

    /* sample normals on the source circle */
    if (BF_LAYER_POT_USES_SRC_NORMALS[helm->layerPot])
      srcNormals = bfCircle2SampleUnitNormals(&srcCirc, prevNumRows);

    /* get the current set of target points */
    tgtPts = bfQuadtreeNodeGetPoints(tgtNode, tgtTree);

    /* get the current target points' unit normals */
    if (BF_LAYER_POT_USES_TGT_NORMALS[helm->layerPot])
      tgtNormals = bfQuadtreeNodeGetUnitNormals(tgtNode, tgtTree);

    mat->super.block[i] = bfHelm2GetKernelMatrix(
      helm, srcCircPts, tgtPts, srcNormals, tgtNormals);
    HANDLE_ERROR();

    BF_ASSERT(mat->super.rowOffset[i + 1] == BF_SIZE_BAD_VALUE);
    mat->super.rowOffset[i + 1] = mat->super.block[i]->numRows;

    BF_ASSERT(mat->super.colOffset[i + 1] == BF_SIZE_BAD_VALUE);
    mat->super.colOffset[i + 1] = mat->super.block[i]->numCols;

    /* hang onto this block's points if we're in debug mode, and free
     * them otherwise */
#if BF_DEBUG
    BfFacAux *facAux = bfMemAlloc(1, sizeof(BfFacAux));
    facAux->srcPts[0] = bfPoints2NewWithDefaultCapacity();
    facAux->srcPts[1] = srcCircPts;
    facAux->tgtPts = tgtPts;
    mat->super.block[i]->aux = facAux;
    mat->super.block[i]->auxDelete = (BfAuxDeleteFunc)facAuxDelete;
#else
    bfPoints2DeinitAndDealloc(&srcCircPts);
    bfPoints2DeinitAndDealloc(&tgtPts);
#endif

    if (BF_LAYER_POT_USES_SRC_NORMALS[helm->layerPot])
      bfVectors2DeinitAndDealloc(&srcNormals);

    if (BF_LAYER_POT_USES_TGT_NORMALS[helm->layerPot])
      bfVectors2DeinitAndDealloc(&tgtNormals);
  }

  /* compute running sum of rows to get row offsets */
  mat->super.rowOffset[0] = 0;
  bfSizeRunningSum(numBlocks + 1, mat->super.rowOffset);

  /* compute running sum of columns to get column offsets */
  mat->super.colOffset[0] = 0;
  bfSizeRunningSum(numBlocks + 1, mat->super.colOffset);

  BF_ERROR_END() {
    bfPoints2DeinitAndDealloc(&srcCircPts);
    bfPoints2DeinitAndDealloc(&tgtPts);
    bfMatBlockDiagDeinitAndDealloc(&mat);
  }

  return bfMatBlockDiagToMat(mat);
}

static bool
allRankEstimatesAreOK(BfHelm2 const *helm,
                      BfQuadtreeNode const *tgtNode,
                      BfPtrArray const *srcLevelNodes)
{
  BfCircle tgtCirc = bfQuadtreeNodeGetBoundingCircle(tgtNode);

  for (BfSize i = 0; i < bfPtrArraySize(srcLevelNodes); ++i) {
    BfQuadtreeNode *srcNode = bfPtrArrayGet(srcLevelNodes, i);

    BfSize numSrcPoints = bfTreeNodeGetNumPoints(&srcNode->super);
    BfCircle srcCirc = bfQuadtreeNodeGetBoundingCircle(srcNode);
    BfSize rank = bfHelm2RankEstForTwoCircles(helm, &tgtCirc, &srcCirc, 1, 1e-15);

    if (rank > numSrcPoints)
      return false;
  }

  return true;
}

/* Prepare a set of quadtree level iterators for computing a butterfly
 * factorization of the kernel matrix mapping from the points in
 * `srcNode` to the points in `tgtNode`. The newly initialized
 * iterators will be pointed to by `srcLevelIter` and `tgtLevelIter`
 * upon successful completion, otherwise they will both equal
 * `NULL`. The number of butterfly factors is also returned. This
 * function will return 0 if it determined that it was not possible to
 * butterfly factorize the kernel matrix.
 *
 * If this function returns successfully, it is the caller's
 * responsibility to clean up `srcLevelIter` and `tgtLevelIter`.
 *
 * (NOTE: the algorithm used to determine whether or not a kernel
 * matrix is butterfliable is very crude at the moment, and may result
 * in false negatives.)
 *
 * TODO: this algorithm is a mess---we should carefully define what we
 * need to do here and implement it more concisely. This function has
 * been a big source of bugs and nuisances. */
BfSize
bfFacHelm2Prepare(BfHelm2 const *helm,
                  BfQuadtreeNode const *srcNode,
                  BfQuadtreeNode const *tgtNode,
                  BfTreeLevelIter *srcLevelIter,
                  BfTreeLevelIter *tgtLevelIter)
{
  BF_ERROR_BEGIN();

  BfSize numFactors = 0;

  BfSize numSrcNodes = bfTreeNodeGetNumPoints(&srcNode->super);
  BfSize numTgtNodes = bfTreeNodeGetNumPoints(&tgtNode->super);

  /* set up the level iterator for the source tree---this iterator
   * goes from the leaves of the tree to the root (in reverse) */
  bfTreeLevelIterInit(
    srcLevelIter,
    BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER,
    bfQuadtreeNodeToTreeNode((BfQuadtreeNode *)srcNode));
  HANDLE_ERROR();

  /* set up the level iterator for the target tree, which goes from
   * the root to the leaves */
  bfTreeLevelIterInit(
    tgtLevelIter,
    BF_TREE_TRAVERSAL_LR_LEVEL_ORDER,
    bfQuadtreeNodeToTreeNode((BfQuadtreeNode *)tgtNode));
  HANDLE_ERROR();

  /* find the deepest, complete, internal level of the target tree
   * using the target level iterator */
  BfSize maxAllowableDepthBelowTgtNode = bfTreeLevelIterCurrentDepth(tgtLevelIter);
  BF_ASSERT(bfTreeLevelIterCurrentLevelIsInternal(tgtLevelIter));
  bfTreeLevelIterNext(tgtLevelIter);
  HANDLE_ERROR();
  while (bfTreeLevelIterCurrentLevelIsInternal(tgtLevelIter)) {
    ++maxAllowableDepthBelowTgtNode;
    bfTreeLevelIterNext(tgtLevelIter);
    HANDLE_ERROR();
  }

  /* reset the target tree level iterator */
  bfTreeLevelIterDeinit(tgtLevelIter);
  bfTreeLevelIterInit(
    tgtLevelIter,
    BF_TREE_TRAVERSAL_LR_LEVEL_ORDER,
    bfQuadtreeNodeToTreeNode((BfQuadtreeNode *)tgtNode));
  HANDLE_ERROR();

  BF_ASSERT(bfTreeLevelIterGetNumPoints(tgtLevelIter) == numTgtNodes);

  /* get the current source and target depths and make sure they're
   * compatible */
  BfSize currentSrcDepth =
    bfTreeLevelIterCurrentDepth(srcLevelIter);
  BfSize currentTgtDepth = bfTreeLevelIterCurrentDepth(tgtLevelIter);
  BF_ASSERT(currentTgtDepth <= currentSrcDepth);

  /* skip source levels until we're no deeper than the maximum depth
   * beneath the target node */
  while (currentSrcDepth > maxAllowableDepthBelowTgtNode) {
    bfTreeLevelIterNext(srcLevelIter);
    --currentSrcDepth;
  }

  /* Skip levels until we're on a complete level of tree */
  while (bfTreeLevelIterGetNumPoints(srcLevelIter) != numSrcNodes) {
    bfTreeLevelIterNext(srcLevelIter);
    --currentSrcDepth;
  }
  BF_ASSERT(bfTreeLevelIterGetNumPoints(srcLevelIter) == numSrcNodes);

  /* Skip levels until we're on an "internal" level of the tree */
  while (!bfTreeLevelIterCurrentLevelIsInternal(srcLevelIter)) {
    bfTreeLevelIterNext(srcLevelIter);
    --currentSrcDepth;
  }

  /* TODO: step the source level iterator until:
   * - the rank estimate between the first target node and each source
   *   node is smaller than corresponding the number of points */

  while (currentSrcDepth > currentTgtDepth &&
         !allRankEstimatesAreOK(helm, tgtNode, srcLevelIter->levelNodes)) {
    bfTreeLevelIterNext(srcLevelIter);
    --currentSrcDepth;
  }

  /* get number of factors in the butterfly factorization... if we
   * can't butterfly this matrix, return 0 to signal this */
  numFactors = allRankEstimatesAreOK(helm, tgtNode, srcLevelIter->levelNodes) ?
      currentSrcDepth - currentTgtDepth + 2 : 0;

  BF_ERROR_END() {
    bfTreeLevelIterDeinit(srcLevelIter);
    bfTreeLevelIterDeinit(tgtLevelIter);
  }

  return numFactors;
}

BfMatProduct *bfFacHelm2Make(BfHelm2 const *helm, BfQuadtree const *srcTree, BfQuadtree const *tgtTree, BfTreeLevelIter *srcLevelIter, BfTreeLevelIter *tgtLevelIter, BfSize numFactors) {
  BF_ERROR_BEGIN();

  /* allocate space for the butterfly factors
   *
   * when multiplying, the factors will be multiplied in the order
   * that they're stored here (e.g., bf_factors[0] is the first factor
   * that will be multiplied) */
  BfMat **factor = bfMemAlloc(1, numFactors*sizeof(BfMat *));
  if (factor == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* make the first factor in the butterfly factorization: this is the
   * factor which initially shifts the charges on the source points to
   * the first level of source circles */
  factor[0] = makeFirstFactor(helm, srcTree, srcLevelIter->levelNodes, tgtLevelIter->levelNodes);
  HANDLE_ERROR();

  for (BfSize i = 1; i < numFactors - 1; ++i) {
    /* go up a level on the source tree */
    bfTreeLevelIterNext(srcLevelIter);

    /* make the next factor */
    factor[i] = makeFactor(
      helm, factor[i - 1], srcLevelIter->levelNodes, tgtLevelIter->levelNodes);
    HANDLE_ERROR();

    /* go down a level on the target tree */
    bfTreeLevelIterNext(tgtLevelIter);
  }

  /* make the last factor in the butterfly factorization: this is the
   * evaluation factor, which computes the potential at each target
   * point due to the charges on the final source circle */
  factor[numFactors - 1] = makeLastFactor(
    helm, tgtTree, factor[numFactors - 2], srcLevelIter->levelNodes,
    tgtLevelIter->levelNodes);
  HANDLE_ERROR();

  BfMatProduct *prod = bfMatProductNew();
  bfMatProductInit(prod);
  for (BfSize i = 0; i < numFactors; ++i)
    bfMatProductPostMultiply(prod, factor[numFactors - 1 - i]);

  BF_ERROR_END() {
    BF_DIE(); // TODO: need to think carefully about how to do this
  }

  bfMemFree(factor);

  return prod;
}

BfMat *bfFacHelm2MakeSingleLevel(BfHelm2 const *helm, BfQuadtreeNode const *srcNode, BfQuadtreeNode const *tgtNode) {
  BF_ERROR_BEGIN();

  BfMatProduct *fac = NULL;

  BfQuadtree const *srcTree = bfQuadtreeNodeGetQuadtreeConst(srcNode);
  BfQuadtree const *tgtTree = bfQuadtreeNodeGetQuadtreeConst(tgtNode);

  BfTreeLevelIter srcLevelIter, tgtLevelIter;
  BfSize numFactors = bfFacHelm2Prepare(helm, srcNode, tgtNode, &srcLevelIter, &tgtLevelIter);
  HANDLE_ERROR();

  fac = bfFacHelm2Make(helm, srcTree, tgtTree, &srcLevelIter, &tgtLevelIter, numFactors);
  HANDLE_ERROR();

  bfTreeLevelIterDeinit(&srcLevelIter);
  bfTreeLevelIterDeinit(&tgtLevelIter);

  BF_ERROR_END() {
    BF_DIE();
  }

  return bfMatProductToMat(fac);
}

/* Return the child nodes of `node` in a `BfPtrArray`. */
static BfPtrArray getChildrenAsPtrArray(BfQuadtreeNode const *node) {
  BfPtrArray childNodes;
  bfInitPtrArray(&childNodes, 4);
  for (BfSize i = 0; i < 4; ++i)
    if (node->super.child[i])
      bfPtrArrayAppend(&childNodes, node->super.child[i]);
  return childNodes;
}

static BfMat *facHelm2MakeMultilevel_dense(BfHelm2 const *helm, BfQuadtree const *srcTree, BfQuadtree const *tgtTree, BfQuadtreeNode const *srcNode, BfQuadtreeNode const *tgtNode) {
  BF_ERROR_BEGIN();

  BfPoints2 *srcPts = NULL, *tgtPts = NULL;
  BfVectors2 *srcNormals = NULL, *tgtNormals = NULL;

  BfMat *Z = NULL;

  srcPts = bfQuadtreeNodeGetPoints(srcNode, srcTree);
  tgtPts = bfQuadtreeNodeGetPoints(tgtNode, tgtTree);

  if (BF_LAYER_POT_USES_SRC_NORMALS[helm->layerPot])
    srcNormals = bfQuadtreeNodeGetUnitNormals(srcNode, srcTree);

  if (BF_LAYER_POT_USES_TGT_NORMALS[helm->layerPot])
    tgtNormals = bfQuadtreeNodeGetUnitNormals(tgtNode, tgtTree);

  Z = bfHelm2GetKernelMatrix(helm, srcPts, tgtPts, srcNormals, tgtNormals);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  bfPoints2DeinitAndDealloc(&srcPts);
  bfPoints2DeinitAndDealloc(&tgtPts);

  if (BF_LAYER_POT_USES_SRC_NORMALS[helm->layerPot])
    bfVectors2DeinitAndDealloc(&srcNormals);

  if (BF_LAYER_POT_USES_TGT_NORMALS[helm->layerPot])
    bfVectors2DeinitAndDealloc(&tgtNormals);

  return Z;
}

static
BfMat *
facHelm2MakeMultilevel_separated(BfHelm2 const *helm,
                                 BfQuadtree const *srcTree,
                                 BfQuadtree const *tgtTree,
                                 BfQuadtreeNode const *srcNode,
                                 BfQuadtreeNode const *tgtNode) {
  BF_ERROR_BEGIN();

  BfMat *fac = NULL;

  BfTreeLevelIter srcLevelIter, tgtLevelIter;
  BfSize numFactors = bfFacHelm2Prepare(helm, srcNode, tgtNode, &srcLevelIter, &tgtLevelIter);

  fac = numFactors == 0 ?
    facHelm2MakeMultilevel_dense(helm, srcTree, tgtTree, srcNode, tgtNode) :
    bfMatProductToMat(bfFacHelm2Make(helm, srcTree, tgtTree, &srcLevelIter, &tgtLevelIter, numFactors));
  HANDLE_ERROR();

  bfTreeLevelIterDeinit(&srcLevelIter);
  bfTreeLevelIterDeinit(&tgtLevelIter);

  BF_ERROR_END() {
    BF_DIE();
  }

  return fac;
}

static void facHelm2MakeMultilevel_rec(BfHelm2 const *helm,
                                       BfQuadtree const *srcTree,
                                       BfQuadtree const *tgtTree,
                                       BfPtrArray const *srcNodes,
                                       BfPtrArray const *tgtNodes,
                                       BfSize level,
                                       BfMatBlockDense *blockMat);

static
BfMat *
facHelm2MakeMultilevel_diag(BfHelm2 const *helm,
                            BfQuadtree const *srcTree,
                            BfQuadtree const *tgtTree,
                            BfQuadtreeNode const *srcNode,
                            BfQuadtreeNode const *tgtNode,
                            BfSize level) {
  BF_ERROR_BEGIN();

  BfPtrArray srcChildNodes, tgtChildNodes;

  BfMat *mat = NULL;
  BfMatBlockDense *childBlockMat = NULL;

  srcChildNodes = getChildrenAsPtrArray(srcNode);
  HANDLE_ERROR();

  tgtChildNodes = getChildrenAsPtrArray(tgtNode);
  HANDLE_ERROR();

  childBlockMat = bfMatBlockDenseNew();
  HANDLE_ERROR();

  mat = bfMatBlockDenseToMat(childBlockMat);

  BfSize numBlockRows = bfPtrArraySize(&tgtChildNodes);
  BfSize numBlockCols = bfPtrArraySize(&srcChildNodes);
  bfMatBlockDenseInit(childBlockMat, numBlockRows, numBlockCols);
  HANDLE_ERROR();

  facHelm2MakeMultilevel_rec(helm, srcTree, tgtTree, &srcChildNodes, &tgtChildNodes,
                             level + 1, childBlockMat);
  HANDLE_ERROR();

  BF_ASSERT(bfMatGetNumRows(mat) == bfTreeNodeGetNumPoints(&tgtNode->super));
  BF_ASSERT(bfMatGetNumCols(mat) == bfTreeNodeGetNumPoints(&srcNode->super));

  BF_ERROR_END() {}

  bfPtrArrayDeinit(&srcChildNodes);
  bfPtrArrayDeinit(&tgtChildNodes);

  return bfMatBlockDenseToMat(childBlockMat);
}

static void facHelm2MakeMultilevel_rec(BfHelm2 const *helm,
                                       BfQuadtree const *srcTree,
                                       BfQuadtree const *tgtTree,
                                       BfPtrArray const *srcNodes,
                                       BfPtrArray const *tgtNodes,
                                       BfSize level,
                                       BfMatBlockDense *blockMat) {
  BF_ERROR_BEGIN();

  BfMatBlock *super = &blockMat->super;
  BfMat *mat = NULL;

  super->rowOffset[0] = 0;
  super->colOffset[0] = 0;

  BfSize numRowBlocks = bfPtrArraySize(tgtNodes);
  BfSize numColBlocks = bfPtrArraySize(srcNodes);

  for (BfSize i = 0; i < numRowBlocks; ++i) {
    BfQuadtreeNode *tgtNode = bfPtrArrayGet(tgtNodes, i);
    BfSize numRows = bfTreeNodeGetNumPoints(&tgtNode->super);

    for (BfSize j = 0; j < numColBlocks; ++j) {
      BfQuadtreeNode *srcNode = bfPtrArrayGet(srcNodes, j);
      BfSize numCols = bfTreeNodeGetNumPoints(&srcNode->super);

      bool separated = bfQuadtreeNodesAreSeparated(srcNode, tgtNode);

      if (numRows*numCols < MAX_DENSE_MATRIX_SIZE)
        mat = facHelm2MakeMultilevel_dense(helm, srcTree, tgtTree, srcNode, tgtNode);
      else if (separated)
        mat = facHelm2MakeMultilevel_separated(helm, srcTree, tgtTree, srcNode, tgtNode);
      else
        /* TODO: we really need to consolidate _rec and _diag (also,
         * "_diag" is a total misnomer) */
        mat = facHelm2MakeMultilevel_diag(helm, srcTree, tgtTree, srcNode, tgtNode, level);

      if (bfMatGetNumRows(mat) != numRows)
        RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

      if (bfMatGetNumCols(mat) != numCols)
        RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

      if (mat == NULL)
        RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

      /* set the current block now that we've computed it
       *
       * TODO: add an init function for MatBlockDense which we can use
       * to do all of this more safely... not great that we need to
       * reach into MatBlockDense to do this right now without hitting
       * the guardrails in SetBlock... */
      // bfMatBlockDenseSetBlock(blockMat, i, j, mat);
      blockMat->super.block[i*blockMat->super.super.numCols + j] = mat;
      HANDLE_ERROR();

      /* update the column block offset if we're making our first pass
       * through a row of blocks */
      if (i == 0)
        super->colOffset[j + 1] = super->colOffset[j] + numCols;
    }

    /* update the row block offset */
    super->rowOffset[i + 1] = super->rowOffset[i] + numRows;
  }

#if BF_DEBUG
  for (BfSize i = 0; i < numRowBlocks; ++i) {
    BfSize m = super->rowOffset[i + 1] - super->rowOffset[i];
    for (BfSize j = 0; j < numColBlocks; ++j) {
      BfSize n = super->colOffset[j + 1] - super->colOffset[j];
      BfMat *block = bfMatBlockDenseGetBlock(blockMat, i, j);
      BF_ASSERT(bfMatGetNumRows(block) == m);
      BF_ASSERT(bfMatGetNumCols(block) == n);
    }
  }
#endif

  BF_ERROR_END() {
    BF_DIE();
  }
}

BfMat *bfFacHelm2MakeMultilevel(BfHelm2 const *helm, BfQuadtree const *srcTree, BfQuadtree const *tgtTree) {
  BF_ERROR_BEGIN();

  // TODO: ugh. Should actually be passing these TreeLevelIters to
  // this function, not the trees themselves. Can fix this up later.
  //
  // An alternative might be to pass the trees and to *not* skip the
  // levels, since this would automatically happen in the recursion
  // anyway (level 2 is the first level with any well-separated
  // nodes). Or maybe it would start at level 1 occassionally? Hm.

  /** Get nodes on second level of source (column) tree: */

  BfTreeLevelIter srcLevelIter;
  bfTreeLevelIterInit(
    &srcLevelIter, BF_TREE_TRAVERSAL_LR_LEVEL_ORDER,
    bfQuadtreeNodeToTreeNode((BfQuadtreeNode *)srcTree->super.root));
  HANDLE_ERROR();

  bfTreeLevelIterNext(&srcLevelIter);
  bfTreeLevelIterNext(&srcLevelIter);
  HANDLE_ERROR();

  BfPtrArray const *srcLevelNodes = srcLevelIter.levelNodes;
  BfSize numSrcNodes = bfPtrArraySize(srcLevelNodes);

  /** Get nodes on second level of target (row) tree: */

  BfTreeLevelIter tgtLevelIter;
  bfTreeLevelIterInit(
    &tgtLevelIter, BF_TREE_TRAVERSAL_LR_LEVEL_ORDER,
    bfQuadtreeNodeToTreeNode((BfQuadtreeNode *)tgtTree->super.root));
  HANDLE_ERROR();

  bfTreeLevelIterNext(&tgtLevelIter);
  bfTreeLevelIterNext(&tgtLevelIter);
  HANDLE_ERROR();

  BfPtrArray const *tgtLevelNodes = tgtLevelIter.levelNodes;
  BfSize numTgtNodes = bfPtrArraySize(tgtLevelNodes);

  /** Build the multilevel butterfly factorization: */

  BfMatBlockDense *matBlockDense = bfMatBlockDenseNew();
  bfMatBlockDenseInit(matBlockDense, numTgtNodes, numSrcNodes);
  HANDLE_ERROR();

  facHelm2MakeMultilevel_rec(helm, srcTree, tgtTree, srcLevelNodes, tgtLevelNodes, 2, matBlockDense);
  HANDLE_ERROR();

  BF_ERROR_END() {
    bfMatBlockDenseDeinitAndDealloc(&matBlockDense);
    BF_DIE();
  }

  bfTreeLevelIterDeinit(&srcLevelIter);
  bfTreeLevelIterDeinit(&tgtLevelIter);

  return bfMatBlockDenseToMat(matBlockDense);
}
