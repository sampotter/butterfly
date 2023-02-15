#include <bf/fac.h>

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <bf/circle.h>
#include <bf/error_macros.h>
#include <bf/helm2.h>
#include <bf/mat_block.h>
#include <bf/mat_block_coo.h>
#include <bf/mat_block_dense.h>
#include <bf/mat_block_diag.h>
#include <bf/mat_product.h>
#include <bf/points.h>
#include <bf/quadtree_node.h>
#include <bf/util.h>

static BfSize const MAX_DENSE_MATRIX_SIZE = 128*128;

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

static BfMat *makeFirstFactor(BfReal K, BfLayerPotential layerPot,
                              BfQuadtree const *tree,
                              BfPtrArray const *srcLevelNodes,
                              BfPtrArray const *tgtLevelNodes) {
  BEGIN_ERROR_HANDLING();

  /* there should only be one target node at the starting level of
   * the target tree when we begin the traversal */
  assert(bfPtrArraySize(tgtLevelNodes) == 1);

  /* the current level of the source node tree shouldn't be empty */
  assert(!bfPtrArrayIsEmpty(srcLevelNodes));

  BfQuadtreeNode const *srcNode = NULL;
  BfQuadtreeNode const *tgtNode = NULL;

  BfLayerPotential proxyLayerPot = BF_PROXY_LAYER_POT[layerPot];

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
  BfPoints2 srcPts, tgtCircPts, srcCircPts;
  for (BfSize i = 0; i < numBlocks; ++i) {
    /* get the current source node and its bounding circle */
    srcNode = bfPtrArrayGet(srcLevelNodes, i);
    BfCircle srcCirc = bfQuadtreeNodeGetBoundingCircle(srcNode);

    /* get the original points in the current source node box */
    srcPts = bfQuadtreeNodeGetPoints(srcNode, tree);
    HANDLE_ERROR();

    /* verify that the source bounding circle contains the points */
    assert(bfCircle2ContainsPoints(&srcCirc, &srcPts));

    /* get the rank estimate for the current pair of source and
     * target bounding circles */
    BfSize p = bfHelm2RankEstForTwoCircles(&srcCirc, &tgtCirc, K, 1, 1e-15);

    /* sample points on the source circle */
    tgtCircPts = bfCircle2SamplePoints(&tgtCirc, p);
    HANDLE_ERROR();

    /* sample points on the target circle */
    srcCircPts = bfCircle2SamplePoints(&srcCirc, p);
    HANDLE_ERROR();

    /* compute the shift matrix and store it in the current block */
    mat->super.block[i] = bfHelm2GetReexpansionMatrix(
      &srcPts, &srcCircPts, &tgtCircPts, K, proxyLayerPot);
    HANDLE_ERROR();

    /* continue initializing the row and column offsets */
    mat->super.rowOffset[i + 1] = bfMatGetNumRows(mat->super.block[i]);
    mat->super.colOffset[i + 1] = bfMatGetNumCols(mat->super.block[i]);

    /* if we're debugging, store this block's points---free them
     * otherwise */
#if BF_DEBUG
    BfFacAux *facAux = malloc(sizeof(BfFacAux));
    facAux->srcPts[0] = srcPts;
    facAux->srcPts[1] = srcCircPts;
    facAux->tgtPts = tgtCircPts;
    mat->super.block[i]->aux = facAux;
#else
    bfFreePoints2(&srcPts);
    bfFreePoints2(&tgtCircPts);
    bfFreePoints2(&srcCircPts);
#endif
  }

  END_ERROR_HANDLING() {
    bfFreePoints2(&srcPts);
    bfFreePoints2(&tgtCircPts);
    bfFreePoints2(&srcCircPts);
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
    assert(!bfTreeNodeIsLeaf(bfPtrArrayGet(levelNodes, i)));
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
  assert(iter->childIndex < iter->numChildren);
  iter->child = iter->children[iter->childIndex];
  iter->childCirc = bfQuadtreeNodeGetBoundingCircle(iter->child);
  return true;
}

/* This function computes an inner factor in a butterfly
 * factorization.
 *
 * ... explain how this works ...
 */
static BfMat *makeFactor(BfMat const *prevMat, BfReal K, BfLayerPotential layerPot, BfPtrArray const *srcLevelNodes, BfPtrArray const *tgtLevelNodes) {
  BEGIN_ERROR_HANDLING();

  /* neither the source nor target levels should be empty */
  assert(!bfPtrArrayIsEmpty(srcLevelNodes));
  assert(!bfPtrArrayIsEmpty(tgtLevelNodes));

  /* get the layer potential we should use for reexpansion */
  BfLayerPotential proxyLayerPot = BF_PROXY_LAYER_POT[layerPot];

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
  assert(numBlockCols == bfMatBlockGetNumRowBlocks(prevMatBlock));

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
      assert(i < numBlockRows);
      assert(j < numBlockCols);

      /* a priori rank estimate for the original circles */
      BfSize rankOr = bfHelm2RankEstForTwoCircles(
        &srcIter.childCirc, &tgtIter.circ, K, 1, 1e-15);

      /* a priori rank estimate for the new circles */
      BfSize rankEq = bfHelm2RankEstForTwoCircles(
        &srcIter.circ, &tgtIter.childCirc, K, 1, 1e-15);

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
      assert(blockIndex < mat->numBlocks);
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

  BfPoints2 srcChildPts, srcPts, tgtChildPts;

  blockIndex = 0;

  /* iterate over all target child nodes */
  resetMakeFactorIter(&tgtIter, tgtLevelNodes);
  do {
    /* iterate over all source child nodes */
    resetMakeFactorIter(&srcIter, srcLevelNodes);
    do {
      BfSize numRows = bfMatBlockGetNumBlockRows(&mat->super, mat->rowInd[blockIndex]);
      BfSize numCols = bfMatBlockGetNumBlockCols(&mat->super, mat->colInd[blockIndex]);
      assert(numRows > 0 && numCols > 0);

      /* sample points on the source child circle */
      srcChildPts = bfCircle2SamplePoints(&srcIter.childCirc, numCols);
      HANDLE_ERROR();

      /* sample points on the source circle */
      srcPts = bfCircle2SamplePoints(&srcIter.circ, numRows);
      HANDLE_ERROR();

      /* sample points on the target child circle */
      tgtChildPts = bfCircle2SamplePoints(&tgtIter.childCirc, numRows);
      HANDLE_ERROR();

      /* compute the shift matrix for this configuration of circles */
      mat->super.block[blockIndex] = bfHelm2GetReexpansionMatrix(
        &srcChildPts, &srcPts, &tgtChildPts, K, proxyLayerPot);
      HANDLE_ERROR();

      /* if we're debugging, store this block's points---free them
       * otherwise */
#if BF_DEBUG
      BfFacAux *facAux = malloc(sizeof(BfFacAux));
      facAux->srcPts[0] = srcChildPts;
      facAux->srcPts[1] = srcPts;
      facAux->tgtPts = tgtChildPts;
      mat->super.block[blockIndex]->aux = facAux;
#else
      bfFreePoints2(&srcChildPts);
      bfFreePoints2(&srcPts);
      bfFreePoints2(&tgtChildPts);
#endif

      ++blockIndex;
    } while (makeFactorIterNext(&srcIter));
  } while (makeFactorIterNext(&tgtIter));

  END_ERROR_HANDLING() {
    bfFreePoints2(&srcChildPts);
    bfFreePoints2(&srcPts);
    bfFreePoints2(&tgtChildPts);
    bfMatBlockCooDeinitAndDealloc(&mat);
  }

  return bfMatBlockCooToMat(mat);
}

static BfMat *makeLastFactor(BfMat const *prevMat, BfReal K, BfLayerPotential layerPot, BfQuadtree const *tree, BfPtrArray const *srcLevelNodes, BfPtrArray const *tgtLevelNodes) {
  BEGIN_ERROR_HANDLING();

  /* the current level of the target node tree shouldn't be empty */
  assert(!bfPtrArrayIsEmpty(tgtLevelNodes));

  /* the should only be one source node when we make the last
   * butterfly factor */
  assert(bfPtrArraySize(srcLevelNodes) == 1);

  /* the number of diagonal blocks in the last factor equals the
   * number of nodes at the current depth of the target tree ... */
  BfSize numBlocks = bfPtrArraySize(tgtLevelNodes);

  /* ... and this number of blocks should match the number of block
   *  rows of the previous factor */
  assert(numBlocks == prevMat->numRows);

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

  BfPoints2 srcCircPts, tgtPts;
  BfVectors2 tgtNormals;

  /* flag to indicate whether we should fetch the unit normals at the
   * target points for layer potential evaluation */
  bool usingTgtNormals = layerPot != BF_LAYER_POTENTIAL_SINGLE;

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

    /* get the current set of target points */
    tgtPts = bfQuadtreeNodeGetPoints(tgtNode, tree);
    HANDLE_ERROR();

    /* get the current target points' unit normals */
    if (usingTgtNormals)
      tgtNormals = bfQuadtreeNodeGetUnitNormals(tgtNode, tree);

    BfVectors2 *tgtNormalsPtr = usingTgtNormals ? &tgtNormals : NULL;

    mat->super.block[i] = bfGetHelm2KernelMatrix(
      &srcCircPts, &tgtPts, tgtNormalsPtr, K, layerPot);
    HANDLE_ERROR();

    assert(mat->super.rowOffset[i + 1] == BF_SIZE_BAD_VALUE);
    mat->super.rowOffset[i + 1] = mat->super.block[i]->numRows;

    assert(mat->super.colOffset[i + 1] == BF_SIZE_BAD_VALUE);
    mat->super.colOffset[i + 1] = mat->super.block[i]->numCols;

    /* hang onto this block's points if we're in debug mode, and free
     * them otherwise */
#if BF_DEBUG
    BfFacAux *facAux = malloc(sizeof(BfFacAux));
    facAux->srcPts[0] = bfGetUninitializedPoints2();
    facAux->srcPts[1] = srcCircPts;
    facAux->tgtPts = tgtPts;
    mat->super.block[i]->aux = facAux;
#else
    bfFreePoints2(&srcCircPts);
    bfFreePoints2(&tgtPts);
#endif
  }

  /* compute running sum of rows to get row offsets */
  mat->super.rowOffset[0] = 0;
  bfSizeRunningSum(numBlocks + 1, mat->super.rowOffset);

  /* compute running sum of columns to get column offsets */
  mat->super.colOffset[0] = 0;
  bfSizeRunningSum(numBlocks + 1, mat->super.colOffset);

  END_ERROR_HANDLING() {
    bfFreePoints2(&srcCircPts);
    bfFreePoints2(&tgtPts);
    bfMatBlockDiagDeinitAndDealloc(&mat);
  }

  return bfMatBlockDiagToMat(mat);
}

static bool
allRankEstimatesAreOK(BfQuadtreeNode const *tgtNode, BfReal K,
                      BfPtrArray const *srcLevelNodes)
{
  BfCircle tgtCirc = bfQuadtreeNodeGetBoundingCircle(tgtNode);

  for (BfSize i = 0; i < bfPtrArraySize(srcLevelNodes); ++i) {
    BfQuadtreeNode *srcNode = bfPtrArrayGet(srcLevelNodes, i);

    BfSize numSrcPoints = bfTreeNodeGetNumPoints(&srcNode->super);
    BfCircle srcCirc = bfQuadtreeNodeGetBoundingCircle(srcNode);
    BfSize rank = bfHelm2RankEstForTwoCircles(&tgtCirc, &srcCirc, K, 1, 1e-15);

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
bfFacHelm2Prepare(BfQuadtreeNode const *srcNode,
                  BfQuadtreeNode const *tgtNode,
                  BfReal K,
                  BfTreeLevelIter *srcLevelIter,
                  BfTreeLevelIter *tgtLevelIter)
{
  BEGIN_ERROR_HANDLING();

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
  assert(bfTreeLevelIterCurrentLevelIsInternal(tgtLevelIter));
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

  assert(bfTreeLevelIterGetNumPoints(tgtLevelIter) == numTgtNodes);

  /* get the current source and target depths and make sure they're
   * compatible */
  BfSize currentSrcDepth =
    bfTreeLevelIterCurrentDepth(srcLevelIter);
  BfSize currentTgtDepth = bfTreeLevelIterCurrentDepth(tgtLevelIter);
  assert(currentTgtDepth <= currentSrcDepth);

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
  assert(bfTreeLevelIterGetNumPoints(srcLevelIter) == numSrcNodes);

  /* Skip levels until we're on an "internal" level of the tree */
  while (!bfTreeLevelIterCurrentLevelIsInternal(srcLevelIter)) {
    bfTreeLevelIterNext(srcLevelIter);
    --currentSrcDepth;
  }

  /* TODO: step the source level iterator until:
   * - the rank estimate between the first target node and each source
   *   node is smaller than corresponding the number of points */

  while (currentSrcDepth > currentTgtDepth &&
         !allRankEstimatesAreOK(tgtNode, K, &srcLevelIter->levelNodes)) {
    bfTreeLevelIterNext(srcLevelIter);
    --currentSrcDepth;
  }

  /* get number of factors in the butterfly factorization... if we
   * can't butterfly this matrix, return 0 to signal this */
  numFactors = allRankEstimatesAreOK(tgtNode, K, &srcLevelIter->levelNodes) ?
      currentSrcDepth - currentTgtDepth + 2 : 0;

  END_ERROR_HANDLING() {
    bfTreeLevelIterDeinit(srcLevelIter);
    bfTreeLevelIterDeinit(tgtLevelIter);
  }

  return numFactors;
}

BfMatProduct *bfFacHelm2Make(BfQuadtree const *tree, BfReal K, BfLayerPotential layerPot, BfTreeLevelIter *srcLevelIter, BfTreeLevelIter *tgtLevelIter, BfSize numFactors) {
  BEGIN_ERROR_HANDLING();

  /* allocate space for the butterfly factors
   *
   * when multiplying, the factors will be multiplied in the order
   * that they're stored here (e.g., bf_factors[0] is the first factor
   * that will be multiplied) */
  BfMat **factor = malloc(numFactors*sizeof(BfMat *));
  if (factor == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* make the first factor in the butterfly factorization: this is the
   * factor which initially shifts the charges on the source points to
   * the first level of source circles */
  factor[0] = makeFirstFactor(
    K, layerPot, tree, &srcLevelIter->levelNodes, &tgtLevelIter->levelNodes);
  HANDLE_ERROR();

  for (BfSize i = 1; i < numFactors - 1; ++i) {
    /* go up a level on the source tree */
    bfTreeLevelIterNext(srcLevelIter);

    /* make the next factor */
    factor[i] = makeFactor(
      factor[i - 1], K, layerPot, &srcLevelIter->levelNodes,
      &tgtLevelIter->levelNodes);
    HANDLE_ERROR();

    /* go down a level on the target tree */
    bfTreeLevelIterNext(tgtLevelIter);
  }

  /* make the last factor in the butterfly factorization: this is the
   * evaluation factor, which computes the potential at each target
   * point due to the charges on the final source circle */
  factor[numFactors - 1] = makeLastFactor(
    factor[numFactors - 2], K, layerPot, tree, &srcLevelIter->levelNodes,
    &tgtLevelIter->levelNodes);
  HANDLE_ERROR();

  BfMatProduct *prod = bfMatProductNew();
  bfMatProductInit(prod);
  for (BfSize i = 0; i < numFactors; ++i)
    bfMatProductPostMultiply(prod, factor[numFactors - 1 - i]);

  END_ERROR_HANDLING() {
    assert(false); // TODO: need to think carefully about how to do this
  }

  return prod;
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

static BfMat *facHelm2MakeMultilevel_dense(BfQuadtree const *tree, BfReal K, BfLayerPotential layerPot, BfQuadtreeNode const *srcNode, BfQuadtreeNode const *tgtNode) {
  BEGIN_ERROR_HANDLING();

  BfPoints2 srcPts, tgtPts;
  BfVectors2 tgtNormals, *tgtNormalsPtr = NULL;
  bool usingTgtNormals = layerPot != BF_LAYER_POTENTIAL_SINGLE;

  BfMat *Z = NULL;

  srcPts = bfQuadtreeNodeGetPoints(srcNode, tree);
  tgtPts = bfQuadtreeNodeGetPoints(tgtNode, tree);

  if (usingTgtNormals) {
    tgtNormals = bfQuadtreeNodeGetUnitNormals(tgtNode, tree);
    tgtNormalsPtr = &tgtNormals;
  }

  Z = bfGetHelm2KernelMatrix(&srcPts, &tgtPts, tgtNormalsPtr, K, layerPot);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfMatDelete(&Z);

  return Z;
}

static
BfMat *
facHelm2MakeMultilevel_separated(BfQuadtree const *tree, BfReal K,
                                 BfLayerPotential layerPot,
                                 BfQuadtreeNode const *srcNode,
                                 BfQuadtreeNode const *tgtNode) {
  BfTreeLevelIter srcLevelIter, tgtLevelIter;
  BfSize numFactors = bfFacHelm2Prepare(
    srcNode, tgtNode, K, &srcLevelIter, &tgtLevelIter);

  if (numFactors == 0)
    return facHelm2MakeMultilevel_dense(tree, K, layerPot, srcNode, tgtNode);

  BfMatProduct *factorization = bfFacHelm2Make(
    tree, K, layerPot, &srcLevelIter, &tgtLevelIter, numFactors);

  return bfMatProductToMat(factorization);
}

static void facHelm2MakeMultilevel_rec(BfQuadtree const *tree, BfReal K,
                                       BfLayerPotential layerPot,
                                       BfPtrArray const *srcNodes,
                                       BfPtrArray const *tgtNodes,
                                       BfSize level,
                                       BfMatBlockDense *blockMat);

static
BfMat *
facHelm2MakeMultilevel_diag(BfQuadtree const *tree, BfReal K,
                            BfLayerPotential layerPot,
                            BfQuadtreeNode const *srcNode,
                            BfQuadtreeNode const *tgtNode,
                            BfSize level) {
  BEGIN_ERROR_HANDLING();

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

  facHelm2MakeMultilevel_rec(tree, K, layerPot, &srcChildNodes, &tgtChildNodes,
                             level + 1, childBlockMat);
  HANDLE_ERROR();

  assert(bfMatGetNumRows(mat) == bfTreeNodeGetNumPoints(&tgtNode->super));
  assert(bfMatGetNumCols(mat) == bfTreeNodeGetNumPoints(&srcNode->super));

  END_ERROR_HANDLING() {}

  bfPtrArrayDeinit(&srcChildNodes);
  bfPtrArrayDeinit(&tgtChildNodes);

  return bfMatBlockDenseToMat(childBlockMat);
}

static void facHelm2MakeMultilevel_rec(BfQuadtree const *tree, BfReal K,
                                       BfLayerPotential layerPot,
                                       BfPtrArray const *srcNodes,
                                       BfPtrArray const *tgtNodes,
                                       BfSize level,
                                       BfMatBlockDense *blockMat) {
  BEGIN_ERROR_HANDLING();

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
        mat = facHelm2MakeMultilevel_dense(tree, K, layerPot, srcNode, tgtNode);
      else if (separated)
        mat = facHelm2MakeMultilevel_separated(tree, K, layerPot, srcNode, tgtNode);
      else
        /* TODO: we really need to consolidate _rec and _diag (also,
         * "_diag" is a total misnomer) */
        mat = facHelm2MakeMultilevel_diag(tree, K, layerPot, srcNode, tgtNode, level);

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
      assert(bfMatGetNumRows(block) == m);
      assert(bfMatGetNumCols(block) == n);
    }
  }
#endif

  END_ERROR_HANDLING() { /* TODO: ... */ }
}

BfMat *bfFacHelm2MakeMultilevel(BfQuadtree const *tree, BfReal K,
                                BfLayerPotential layerPot) {
  BEGIN_ERROR_HANDLING();

  /* Iterate over the quadtree level by level to generate pairs of
   * nodes to compress */
  BfTreeLevelIter levelIter;
  bfTreeLevelIterInit(
    &levelIter, BF_TREE_TRAVERSAL_LR_LEVEL_ORDER,
    bfQuadtreeNodeToTreeNode((BfQuadtreeNode *)tree->super.root));
  HANDLE_ERROR();

  /* Skip to level 2 of quadtree */
  bfTreeLevelIterNext(&levelIter);
  bfTreeLevelIterNext(&levelIter);
  HANDLE_ERROR();

  /* Get the nodes at the current level along with their count */
  BfPtrArray const *levelNodes = &levelIter.levelNodes;
  BfSize numNodes = bfPtrArraySize(levelNodes);

  /* Create a new dense block matrix to store the HODBF matrix */
  BfMatBlockDense *mat = bfMatBlockDenseNew();
  bfMatBlockDenseInit(mat, numNodes, numNodes);

  /* Build the multilevel butterfly factorization */
  facHelm2MakeMultilevel_rec(tree, K, layerPot, levelNodes, levelNodes, 2, mat);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfMatBlockDenseDeinitAndDealloc(&mat);

  bfTreeLevelIterDeinit(&levelIter);

  return bfMatBlockDenseToMat(mat);
}
