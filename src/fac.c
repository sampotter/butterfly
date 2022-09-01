#include <bf/fac.h>

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include <bf/error_macros.h>
#include <bf/helm2.h>
#include <bf/mat_block.h>
#include <bf/mat_block_coo.h>
#include <bf/mat_block_dense.h>
#include <bf/mat_block_diag.h>
#include <bf/mat_product.h>
#include <bf/util.h>

#define MAX_DENSE_MATRIX_SIZE 16384 // == 128*128

static BfSize
getChildren(BfQuadtreeNode const *node, BfQuadtreeNode const *child[4]) {
  for (BfSize i = 0; i < 4; ++i)
    child[i] = NULL;

  BfSize numChildren = 0;
  for (BfSize i = 0; i < 4; ++i)
    if (node->child[i] != NULL)
      child[numChildren++] = node->child[i];
  return numChildren;
}

static BfSize
getNumChildren(BfQuadtreeNode const *node) {
  BfSize numChildren = 0;
  for (BfSize i = 0; i < 4; ++i)
    numChildren += node->child[i] != NULL;
  return numChildren;
}

static BfMatBlockDiag *makeFirstFactor(BfReal K,
                                       BfQuadtree const *tree,
                                       BfPtrArray const *srcLevelNodes,
                                       BfPtrArray const *tgtLevelNodes)
{
  BEGIN_ERROR_HANDLING();

  /* there should only be one target node at the starting level of
   * the target tree when we begin the traversal */
  assert(bfPtrArraySize(tgtLevelNodes) == 1);

  /* the current level of the source node tree shouldn't be empty */
  assert(!bfPtrArrayIsEmpty(srcLevelNodes));

  BfQuadtreeNode const *srcNode = NULL;
  BfQuadtreeNode const *tgtNode = NULL;

  /* get the lone target node and its bounding circle */
  bfPtrArrayGetFirst(tgtLevelNodes, (BfPtr *)&tgtNode);
  BfCircle2 tgtCirc = bfGetQuadtreeNodeBoundingCircle(tgtNode);

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
    BfCircle2 srcCirc = bfGetQuadtreeNodeBoundingCircle(srcNode);

    /* get the original points in the current source node box */
    bfGetQuadtreeNodePoints(tree, srcNode, &srcPts);
    HANDLE_ERROR();

    /* verify that the source bounding circle contains the points */
    assert(bfCircle2ContainsPoints(&srcCirc, &srcPts));

    /* get the rank estimate for the current pair of source and
`     * target bounding circles */
    BfSize p = bfHelm2RankEstForTwoCircles(&srcCirc, &tgtCirc, K, 1, 1e-15);

    /* sample points on the source circle */
    tgtCircPts = bfSamplePointsOnCircle2(&tgtCirc, p);
    HANDLE_ERROR();

    /* sample points on the target circle */
    srcCircPts = bfSamplePointsOnCircle2(&srcCirc, p);
    HANDLE_ERROR();

    /* compute the shift matrix and store it in the current block */
    mat->super.block[i] = bfMatDenseComplexGetMatPtr(
      bfHelm2GetShiftMatrix(&srcPts, &srcCircPts, &tgtCircPts, K));
    HANDLE_ERROR();

    /* continue initializing the row and column offsets */
    mat->super.rowOffset[i + 1] = bfMatGetNumRows(mat->super.block[i]);
    mat->super.colOffset[i + 1] = bfMatGetNumCols(mat->super.block[i]);

    /* if we're debugging, store this block's points---free them
     * otherwise */
#if BF_DEBUG
    BfPoints2 *auxPts = malloc(3*sizeof(BfPoints2));
    auxPts[0] = srcPts;
    auxPts[1] = srcCircPts;
    auxPts[2] = tgtCircPts;
    mat->super.block[i]->aux = auxPts;
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
    bfMatBlockDiagDeinitAndDelete(&mat);
  }

  /* finish initializing row and column offsets */
  mat->super.rowOffset[0] = mat->super.colOffset[0] = 0;
  bfSizeRunningSum(numBlocks + 1, mat->super.rowOffset);
  bfSizeRunningSum(numBlocks + 1, mat->super.colOffset);

  return mat;
}

static BfSize getTotalNumChildren(BfPtrArray const *levelNodes) {
  BfSize totalNumChildren = 0;
  for (BfSize i = 0; i < bfPtrArraySize(levelNodes); ++i) {
    BfQuadtreeNode const *node = bfPtrArrayGet(levelNodes, i);
    totalNumChildren += getNumChildren(node);
  }
  return totalNumChildren;
}

typedef struct {
  BfPtrArray const *levelNodes;
  BfSize nodeIndex;
  BfSize childIndex;
  BfSize numChildren;
  BfQuadtreeNode const *node;
  BfQuadtreeNode const *children[4];
  BfQuadtreeNode const *child;
  BfCircle2 circ;
  BfCircle2 childCirc;
} MakeFactorIter;

static void resetMakeFactorIter(MakeFactorIter *iter, BfPtrArray const *levelNodes) {
  iter->levelNodes = levelNodes;

  iter->nodeIndex = iter->childIndex = 0;

  iter->node = bfPtrArrayGet(iter->levelNodes, iter->nodeIndex);
  iter->numChildren = getChildren(iter->node, iter->children);
  iter->child = iter->children[iter->childIndex];

  iter->circ = bfGetQuadtreeNodeBoundingCircle(iter->node);
  iter->childCirc = bfGetQuadtreeNodeBoundingCircle(iter->child);
}

static bool makeFactorIterNext(MakeFactorIter *iter) {
  if (++iter->childIndex == iter->numChildren) {
    if (++iter->nodeIndex == bfPtrArraySize(iter->levelNodes))
      return false;

    iter->childIndex = 0;

    iter->node = bfPtrArrayGet(iter->levelNodes, iter->nodeIndex);

    iter->numChildren = getChildren(iter->node, iter->children);

    iter->circ = bfGetQuadtreeNodeBoundingCircle(iter->node);
  }

  assert(iter->childIndex < iter->numChildren);
  iter->child = iter->children[iter->childIndex];

  iter->childCirc = bfGetQuadtreeNodeBoundingCircle(iter->child);

  return true;
}

/* This function computes an inner factor in a butterfly
 * factorization.
 *
 * ... explain how this works ...
 */
static BfMatBlockCoo *
makeFactor(BfMatBlock const *prevMat, BfReal K,
           BfPtrArray const *srcLevelNodes, BfPtrArray const *tgtLevelNodes)
{
  BEGIN_ERROR_HANDLING();

  /* neither the source nor target levels should be empty */
  assert(!bfPtrArrayIsEmpty(srcLevelNodes));
  assert(!bfPtrArrayIsEmpty(tgtLevelNodes));

  /* count the total number of target children on this level */
  BfSize totalNumTgtChildren = getTotalNumChildren(tgtLevelNodes);

  /* count the total number of source children on this level */
  BfSize totalNumSrcChildren = getTotalNumChildren(srcLevelNodes);

  /* compute the number of block rows and columns as well as the total
   * number of blocks from this level's layout */
  BfSize numBlocks = totalNumSrcChildren*totalNumTgtChildren;
  BfSize numBlockRows = totalNumTgtChildren*bfPtrArraySize(srcLevelNodes);
  BfSize numBlockCols = prevMat->super.numRows;

  BfMatBlockCoo *mat = bfMatBlockCooNew();
  bfMatBlockCooInit(mat, numBlockRows, numBlockCols, numBlocks);
  HANDLE_ERROR();

  /* set the number of columns in each block column to equal the
   * number of rows in each block row of the previous factor  */
  for (BfSize j = 0; j < numBlockCols; ++j)
    mat->super.colOffset[j + 1] = bfMatBlockGetNumBlockRows(prevMat, j);

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
    j = tgtIter.nodeIndex*srcIter.numChildren;
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
      srcChildPts = bfSamplePointsOnCircle2(&srcIter.childCirc, numCols);
      HANDLE_ERROR();

      /* sample points on the source circle */
      srcPts = bfSamplePointsOnCircle2(&srcIter.circ, numRows);
      HANDLE_ERROR();

      /* sample points on the target child circle */
      tgtChildPts = bfSamplePointsOnCircle2(&tgtIter.childCirc, numRows);
      HANDLE_ERROR();

      /* compute the shift matrix for this configuration of circles */
      mat->super.block[blockIndex] =
        (BfMat *)bfHelm2GetShiftMatrix(&srcChildPts, &srcPts, &tgtChildPts, K);
      HANDLE_ERROR();

      /* if we're debugging, store this block's points---free them
       * otherwise */
#if BF_DEBUG
      BfPoints2 *auxPts = malloc(3*sizeof(BfPoints2));
      auxPts[0] = srcChildPts;
      auxPts[1] = srcPts;
      auxPts[2] = tgtChildPts;
      mat->super.block[blockIndex]->aux = auxPts;
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
    bfMatBlockCooDeinitAndDelete(&mat);
  }

  return mat;
}

static BfMatBlockDiag *makeLastFactor(BfMatBlock const *prevMat, BfReal K,
                                      BfQuadtree const *tree,
                                      BfPtrArray const *srcLevelNodes,
                                      BfPtrArray const *tgtLevelNodes)
{
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
  assert(numBlocks == prevMat->super.numRows);

  BfQuadtreeNode const *srcNode = NULL;
  BfQuadtreeNode const *tgtNode = NULL;

  /* get the lone source node and its bounding circle */
  bfPtrArrayGetFirst(srcLevelNodes, (BfPtr *)&srcNode);
  BfCircle2 srcCirc = bfGetQuadtreeNodeBoundingCircle(srcNode);

  /* initialize the last butterfly factor */
  BfMatBlockDiag *mat = bfMatBlockDiagNew();
  bfMatBlockDiagInit(mat, numBlocks, numBlocks);
  HANDLE_ERROR();

  BfPoints2 srcCircPts, tgtPts;

  /* iterate over each node of the final level of the target node tree
   * and compute the matrix which will evaluate the potential at each
   * target point due to the charges on the source proxy points */
  for (BfSize i = 0; i < numBlocks; ++i) {
    BfSize prevNumRows = bfMatBlockGetNumBlockRows(prevMat, i);

    /* get the current target node */
    tgtNode = bfPtrArrayGet(tgtLevelNodes, i);

    /* get the proxy points on the source circle for the current
     * source block */
    srcCircPts = bfSamplePointsOnCircle2(&srcCirc, prevNumRows);
    HANDLE_ERROR();

    /* get the current set of target points */
    bfGetQuadtreeNodePoints(tree, tgtNode, &tgtPts);
    HANDLE_ERROR();

    mat->super.block[i] = (BfMat *)bfGetHelm2KernelMatrix(&srcCircPts, &tgtPts, K);
    HANDLE_ERROR();

    assert(mat->super.rowOffset[i + 1] == BF_SIZE_BAD_VALUE);
    mat->super.rowOffset[i + 1] = mat->super.block[i]->numRows;

    assert(mat->super.colOffset[i + 1] == BF_SIZE_BAD_VALUE);
    mat->super.colOffset[i + 1] = mat->super.block[i]->numCols;

    /* hang onto this block's points if we're in debug mode, and free
     * them otherwise */
#if BF_DEBUG
    BfPoints2 *auxPts = malloc(3*sizeof(BfPoints2));
    auxPts[0] = bfGetUninitializedPoints2();
    auxPts[1] = srcCircPts;
    auxPts[2] = tgtPts;
    mat->super.block[i]->aux = auxPts;
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
    bfMatBlockDiagDeinitAndDelete(&mat);
  }

  return mat;
}

static bool
allSourceNodesHaveChildren(BfPtrArray const *srcLevelNodes)
{
  for (BfSize i = 0; i < bfPtrArraySize(srcLevelNodes); ++i) {
    BfQuadtreeNode *srcNode = bfPtrArrayGet(srcLevelNodes, i);
    if (bfQuadtreeNodeIsLeaf(srcNode))
      return false;
  }
  return true;
}

static bool
allRankEstimatesAreOK(BfQuadtreeNode const *tgtNode, BfReal K,
                      BfPtrArray const *srcLevelNodes)
{
  BfCircle2 tgtCirc = bfGetQuadtreeNodeBoundingCircle(tgtNode);

  for (BfSize i = 0; i < bfPtrArraySize(srcLevelNodes); ++i) {
    BfQuadtreeNode *srcNode = bfPtrArrayGet(srcLevelNodes, i);

    BfSize numSrcPoints = bfQuadtreeNodeNumPoints(srcNode);
    BfCircle2 srcCirc = bfGetQuadtreeNodeBoundingCircle(srcNode);
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
 * in false negatives.) */
BfSize
bfFacHelm2Prepare(BfQuadtree const *tree, BfQuadtreeNode const *srcNode,
                  BfQuadtreeNode const *tgtNode, BfReal K,
                  BfQuadtreeLevelIter *srcLevelIter,
                  BfQuadtreeLevelIter *tgtLevelIter)
{
  BEGIN_ERROR_HANDLING();

  BfSize numFactors = 0;

  /* set up the level iterator for the source tree---this iterator
   * goes from the leaves of the tree to the root (in reverse) */
  *srcLevelIter = bfInitQuadtreeLevelIter(
    BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER, (BfQuadtreeNode *)srcNode);
  HANDLE_ERROR();

  /* set up the level iterator for the target tree, which goes from
   * the root to the leaves */
  *tgtLevelIter = bfInitQuadtreeLevelIter(
    BF_TREE_TRAVERSAL_LR_LEVEL_ORDER, (BfQuadtreeNode *)tgtNode);
  HANDLE_ERROR();

  /* get the current source and target depths and make sure they're
   * compatible */
  BfSize currentSrcDepth = bfQuadtreeLevelIterCurrentDepth(srcLevelIter);
  BfSize currentTgtDepth = bfQuadtreeLevelIterCurrentDepth(tgtLevelIter);
  assert(currentTgtDepth <= currentSrcDepth);

  /* skip source levels until we're no deeper than the maximum depth
   * beneath the target node */
  BfSize maxDepthBelowTgtNode = bfGetMaxDepthBelowQuadtreeNode(tgtNode);
  while (currentSrcDepth > maxDepthBelowTgtNode) {
    bfQuadtreeLevelIterNext(srcLevelIter);
    --currentSrcDepth;
  }

  /* skip source levels until we reach a level where each node has children
   *
   * TODO: I am *NOT AT ALL* sure this is right... but it should at
   * get me unstuck for now at least. */
  while (currentSrcDepth > currentTgtDepth &&
         !allSourceNodesHaveChildren(&srcLevelIter->levelNodes)) {
    bfQuadtreeLevelIterNext(srcLevelIter);
    --currentSrcDepth;
  }
  assert(allSourceNodesHaveChildren(&srcLevelIter->levelNodes));

  /* TODO: step the source level iterator until:
   * - the rank estimate between the first target node and each source
   *   node is smaller than corresponding the number of points */

  while (currentSrcDepth > currentTgtDepth &&
         !allRankEstimatesAreOK(tgtNode, K, &srcLevelIter->levelNodes)) {
    bfQuadtreeLevelIterNext(srcLevelIter);
    --currentSrcDepth;
  }

  /* get number of factors in the butterfly factorization... if we
   * can't butterfly this matrix, return 0 to signal this */
  numFactors = allRankEstimatesAreOK(tgtNode, K, &srcLevelIter->levelNodes) ?
      currentSrcDepth - currentTgtDepth + 2 : 0;

  END_ERROR_HANDLING() {
    bfFreeQuadtreeLevelIter(srcLevelIter);
    bfFreeQuadtreeLevelIter(tgtLevelIter);
  }

  return numFactors;
}

BfMatProduct *
bfFacHelm2Make(BfQuadtree const *tree, BfQuadtreeNode const *srcNode,
               BfQuadtreeNode const *tgtNode, BfReal K,
               BfQuadtreeLevelIter *srcLevelIter,
               BfQuadtreeLevelIter *tgtLevelIter,
               BfSize numFactors)
{
  BEGIN_ERROR_HANDLING();

  BfMatBlock *prevMatBlock = NULL;
  BfMatBlockCoo *matBlockCoo = NULL;
  BfMatBlockDiag *matBlockDiag = NULL;

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
  matBlockDiag = makeFirstFactor(
    K, tree, &srcLevelIter->levelNodes, &tgtLevelIter->levelNodes);
  factor[0] = bfMatBlockDiagGetMatPtr(matBlockDiag);
  HANDLE_ERROR();

  for (BfSize i = 1; i < numFactors - 1; ++i) {
    /* go up a level on the source tree */
    bfQuadtreeLevelIterNext(srcLevelIter);

    /* make the next factor */
    prevMatBlock = (BfMatBlock *)factor[i - 1];
    matBlockCoo = makeFactor(
      prevMatBlock, K, &srcLevelIter->levelNodes, &tgtLevelIter->levelNodes);
    factor[i] = bfMatBlockCooGetMatPtr(matBlockCoo);
    HANDLE_ERROR();

    /* go down a level on the target tree */
    bfQuadtreeLevelIterNext(tgtLevelIter);
  }

  /* make the last factor in the butterfly factorization: this is the
   * evaluation factor, which computes the potential at each target
   * point due to the charges on the final source circle */
  prevMatBlock = (BfMatBlock *)factor[numFactors - 2];
  matBlockDiag = makeLastFactor(
    prevMatBlock, K, tree, &srcLevelIter->levelNodes, &tgtLevelIter->levelNodes);
  factor[numFactors - 1] = bfMatBlockDiagGetMatPtr(matBlockDiag);
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
    if (node->child[i])
      bfPtrArrayAppend(&childNodes, node->child[i]);
  return childNodes;
}

static
BfMat *
facHelm2MakeMultilevel_dense(BfQuadtree const *tree, BfReal K,
                             BfQuadtreeNode const *srcNode,
                             BfQuadtreeNode const *tgtNode) {
  BEGIN_ERROR_HANDLING();

  BfPoints2 srcPts, tgtPts;
  BfMatDenseComplex *Z = NULL;

  bfGetQuadtreeNodePoints(tree, srcNode, &srcPts);
  bfGetQuadtreeNodePoints(tree, tgtNode, &tgtPts);

  Z = bfGetHelm2KernelMatrix(&srcPts, &tgtPts, K);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDelete(&Z);

  return bfMatDenseComplexGetMatPtr(Z);
}

static
BfMat *
facHelm2MakeMultilevel_separated(BfQuadtree const *tree, BfReal K,
                                 BfQuadtreeNode const *srcNode,
                                 BfQuadtreeNode const *tgtNode) {
  BfQuadtreeLevelIter srcLevelIter, tgtLevelIter;
  BfSize numFactors = bfFacHelm2Prepare(
    tree, srcNode, tgtNode, K, &srcLevelIter, &tgtLevelIter);

  if (numFactors == 0)
    return facHelm2MakeMultilevel_dense(tree, K, srcNode, tgtNode);

  BfMatProduct *factorization = bfFacHelm2Make(
    tree, srcNode, tgtNode, K, &srcLevelIter, &tgtLevelIter, numFactors);

  return bfMatProductGetMatPtr(factorization);
}

static void facHelm2MakeMultilevel_rec(BfQuadtree const *tree, BfReal K,
                                       BfPtrArray const *srcNodes,
                                       BfPtrArray const *tgtNodes,
                                       BfSize level,
                                       BfMatBlockDense *blockMat);

static
BfMat *
facHelm2MakeMultilevel_diag(BfQuadtree const *tree, BfReal K,
                            BfQuadtreeNode const *srcNode,
                            BfQuadtreeNode const *tgtNode,
                            BfSize level) {
  BEGIN_ERROR_HANDLING();

  BfPtrArray srcChildNodes, tgtChildNodes;
  BfMatBlockDense *childBlockMat = NULL;

  srcChildNodes = getChildrenAsPtrArray(srcNode);
  HANDLE_ERROR();

  tgtChildNodes = getChildrenAsPtrArray(tgtNode);
  HANDLE_ERROR();

  childBlockMat = bfMatBlockDenseNew();
  HANDLE_ERROR();

  BfSize numBlockRows = bfPtrArraySize(&tgtChildNodes);
  BfSize numBlockCols = bfPtrArraySize(&srcChildNodes);
  bfMatBlockDenseInit(childBlockMat, numBlockRows, numBlockCols);
  HANDLE_ERROR();

  facHelm2MakeMultilevel_rec(
    tree, K, &srcChildNodes, &tgtChildNodes, level + 1, childBlockMat);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {}

  bfFreePtrArray(&srcChildNodes);
  bfFreePtrArray(&tgtChildNodes);

  return bfMatBlockDenseGetMatPtr(childBlockMat);
}

static void facHelm2MakeMultilevel_rec(BfQuadtree const *tree, BfReal K,
                                       BfPtrArray const *srcNodes,
                                       BfPtrArray const *tgtNodes,
                                       BfSize level,
                                       BfMatBlockDense *blockMat) {
  BEGIN_ERROR_HANDLING();

  BfMatBlock *super = &blockMat->super;
  BfMat *mat = NULL;

  super->rowOffset[0] = 0;
  super->colOffset[0] = 0;

  for (BfSize i = 0; i < bfPtrArraySize(tgtNodes); ++i) {
    BfQuadtreeNode *tgtNode = bfPtrArrayGet(tgtNodes, i);
    BfSize numRows = bfQuadtreeNodeNumPoints(tgtNode);

    for (BfSize j = 0; j < bfPtrArraySize(srcNodes); ++j) {
      BfQuadtreeNode *srcNode = bfPtrArrayGet(srcNodes, j);
      BfSize numCols = bfQuadtreeNodeNumPoints(srcNode);

      bool separated = bfQuadtreeNodesAreSeparated(srcNode, tgtNode);

      if (numRows*numCols < MAX_DENSE_MATRIX_SIZE)
        mat = facHelm2MakeMultilevel_dense(tree, K, srcNode, tgtNode);
      else if (separated)
        mat = facHelm2MakeMultilevel_separated(tree, K, srcNode, tgtNode);
      else
        /* TODO: we really need to consolidate _rec and _diag (also,
         * "_diag" is a total misnomer) */
        mat = facHelm2MakeMultilevel_diag(tree, K, srcNode, tgtNode, level);

      if (mat == NULL)
        RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

      /* set the current block now that we've computed it */
      bfMatBlockDenseSetBlock(blockMat, i, j, mat);
      HANDLE_ERROR();

      /* update the column block offset if we're making our first pass
       * through a row of blocks */
      if (i == 0)
        super->colOffset[j + 1] = super->colOffset[j] + numCols;
    }

    /* update the row block offset */
    super->rowOffset[i + 1] = super->rowOffset[i] + numRows;
  }

  END_ERROR_HANDLING() {}
}

BfMatBlockDense *bfFacHelm2MakeMultilevel(BfQuadtree const *tree, BfReal K) {
  BEGIN_ERROR_HANDLING();

  /* Use an in-order quadtree level iterator to enumerate the pairs of
   * nodes which we'll make butterfly factorizations for */
  BfQuadtreeLevelIter levelIter = bfInitQuadtreeLevelIter(
    BF_TREE_TRAVERSAL_LR_LEVEL_ORDER, (BfQuadtreeNode *)tree->root);
  HANDLE_ERROR();

  /* Skip to level 2 of quadtree */
  bfQuadtreeLevelIterNext(&levelIter);
  bfQuadtreeLevelIterNext(&levelIter);
  HANDLE_ERROR();

  BfPtrArray const *levelNodes = &levelIter.levelNodes;
  BfSize numNodes = bfPtrArraySize(levelNodes);

  BfMatBlockDense *mat = bfMatBlockDenseNew();
  bfMatBlockDenseInit(mat, numNodes, numNodes);

  /* Build the multilevel butterfly factorization */
  facHelm2MakeMultilevel_rec(tree, K, levelNodes, levelNodes, 2, mat);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfMatBlockDenseDeinitAndDelete(&mat);

  bfFreeQuadtreeLevelIter(&levelIter);

  return mat;
}
