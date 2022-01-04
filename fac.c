#include "fac.h"

#include <assert.h>
#include <math.h>

#include "error_macros.h"
#include "helm2.h"

static void
initEmptyFactor(BfFactor *factor, BfSize numBlockRows, BfSize numBlockCols,
                BfSize numBlocks)
{
  BEGIN_ERROR_HANDLING();

  factor->numBlockRows = numBlockRows;
  factor->numBlockCols = numBlockCols;
  factor->numBlocks = numBlocks;

  /* allocate and initialize row indices to dummy values */
  factor->rowInd = malloc(numBlocks*sizeof(BfSize));
  if (factor->rowInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  for (BfSize i = 0; i < numBlocks; ++i)
    factor->rowInd[i] = BF_SIZE_BAD_VALUE;

  /* allocate and initialize column indices to dummy values */
  factor->colInd = malloc(numBlocks*sizeof(BfSize));
  if (factor->colInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  for (BfSize i = 0; i < numBlocks; ++i)
    factor->colInd[i] = BF_SIZE_BAD_VALUE;

  /* alloc and init the row index offsets for each block row */
  factor->rowOffset = malloc((numBlockRows + 1)*sizeof(BfSize));
  if (factor->rowOffset == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  for (BfSize i = 0; i < numBlockRows + 1; ++i)
    factor->rowOffset[i] = BF_SIZE_BAD_VALUE;

  /* alloc and init the column index offsets for each block column */
  factor->colOffset = malloc((numBlockCols + 1)*sizeof(BfSize));
  if (factor->colOffset == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  for (BfSize i = 0; i < numBlockCols + 1; ++i)
    factor->colOffset[i] = BF_SIZE_BAD_VALUE;

  /* allocate and initialize each nonzero (diagonal) block */
  factor->block = malloc(numBlocks*sizeof(BfMat));
  if (factor->block == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  for (BfSize i = 0; i < numBlocks; ++i)
    factor->block[i] = bfGetUninitializedMat();

  /* if we're debugging, allocate space in this factor to store each
   * block's points */
#if BF_DEBUG
  factor->srcPtsOrig = malloc(numBlocks*sizeof(BfPoints2));
  if (factor->srcPtsOrig == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  for (BfSize i = 0; i < numBlocks; ++i)
    factor->srcPtsOrig[i] = bfGetUninitializedPoints2();

  factor->srcPtsEquiv = malloc(numBlocks*sizeof(BfPoints2));
  if (factor->srcPtsEquiv == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  for (BfSize i = 0; i < numBlocks; ++i)
    factor->srcPtsEquiv[i] = bfGetUninitializedPoints2();

  factor->tgtPts = malloc(numBlocks*sizeof(BfPoints2));
  if (factor->tgtPts == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  for (BfSize i = 0; i < numBlocks; ++i)
    factor->tgtPts[i] = bfGetUninitializedPoints2();
#endif

  END_ERROR_HANDLING() {
    bfFreeFactor(factor);
  }
}

static void initDiagonalFactor(BfFactor *factor, BfSize numBlocks) {
  BEGIN_ERROR_HANDLING();

  initEmptyFactor(factor, numBlocks, numBlocks, numBlocks);
  HANDLE_ERROR();

  for (BfSize i = 0; i < numBlocks; ++i)
    factor->rowInd[i] = i;

  for (BfSize i = 0; i < numBlocks; ++i)
    factor->colInd[i] = i;

  END_ERROR_HANDLING() {
    bfFreeFactor(factor);
  }
}

void bfFreeFactor(BfFactor *factor) {
  free(factor->rowInd);
  free(factor->colInd);
  free(factor->rowOffset);
  free(factor->colOffset);

  for (BfSize i = 0; i < factor->numBlocks; ++i)
    bfFreeMat(&factor->block[i]);
  free(factor->block);

#if BF_DEBUG
  for (BfSize i = 0; i < factor->numBlocks; ++i)
    bfFreePoints2(&factor->srcPtsOrig[i]);
  free(factor->srcPtsOrig);

  /* note that the final factor doesn't have equivalent source points,
   * so we don't need to free them here */
  for (BfSize i = 0; i < factor->numBlocks; ++i)
    if (bfPoints2Initialized(&factor->srcPtsEquiv[i]))
      bfFreePoints2(&factor->srcPtsEquiv[i]);
  free(factor->srcPtsEquiv);

  for (BfSize i = 0; i < factor->numBlocks; ++i)
    bfFreePoints2(&factor->tgtPts[i]);
  free(factor->tgtPts);
#endif
}

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

static void
getShiftMat(BfPoints2 const *srcPtsOrig, BfPoints2 const *srcPtsEquiv,
            BfPoints2 const *tgtPts, BfReal K, BfMat *Z_shift)
{
  BEGIN_ERROR_HANDLING();

  /* compute the kernel matrix mapping charges on the original sources
   * points to potentials on the original target points */
  BfMat Z_orig = bfGetUninitializedMat();
  bfGetHelm2KernelMatrix(&Z_orig, srcPtsOrig, tgtPts, K);
  HANDLE_ERROR();

  /* compute the kernel matrix mapping charges on the source
   * circle to potentials on the target circle */
  BfMat Z_equiv = bfGetUninitializedMat();
  bfGetHelm2KernelMatrix(&Z_equiv, srcPtsEquiv, tgtPts, K);
  HANDLE_ERROR();

  /* set the "shift matrix" to Z_equiv\Z_orig */
  *Z_shift = bfGetUninitializedMat();
  bfMatLstSq(&Z_equiv, &Z_orig, Z_shift);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    bfFreeMat(Z_shift);
  }

  bfFreeMat(&Z_orig);
  bfFreeMat(&Z_equiv);
}

static BfSize getRows(BfFactor const *factor, BfSize i) {
  assert(i < factor->numBlockRows);

  return factor->rowOffset[i + 1] - factor->rowOffset[i];
}

static BfSize getCols(BfFactor const *factor, BfSize j) {
  assert(j < factor->numBlockCols);

  return factor->colOffset[j + 1] - factor->colOffset[j];
}

static BfSize getNumRows(BfFactor const *factor) {
  return factor->rowOffset[factor->numBlockRows];
}

static BfSize getNumCols(BfFactor const *factor) {
  return factor->colOffset[factor->numBlockCols];
}

/* When initializing `factor`, we first set `factor->rowOffset[i + 1]`
 * and `factor->colOffset[j + 1]` for each `i` and `j` to contain the
 * number of rows and columns in the `i`th block row and `j`th block
 * column, respectively. Afterwards, to transform these arrays into
 * the actual block row and column offsets, we set the first value of
 * each array to zero and take their running sums in place. */
static void cumSumRowAndColOffsets(BfFactor *factor) {
  /* replace `rowOffset` with its cumulative sum */
  factor->rowOffset[0] = 0;
  for (BfSize i = 0; i < factor->numBlockRows; ++i) {
    assert(factor->rowOffset[i + 1] > 0);
    factor->rowOffset[i + 1] += factor->rowOffset[i];
  }

  /* replace `colOffset` with its cumulative sum */
  factor->colOffset[0] = 0;
  for (BfSize i = 0; i < factor->numBlockCols; ++i) {
    assert(factor->colOffset[i + 1] > 0);
    factor->colOffset[i + 1] += factor->colOffset[i];
  }
}

static void
makeFirstFactor(BfFactor *factor, BfReal K,
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
  initDiagonalFactor(factor, numBlocks);
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
    bfPtrArrayGet(srcLevelNodes, i, (BfPtr *)&srcNode);
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
    getShiftMat(&srcPts, &srcCircPts, &tgtCircPts, K, &factor->block[i]);
    HANDLE_ERROR();

    /* set the next row offset to current block's number of rows */
    assert(factor->block[i].numRows > 0);
    factor->rowOffset[i + 1] = factor->block[i].numRows;

    /* set the next column offset to current block's number of column */
    assert(factor->block[i].numCols > 0);
    factor->colOffset[i + 1] = factor->block[i].numCols;

    /* if we're debugging, store this block's points---free them
     * otherwise */
#if BF_DEBUG
    factor->srcPtsOrig[i] = srcPts;
    factor->srcPtsEquiv[i] = srcCircPts;
    factor->tgtPts[i] = tgtCircPts;
#else
    bfFreePoints2(&srcPts);
    bfFreePoints2(&tgtCircPts);
    bfFreePoints2(&srcCircPts);
#endif
  }

  cumSumRowAndColOffsets(factor);

  END_ERROR_HANDLING() {
    bfFreePoints2(&srcPts);
    bfFreePoints2(&tgtCircPts);
    bfFreePoints2(&srcCircPts);
  }
}

static BfSize getTotalNumChildren(BfPtrArray const *levelNodes) {
  BfSize totalNumChildren = 0;
  for (BfSize i = 0; i < bfPtrArraySize(levelNodes); ++i) {
    BfQuadtreeNode const *node;
    bfPtrArrayGet(levelNodes, i, (BfPtr *)&node);
    totalNumChildren += getNumChildren(node);
  }
  return totalNumChildren;
}

/* This function computes an inner factor in a butterfly
 * factorization.
 *
 * ... explain how this works ...
 */
static void
makeFactor(BfFactor *factor, BfFactor const *prevFactor, BfReal K,
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
  BfSize numBlockCols = prevFactor->numBlockRows;

  initEmptyFactor(factor, numBlockRows, numBlockCols, numBlocks);
  HANDLE_ERROR();

  for (BfSize i = 1; i <= numBlockRows; ++i)
    factor->rowOffset[i] = 0;

  /* set the number of columns in each block column to equal the
   * number of rows in each block row of the previous factor */
  factor->colOffset[0] = 0;
  for (BfSize j = 0; j < numBlockCols; ++j)
    factor->colOffset[j + 1] = getRows(prevFactor, j);

  /* bookkeeping variables for the next two sections */
  BfSize p0, p1, q0, q1;
  BfSize i_offset, j;
  BfSize numTgtChildren, numSrcChildren;
  BfQuadtreeNode const *tgtNode, *tgtChild[4];
  BfQuadtreeNode const *srcNode, *srcChild[4];
  BfCircle2 tgtCirc, tgtChildCirc;
  BfCircle2 srcCirc, srcChildCirc;

  /* next, we set the number of block rows and block columns
   *
   * we need to do this as an intermediate step since the individual
   * rank estimates for each pair of source and target circles could
   * vary slightly---we need the number of rows and columns of each
   * block to be compatible. hence, we use the corresponding maximum
   * rank estimate. */

  BfSize blockIndex = 0;

  i_offset = 0;

  /* iterate over all target children */
  for (p0 = 0; p0 < bfPtrArraySize(tgtLevelNodes); ++p0) {
    bfPtrArrayGet(tgtLevelNodes, p0, (BfPtr *)&tgtNode);
    numTgtChildren = getChildren(tgtNode, tgtChild);
    tgtCirc = bfGetQuadtreeNodeBoundingCircle(tgtNode);
    for (p1 = 0; p1 < numTgtChildren; ++p1) {
      tgtChildCirc = bfGetQuadtreeNodeBoundingCircle(tgtChild[p1]);

      j = p0*numSrcChildren;

      /* iterate over all source children */
      for (q0 = 0; q0 < bfPtrArraySize(srcLevelNodes); ++q0) {
        bfPtrArrayGet(srcLevelNodes, q0, (BfPtr *)&srcNode);
        numSrcChildren = getChildren(srcNode, srcChild);
        srcCirc = bfGetQuadtreeNodeBoundingCircle(srcNode);
        for (q1 = 0; q1 < numSrcChildren; ++q1) {
          srcChildCirc = bfGetQuadtreeNodeBoundingCircle(srcChild[q1]);

          BfSize i = i_offset + q0;

          assert(i < numBlockRows);
          assert(j < numBlockCols);

          /* a priori rank estimate for the original circles */
          BfSize rankOr = bfHelm2RankEstForTwoCircles(
            &srcChildCirc, &tgtCirc, K, 1, 1e-15);

          /* a priori rank estimate for the new circles */
          BfSize rankEq = bfHelm2RankEstForTwoCircles(
            &srcCirc, &tgtChildCirc, K, 1, 1e-15);

          /* use the larger of the two rank estimates... not sure if
           * this is totally necessary, probably being a little
           * paranoid... they should be nearly the same, but might
           * differ a little due to rounding */
          BfSize rank = rankOr > rankEq ? rankOr : rankEq;

          /* update number of rows for current block */
          if (rank > factor->rowOffset[i + 1])
            factor->rowOffset[i + 1] = rank;

          /* set block row and column indices */
          assert(blockIndex < factor->numBlocks);
          factor->rowInd[blockIndex] = i;
          factor->colInd[blockIndex] = j;

          ++blockIndex;

          ++j;

        }
      } /* end loop over all source children */

      i_offset += bfPtrArraySize(srcLevelNodes);

    }
  } /* end loop over all target children */

  cumSumRowAndColOffsets(factor);

  // foreach tgtNode:
  //     foreach tgtChildNode:
  //         foreach srcParentNode:
  //             foreach srcChildNode:
  //                 set corresponding block to shift matrix

  BfPoints2 srcChildPts, srcPts, tgtChildPts;

  blockIndex = 0;

  for (p0 = 0; p0 < bfPtrArraySize(tgtLevelNodes); ++p0) {
    bfPtrArrayGet(tgtLevelNodes, p0, (BfPtr *)&tgtNode);
    numTgtChildren = getChildren(tgtNode, tgtChild);
    // tgtCirc = bfGetQuadtreeNodeBoundingCircle(tgtNode);
    for (p1 = 0; p1 < numTgtChildren; ++p1) {
      tgtChildCirc = bfGetQuadtreeNodeBoundingCircle(tgtChild[p1]);

      for (q0 = 0; q0 < bfPtrArraySize(srcLevelNodes); ++q0) {
        bfPtrArrayGet(srcLevelNodes, q0, (BfPtr *)&srcNode);
        numSrcChildren = getChildren(srcNode, srcChild);
        srcCirc = bfGetQuadtreeNodeBoundingCircle(srcNode);
        for (q1 = 0; q1 < numSrcChildren; ++q1) {
          srcChildCirc = bfGetQuadtreeNodeBoundingCircle(srcChild[q1]);

          BfSize numRows = getRows(factor, factor->rowInd[blockIndex]);
          BfSize numCols = getCols(factor, factor->colInd[blockIndex]);
          assert(numRows > 0 && numCols > 0);

          /* sample points on the source child circle */
          srcChildPts = bfSamplePointsOnCircle2(&srcChildCirc, numCols);
          HANDLE_ERROR();

          /* sample points on the source circle */
          srcPts = bfSamplePointsOnCircle2(&srcCirc, numRows);
          HANDLE_ERROR();

          /* sample points on the target child circle */
          tgtChildPts = bfSamplePointsOnCircle2(&tgtChildCirc, numRows);
          HANDLE_ERROR();

          /* compute the shift matrix for this configuration of circles */
          getShiftMat(&srcChildPts, &srcPts, &tgtChildPts, K,
                      &factor->block[blockIndex]);
          HANDLE_ERROR();

          /* if we're debugging, store this block's points---free them
           * otherwise */
#if BF_DEBUG
          factor->srcPtsOrig[blockIndex] = srcChildPts;
          factor->srcPtsEquiv[blockIndex] = srcPts;
          factor->tgtPts[blockIndex] = tgtChildPts;
#else
          bfFreePoints2(&srcChildPts);
          bfFreePoints2(&srcPts);
          bfFreePoints2(&tgtChildPts);
#endif

          ++blockIndex;

        } /* end loop over all source children */
      }

    }
  } /* end loop over all target children */

  END_ERROR_HANDLING() {
    bfFreePoints2(&srcChildPts);
    bfFreePoints2(&srcPts);
    bfFreePoints2(&tgtChildPts);
    bfFreeFactor(factor);
  }
}

static void
makeLastFactor(BfFactor *factor, BfFactor const *prevFactor, BfReal K,
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
  assert(numBlocks == prevFactor->numBlockRows);

  BfQuadtreeNode const *srcNode = NULL;
  BfQuadtreeNode const *tgtNode = NULL;

  /* get the lone source node and its bounding circle */
  bfPtrArrayGetFirst(srcLevelNodes, (BfPtr *)&srcNode);
  BfCircle2 srcCirc = bfGetQuadtreeNodeBoundingCircle(srcNode);

  /* initialize the last butterfly factor */
  initDiagonalFactor(factor, numBlocks);
  HANDLE_ERROR();

  BfPoints2 srcCircPts, tgtPts;

  /* iterate over each node of the final level of the target node tree
   * and compute the matrix which will evaluate the potential at each
   * target point due to the charges on the source proxy points */
  for (BfSize i = 0; i < numBlocks; ++i) {
    BfSize prevNumRows = getRows(prevFactor, i);

    /* get the current target node */
    bfPtrArrayGet(tgtLevelNodes, i, (BfPtr *)&tgtNode);

    /* get the proxy points on the source circle for the current
     * source block */
    srcCircPts = bfSamplePointsOnCircle2(&srcCirc, prevNumRows);
    HANDLE_ERROR();

    /* get the current set of target points */
    bfGetQuadtreeNodePoints(tree, tgtNode, &tgtPts);
    HANDLE_ERROR();

    bfGetHelm2KernelMatrix(&factor->block[i], &srcCircPts, &tgtPts, K);
    HANDLE_ERROR();

    assert(factor->rowOffset[i + 1] == BF_SIZE_BAD_VALUE);
    factor->rowOffset[i + 1] = factor->block[i].numRows;

    assert(factor->colOffset[i + 1] == BF_SIZE_BAD_VALUE);
    factor->colOffset[i + 1] = factor->block[i].numCols;

    /* hang onto this block's points if we're in debug mode, and free
     * them otherwise */
#if BF_DEBUG
    factor->srcPtsOrig[i] = srcCircPts;
    factor->tgtPts[i] = tgtPts;
#else
    bfFreePoints2(&srcCircPts);
    bfFreePoints2(&tgtPts);
#endif
  }

  cumSumRowAndColOffsets(factor);

  END_ERROR_HANDLING() {
    bfFreePoints2(&srcCircPts);
    bfFreePoints2(&tgtPts);
    bfFreeFactor(factor);
  }
}

static bool
allSourceNodesHaveChildren(BfPtrArray const *srcLevelNodes)
{
  for (BfSize i = 0; i < bfPtrArraySize(srcLevelNodes); ++i) {
    BfQuadtreeNode *srcNode;
    bfPtrArrayGet(srcLevelNodes, i, (BfPtr *)&srcNode);
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
    BfQuadtreeNode *srcNode;
    bfPtrArrayGet(srcLevelNodes, i, (BfPtr *)&srcNode);

    BfSize numSrcPoints = bfQuadtreeNodeNumPoints(srcNode);
    BfCircle2 srcCirc = bfGetQuadtreeNodeBoundingCircle(srcNode);
    BfSize rank = bfHelm2RankEstForTwoCircles(&tgtCirc, &srcCirc, K, 1, 1e-15);

    if (rank > numSrcPoints)
      return false;
  }

  return true;
}

void
bfMakeFac(BfQuadtree const *tree,
          BfQuadtreeNode const *srcNode, BfQuadtreeNode const *tgtNode,
          BfReal K, BfSize *numFactors, BfFactor **factors)
{
  BEGIN_ERROR_HANDLING();

  /* set up the level iterator for the source tree---this iterator
   * goes from the leaves of the tree to the root (in reverse) */
  BfQuadtreeLevelIter srcLevelIter = bfInitQuadtreeLevelIter(
    BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER, (BfQuadtreeNode *)srcNode);
  HANDLE_ERROR();

  /* set up the level iterator for the target tree, which goes from
   * the root to the leaves */
  BfQuadtreeLevelIter tgtLevelIter = bfInitQuadtreeLevelIter(
    BF_TREE_TRAVERSAL_LR_LEVEL_ORDER, (BfQuadtreeNode *)tgtNode);
  HANDLE_ERROR();

  /* get the current source and target depths and make sure they're
   * compatible */
  BfSize currentSrcDepth = bfQuadtreeLevelIterCurrentDepth(&srcLevelIter);
  BfSize currentTgtDepth = bfQuadtreeLevelIterCurrentDepth(&tgtLevelIter);
  assert(currentTgtDepth <= currentSrcDepth);

  /* skip source levels until we're no deeper than the maximum depth
   * beneath the target node */
  BfSize maxDepthBelowTgtNode = bfGetMaxDepthBelowQuadtreeNode(tgtNode);
  while (currentSrcDepth > maxDepthBelowTgtNode) {
    bfQuadtreeLevelIterNext(&srcLevelIter);
    --currentSrcDepth;
  }

  /* skip source levels until we reach a level where each node has children
   *
   * TODO: I am *NOT AT ALL* sure this is right... but it should at
   * get me unstuck for now at least. */
  while (currentSrcDepth > currentTgtDepth &&
         !allSourceNodesHaveChildren(&srcLevelIter.levelNodes)) {
    bfQuadtreeLevelIterNext(&srcLevelIter);
    --currentSrcDepth;
  }
  assert(allSourceNodesHaveChildren(&srcLevelIter.levelNodes));

  /* TODO: step the source level iterator until:
   * - the rank estimate between the first target node and each source
   *   node is smaller than corresponding the number of points */

  while (currentSrcDepth > currentTgtDepth &&
         !allRankEstimatesAreOK(tgtNode, K, &srcLevelIter.levelNodes)) {
    bfQuadtreeLevelIterNext(&srcLevelIter);
    --currentSrcDepth;
  }
  assert(allRankEstimatesAreOK(tgtNode, K, &srcLevelIter.levelNodes));

  /* get number of factors in the butterfly factorization */
  *numFactors = currentSrcDepth - currentTgtDepth + 2;

  /* allocate space for the butterfly factors
   *
   * when multiplying, the factors will be multiplied in the order
   * that they're stored here (e.g., bf_factors[0] is the first factor
   * that will be multiplied) */
  *factors = malloc(*numFactors*sizeof(BfFactor));
  if (*factors == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfFactor *factor = *factors;

  /* make the first factor in the butterfly factorization: this is the
   * factor which initially shifts the charges on the source points to
   * the first level of source circles */
  makeFirstFactor(&factor[0], K, tree,
                  &srcLevelIter.levelNodes, &tgtLevelIter.levelNodes);
  HANDLE_ERROR();

  for (BfSize i = 1; i < *numFactors - 1; ++i) {
    /* go up a level on the source tree */
    bfQuadtreeLevelIterNext(&srcLevelIter);

    /* make the next factor */
    makeFactor(&factor[i], &factor[i - 1], K,
               &srcLevelIter.levelNodes, &tgtLevelIter.levelNodes);
    HANDLE_ERROR();

    /* go down a level on the target tree */
    bfQuadtreeLevelIterNext(&tgtLevelIter);
  }

  /* make the last factor in the butterfly factorization: this is the
   * evaluation factor, which computes the potential at each target
   * point due to the charges on the final source circle */
  makeLastFactor(&factor[*numFactors - 1], &factor[*numFactors - 2], K,
                 tree, &srcLevelIter.levelNodes, &tgtLevelIter.levelNodes);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    bfFreeFac(*numFactors, factors);
  }

  bfFreeQuadtreeLevelIter(&srcLevelIter);
  bfFreeQuadtreeLevelIter(&tgtLevelIter);
}

void bfFreeFac(BfSize numFactors, BfFactor **factorPtr) {
  BfFactor *factor = *factorPtr;
  for (BfSize i = 0; i < numFactors; ++i)
    bfFreeFactor(&factor[i]);
  free(factor);
  *factorPtr = NULL;
}

void bfMulFac(BfFactor const *factor, BfMat const *X, BfMat *Y) {
  BEGIN_ERROR_HANDLING();

  BfSize p = getNumCols(factor);
  if (X->numRows != p)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize m = getNumRows(factor), n = X->numCols;

  *Y = bfGetUninitializedMat();
  bfMatZeros(Y, BF_DTYPE_COMPLEX, m, n);
  HANDLE_ERROR();

  BfMat tmp;

  for (BfSize k = 0; k < factor->numBlocks; ++k) {
    BfSize i0 = factor->rowOffset[factor->rowInd[k]];
    BfSize i1 = factor->rowOffset[factor->rowInd[k] + 1];

    BfSize j0 = factor->colOffset[factor->colInd[k]];
    BfSize j1 = factor->colOffset[factor->colInd[k] + 1];

    BfMat Xrows = bfGetMatRowRange(X, j0, j1);

    tmp = bfGetUninitializedMat();
    bfMatMul(&factor->block[k], &Xrows, &tmp);
    HANDLE_ERROR();

    BfMat Yrows = bfGetMatRowRange(Y, i0, i1);
    bfMatAddInplace(&Yrows, &tmp);
    HANDLE_ERROR();

    bfFreeMat(&tmp);
  }

  END_ERROR_HANDLING() {
    bfFreeMat(&tmp);
    bfFreeMat(Y);
  }
}
