#include "fac.h"

#include <assert.h>
#include <math.h>

#include "helm2.h"

static enum BfError
initEmptyFactor(BfFactor *factor, BfSize numBlockRows, BfSize numBlockCols,
                BfSize numBlocks)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  factor->numBlockRows = numBlockRows;
  factor->numBlockCols = numBlockCols;
  factor->numBlocks = numBlocks;

  /* allocate and initialize row indices to dummy values */
  factor->rowInd = malloc(numBlocks*sizeof(BfSize));
  if (factor->rowInd == NULL) {
    error = BF_ERROR_MEMORY_ERROR;
    goto cleanup;
  }
  for (BfSize i = 0; i < numBlocks; ++i)
    factor->rowInd[i] = BF_SIZE_BAD_VALUE;

  /* allocate and initialize column indices to dummy values */
  factor->colInd = malloc(numBlocks*sizeof(BfSize));
  if (factor->colInd == NULL) {
    error = BF_ERROR_MEMORY_ERROR;
    goto cleanup;
  }
  for (BfSize i = 0; i < numBlocks; ++i)
    factor->colInd[i] = BF_SIZE_BAD_VALUE;

  /* alloc and init the row index offsets for each block row */
  factor->rowOffset = malloc((numBlockRows + 1)*sizeof(BfSize));
  if (factor->rowOffset == NULL) {
    error = BF_ERROR_MEMORY_ERROR;
    goto cleanup;
  }
  for (BfSize i = 0; i < numBlockRows + 1; ++i)
    factor->rowOffset[i] = BF_SIZE_BAD_VALUE;

  /* alloc and init the column index offsets for each block column */
  factor->colOffset = malloc((numBlockCols + 1)*sizeof(BfSize));
  if (factor->colOffset == NULL) {
    error = BF_ERROR_MEMORY_ERROR;
    goto cleanup;
  }
  for (BfSize i = 0; i < numBlockCols + 1; ++i)
    factor->colOffset[i] = BF_SIZE_BAD_VALUE;

  /* allocate and initialize each nonzero (diagonal) block */
  factor->block = malloc(numBlocks*sizeof(BfMat));
  if (factor->block == NULL) {
    error = BF_ERROR_MEMORY_ERROR;
    goto cleanup;
  }
  for (BfSize i = 0; i < numBlocks; ++i)
    factor->block[i] = bfGetUninitializedMat();

  /* if we're debugging, allocate space in this factor to store each
   * block's points */
#if BF_DEBUG
  factor->srcPtsOrig = malloc(numBlocks*sizeof(BfPoints2));
  for (BfSize i = 0; i < numBlocks; ++i)
    factor->srcPtsOrig[i] = bfGetUninitializedPoints2();

  factor->srcPtsEquiv = malloc(numBlocks*sizeof(BfPoints2));
  for (BfSize i = 0; i < numBlocks; ++i)
    factor->srcPtsEquiv[i] = bfGetUninitializedPoints2();

  factor->tgtPts = malloc(numBlocks*sizeof(BfPoints2));
  for (BfSize i = 0; i < numBlocks; ++i)
    factor->tgtPts[i] = bfGetUninitializedPoints2();
#endif

cleanup:
  if (error)
    bfFreeFactor(factor);

  return error;
}

static enum BfError initDiagonalFactor(BfFactor *factor, BfSize numBlocks) {
  enum BfError error = BF_ERROR_NO_ERROR;

  error = initEmptyFactor(factor, numBlocks, numBlocks, numBlocks);
  if (error)
    goto cleanup;

  for (BfSize i = 0; i < numBlocks; ++i)
    factor->rowInd[i] = i;

  for (BfSize i = 0; i < numBlocks; ++i)
    factor->colInd[i] = i;

cleanup:
  if (error)
    bfFreeFactor(factor);

  return error;
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

static
enum BfError
getShiftMat(BfPoints2 const *srcPtsOrig, BfPoints2 const *srcPtsEquiv,
            BfPoints2 const *tgtPts, BfReal K, BfMat *Z_shift)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  /* compute the kernel matrix mapping charges on the original sources
   * points to potentials on the original target points */
  BfMat Z_orig = bfGetUninitializedMat();
  error = bfGetHelm2KernelMatrix(&Z_orig, srcPtsOrig, tgtPts, K);
  if (error) {
    bfFreeMat(&Z_orig);
    return error;
  }

  /* compute the kernel matrix mapping charges on the source
   * circle to potentials on the target circle */
  BfMat Z_equiv = bfGetUninitializedMat();
  error = bfGetHelm2KernelMatrix(&Z_equiv, srcPtsEquiv, tgtPts, K);
  if (error) {
    bfFreeMat(&Z_orig);
    bfFreeMat(&Z_equiv);
    return error;
  }

  /* set the "shift matrix" to Z_equiv\Z_orig */
  *Z_shift = bfGetUninitializedMat();
  error = bfMatLstSq(&Z_equiv, &Z_orig, Z_shift);
  if (error) {
    bfFreeMat(&Z_orig);
    bfFreeMat(&Z_equiv);
    bfFreeMat(Z_shift);
    return error;
  }

  bfFreeMat(&Z_orig);
  bfFreeMat(&Z_equiv);

  return error;
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

static enum BfError
makeFirstFactor(BfFactor *factor, BfReal K,
                BfQuadtree const *tree,
                BfPtrArray const *srcLevelNodes,
                BfPtrArray const *tgtLevelNodes)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  /* there should only be one target node at the starting level of
   * the target tree when we begin the traversal */
  if (bfPtrArraySize(tgtLevelNodes) != 1)
    return BF_ERROR_RUNTIME_ERROR;

  /* get the lone target node and its bounding circle */
  BfQuadtreeNode const *tgtNode;
  bfPtrArrayGetFirst(tgtLevelNodes, (BfPtr *)&tgtNode);
  BfCircle2 tgtCirc = bfGetQuadtreeNodeBoundingCircle(tgtNode);

  /* the current level of the source node tree shouldn't be empty */
  if (bfPtrArrayIsEmpty(srcLevelNodes))
    return BF_ERROR_RUNTIME_ERROR;

  /* the number of diagonal blocks in the first factor equals the
   * number of nodes at the current depth in the source tree */
  BfSize numBlocks = bfPtrArraySize(srcLevelNodes);

  /* initialize the first block diagonal butterfly factor */
  error = initDiagonalFactor(factor, numBlocks);
  assert(!error);

  /* iterate over each node at the starting level of the source tree
   * and compute the matrix which will maps the charges at the
   * original points in each source node to equivalent sources on a
   * circle surrounding the source box
   *
   * in the following loop, we initially set `rowOffset[i + 1]` and
   * `colOffset[i + 1]` to the number of rows in each block row and
   * column */
  for (BfSize i = 0; i < numBlocks; ++i) {
    /* get the current source node and its bounding circle */
    BfQuadtreeNode const *srcNode;
    bfPtrArrayGet(srcLevelNodes, i, (BfPtr *)&srcNode);
    BfCircle2 srcCirc = bfGetQuadtreeNodeBoundingCircle(srcNode);

    /* get the original points in the current source node box */
    BfPoints2 srcPts;
    bfGetQuadtreeNodePoints(tree, srcNode, &srcPts);

    /* verify that the source bounding circle contains the points */
    assert(bfCircle2ContainsPoints(&srcCirc, &srcPts));

    /* get the rank estimate for the current pair of source and
`     * target bounding circles */
    BfSize p = bfHelm2RankEstForTwoCircles(&srcCirc, &tgtCirc, K, 1, 1e-15);

    /* sample points on the source and target circles */
    BfPoints2 tgtCircPts = bfSamplePointsOnCircle2(&tgtCirc, p);
    BfPoints2 srcCircPts = bfSamplePointsOnCircle2(&srcCirc, p);

    /* compute the shift matrix and store it in the current block */
    error = getShiftMat(&srcPts, &srcCircPts, &tgtCircPts, K, &factor->block[i]);

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

    if (error)
      break;
  }

  cumSumRowAndColOffsets(factor);

  assert(!error); // TODO: if (error) { free blocks }

  return error;
}

/* This function computes an inner factor in a butterfly
 * factorization.
 *
 * ... explain how this works ...
 */
static
enum BfError
makeFactor(BfFactor *factor, BfFactor const *prevFactor, BfReal K,
           BfPtrArray const *srcLevelNodes, BfPtrArray const *tgtLevelNodes)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  /* first, we need to determine the number of blocks in each row and
   * column
   *
   * at the same time, we count the total number of blocks so we can
   * preallocate them */

  BfSize numBlockRows = 0, numBlockCols = 0, numBlocks = 0;

  for (BfSize p0 = 0, p = 0, i = 0, dj = 0; p0 < bfPtrArraySize(tgtLevelNodes); ++p0) {
    BfQuadtreeNode const *tgtNode;
    bfPtrArrayGet(tgtLevelNodes, p0, (BfPtr *)&tgtNode);

    BfQuadtreeNode const *tgtChild[4];
    BfSize numTgtChildren = getChildren(tgtNode, tgtChild);

    BfSize qmax = 0;

    for (BfSize p1 = 0; p1 < numTgtChildren; ++p1) {
      for (BfSize q0 = 0, q = 0; q0 < bfPtrArraySize(srcLevelNodes); ++q0) {
        BfQuadtreeNode const *srcNode;
        bfPtrArrayGet(srcLevelNodes, q0, (BfPtr *)&srcNode);

        BfQuadtreeNode const *srcChild[4]; // TODO: don't actually
                                           // need the source children
                                           // here, just their number
        BfSize numSrcChildren = getChildren(srcNode, srcChild);
        assert(numSrcChildren > 0);

        for (BfSize q1 = 0; q1 < numSrcChildren; ++q1) {
          BfSize j = q + dj;

          /* update the number of block columns and the total number
           * of blocks */
          numBlockCols = j + 1 > numBlockCols ? j + 1 : numBlockCols;

          ++numBlocks;

          ++q;
          qmax = q > qmax ? q : qmax;
        }

        /* update the number of block rows */
        numBlockRows = i + 1 > numBlockRows ? i + 1 : numBlockRows;

        ++i;
      }
      ++p;
    }
    dj += qmax;
  }

  initEmptyFactor(factor, numBlockRows, numBlockCols, numBlocks);

  for (BfSize i = 1; i <= numBlockRows; ++i)
    factor->rowOffset[i] = 0;

  for (BfSize j = 1; j <= numBlockCols; ++j)
    factor->colOffset[j] = 0;

  /* next, we set the number of rows and columns in each row or column of
   * blocks */

  for (BfSize p0 = 0, p = 0, i = 0, dj = 0; p0 < bfPtrArraySize(tgtLevelNodes); ++p0) {
    BfQuadtreeNode const *tgtNode;
    bfPtrArrayGet(tgtLevelNodes, p0, (BfPtr *)&tgtNode);

    BfQuadtreeNode const *tgtChild[4];
    BfSize numTgtChildren = getChildren(tgtNode, tgtChild);

    BfSize qmax = 0;

    BfCircle2 tgtCirc = bfGetQuadtreeNodeBoundingCircle(tgtNode);

    for (BfSize p1 = 0; p1 < numTgtChildren; ++p1) {
      BfCircle2 tgtChildCirc = bfGetQuadtreeNodeBoundingCircle(tgtChild[p1]);

      for (BfSize q0 = 0, q = 0; q0 < bfPtrArraySize(srcLevelNodes); ++q0) {
        BfQuadtreeNode const *srcNode;
        bfPtrArrayGet(srcLevelNodes, q0, (BfPtr *)&srcNode);

        BfCircle2 srcCirc = bfGetQuadtreeNodeBoundingCircle(srcNode);

        BfQuadtreeNode const *srcChild[4];
        BfSize numSrcChildren = getChildren(srcNode, srcChild);

        for (BfSize q1 = 0; q1 < numSrcChildren; ++q1) {
          BfSize j = q + dj;

          /* the number of original source points matches the number
           * of rows in the jth block of the previous factor */
          BfSize numSrcChildPts = getRows(prevFactor, j);

          /* set number of columns for current block: we find the
           * maximum over each number of source child points */
          assert(j < numBlockCols);
          if (numSrcChildPts > factor->colOffset[j + 1])
            factor->colOffset[j + 1] = numSrcChildPts;

          /* get the bounding circle for the current source child */
          BfCircle2 srcChildCirc = bfGetQuadtreeNodeBoundingCircle(srcChild[q1]);

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
          assert(i < factor->numBlockRows);
          if (rank > factor->rowOffset[i + 1])
            factor->rowOffset[i + 1] = rank;

          ++q;
          qmax = q > qmax ? q : qmax;
        }
        ++i;
      }
      ++p;
    }
    dj += qmax;
  }

  cumSumRowAndColOffsets(factor);

  // foreach tgtNode:
  //     foreach tgtChildNode:
  //         foreach srcParentNode:
  //             foreach srcChildNode:
  //                 set corresponding block to shift matrix

  BfSize blockIndex = 0;

  for (BfSize p0 = 0, p = 0, i = 0, dj = 0; p0 < bfPtrArraySize(tgtLevelNodes); ++p0) {
    BfQuadtreeNode const *tgtNode;
    bfPtrArrayGet(tgtLevelNodes, p0, (BfPtr *)&tgtNode);

    BfQuadtreeNode const *tgtChild[4];
    BfSize numTgtChildren = getChildren(tgtNode, tgtChild);

    BfSize qmax = 0;

    for (BfSize p1 = 0; p1 < numTgtChildren; ++p1) {
      BfCircle2 tgtChildCirc = bfGetQuadtreeNodeBoundingCircle(tgtChild[p1]);

      for (BfSize q0 = 0, q = 0; q0 < bfPtrArraySize(srcLevelNodes); ++q0) {
        BfQuadtreeNode const *srcNode;
        bfPtrArrayGet(srcLevelNodes, q0, (BfPtr *)&srcNode);

        BfCircle2 srcCirc = bfGetQuadtreeNodeBoundingCircle(srcNode);

        BfQuadtreeNode const *srcChild[4];
        BfSize numSrcChildren = getChildren(srcNode, srcChild);

        for (BfSize q1 = 0; q1 < numSrcChildren; ++q1) {
          BfSize j = q + dj;

          BfCircle2 srcChildCirc = bfGetQuadtreeNodeBoundingCircle(srcChild[q1]);

          BfSize numRows = getRows(factor, i);
          BfSize numCols = getCols(factor, j);
          assert(numRows > 0 && numCols > 0);

          /* sample points on each of the circles */
          BfPoints2 srcChildPts = bfSamplePointsOnCircle2(&srcChildCirc, numCols);
          BfPoints2 srcPts = bfSamplePointsOnCircle2(&srcCirc, numRows);
          BfPoints2 tgtChildPts = bfSamplePointsOnCircle2(&tgtChildCirc, numRows);

          /* compute the shift matrix for this configuration of circles */
          BfMat *Z_shift = &factor->block[blockIndex];
          error = getShiftMat(&srcChildPts, &srcPts, &tgtChildPts, K, Z_shift);
          assert(!error);

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

          /* set block row and column indices */
          assert(blockIndex < factor->numBlocks);
          factor->rowInd[blockIndex] = i;
          factor->colInd[blockIndex] = j;

          ++blockIndex;
          ++q;
          qmax = q > qmax ? q : qmax;
        }
        ++i;
      }
      ++p;
    }
    dj += qmax;
  }

  assert(!error); // TODO: if (error) { free stuff }

  return error;
}

static enum BfError
makeLastFactor(BfFactor *factor, BfFactor const *prevFactor, BfReal K,
               BfQuadtree const *tree,
               BfPtrArray const *srcLevelNodes,
               BfPtrArray const *tgtLevelNodes)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  /* the should only be one source node when we make the last
   * butterfly factor */
  if (bfPtrArraySize(srcLevelNodes) != 1)
    return BF_ERROR_RUNTIME_ERROR;

  /* get the lone source node and its bounding circle */
  BfQuadtreeNode const *srcNode;
  bfPtrArrayGetFirst(srcLevelNodes, (BfPtr *)&srcNode);
  BfCircle2 srcCirc = bfGetQuadtreeNodeBoundingCircle(srcNode);

  /* the current level of the target node tree shouldn't be empty */
  if (bfPtrArrayIsEmpty(tgtLevelNodes))
    return BF_ERROR_RUNTIME_ERROR;

  /* the number of diagonal blocks in the last factor equals the
   * number of nodes at the current depth of the target tree ... */
  BfSize numBlocks = bfPtrArraySize(tgtLevelNodes);

  /* ... and this number of blocks should match the number of block
   *  rows of the previous factor */
  if (numBlocks != prevFactor->numBlockRows)
    return BF_ERROR_RUNTIME_ERROR;

  /* initialize the last butterfly factor */
  error = initDiagonalFactor(factor, numBlocks);
  assert(!error);

  /* iterate over each node of the final level of the target node tree
   * and compute the matrix which will evaluate the potential at each
   * target point due to the charges on the source proxy points */
  for (BfSize i = 0; i < numBlocks; ++i) {
    BfSize prevNumRows = getRows(prevFactor, i);

    /* get the proxy points on the source circle for the current
     * source block */
    BfPoints2 srcCircPts = bfSamplePointsOnCircle2(&srcCirc, prevNumRows);

    /* get the current target node */
    BfQuadtreeNode const *tgtNode;
    bfPtrArrayGet(tgtLevelNodes, i, (BfPtr *)&tgtNode);

    /* get the current set of target points */
    BfPoints2 tgtPts;
    bfGetQuadtreeNodePoints(tree, tgtNode, &tgtPts);

    error = bfGetHelm2KernelMatrix(&factor->block[i], &srcCircPts, &tgtPts, K);
    assert(!error);

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

    if (error)
      break;
  }

  cumSumRowAndColOffsets(factor);

  assert(!error); // TODO: if (error) { free stuff }

  return error;
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

enum BfError
bfMakeFac(BfQuadtree const *tree,
          BfQuadtreeNode const *srcNode, BfQuadtreeNode const *tgtNode,
          BfReal K, BfSize *numFactors, BfFactor **factors)
{
  /* Let's try manually doing a couple levels */

  BfQuadtreeLevelIter srcLevelIter;
  bfInitQuadtreeLevelIter(
    &srcLevelIter, BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER,
    (BfQuadtreeNode *)srcNode);

  BfQuadtreeLevelIter tgtLevelIter;
  bfInitQuadtreeLevelIter(
    &tgtLevelIter, BF_TREE_TRAVERSAL_LR_LEVEL_ORDER,
    (BfQuadtreeNode *)tgtNode);

  BfSize currentSrcDepth;
  bfQuadtreeLevelIterCurrentDepth(&srcLevelIter, &currentSrcDepth);

  BfSize currentTgtDepth;
  bfQuadtreeLevelIterCurrentDepth(&tgtLevelIter, &currentTgtDepth);

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
  while (currentSrcDepth > 0 &&
         !allSourceNodesHaveChildren(&srcLevelIter.levelNodes)) {
    bfQuadtreeLevelIterNext(&srcLevelIter);
    --currentSrcDepth;
  }
  assert(allSourceNodesHaveChildren(&srcLevelIter.levelNodes));

  /* TODO: step the source level iterator until:
   * - the rank estimate between the first target node and each source
   *   node is smaller than corresponding the number of points */

  while (currentSrcDepth > 0 &&
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

  BfFactor *factor = *factors;

  makeFirstFactor(&factor[0], K, tree,
                  &srcLevelIter.levelNodes,
                  &tgtLevelIter.levelNodes);

  for (BfSize i = 1; i < *numFactors - 1; ++i) {
    /* go up a level on the source tree */
    bfQuadtreeLevelIterNext(&srcLevelIter);

    /* make the next factor */
    makeFactor(&factor[i], &factor[i - 1], K,
               &srcLevelIter.levelNodes, &tgtLevelIter.levelNodes);

    /* go down a level on the target tree */
    bfQuadtreeLevelIterNext(&tgtLevelIter);
  }

  makeLastFactor(&factor[*numFactors - 1], &factor[*numFactors - 2], K,
                 tree, &srcLevelIter.levelNodes, &tgtLevelIter.levelNodes);

  bfFreeQuadtreeLevelIter(&srcLevelIter);
  bfFreeQuadtreeLevelIter(&tgtLevelIter);

  return BF_ERROR_NO_ERROR;
}

void bfFreeFac(BfSize numFactors, BfFactor **factorPtr) {
  BfFactor *factor = *factorPtr;
  for (BfSize i = 0; i < numFactors; ++i)
    bfFreeFactor(&factor[i]);
  free(factor);
  *factorPtr = NULL;
}

enum BfError
bfMulFac(BfFactor const *factor, BfMat const *X, BfMat *Y)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  BfSize p = getNumCols(factor);

  if (X->numRows != p) {
    error = BF_ERROR_BAD_SHAPE;
    goto cleanup;
  }

  BfSize m = getNumRows(factor), n = X->numCols;

  *Y = bfGetUninitializedMat();
  error = bfMatZeros(Y, BF_DTYPE_COMPLEX, m, n);
  assert(!error);

  for (BfSize k = 0; k < factor->numBlocks; ++k) {
    BfSize i0 = factor->rowOffset[factor->rowInd[k]];
    BfSize i1 = factor->rowOffset[factor->rowInd[k] + 1];

    BfSize j0 = factor->colOffset[factor->colInd[k]];
    BfSize j1 = factor->colOffset[factor->colInd[k] + 1];

    BfMat Xrows = bfGetMatRowRange(X, j0, j1);

    BfMat tmp = bfGetUninitializedMat();
    error = bfMatMul(&factor->block[k], &Xrows, &tmp);
    assert(!error);

    BfMat Yrows = bfGetMatRowRange(Y, i0, i1);
    error = bfMatAddInplace(&Yrows, &tmp);
    assert(!error);

    bfFreeMat(&tmp);
  }

cleanup:
  if (error)
    bfFreeMat(Y);

  return error;
}
