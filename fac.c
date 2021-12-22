#include "fac.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

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

  /* allocate and initialize the number of rows corresponding to each
   * row of blocks */
  factor->numRows = malloc(numBlockRows*sizeof(BfSize));
  if (factor->numRows == NULL) {
    error = BF_ERROR_MEMORY_ERROR;
    goto cleanup;
  }
  for (BfSize i = 0; i < numBlockRows; ++i)
    factor->numRows[i] = 0;

  /* allocate and initialize the number of columns corresponding to
   * column of blocks */
  factor->numCols = malloc(numBlockCols*sizeof(BfSize));
  if (factor->numCols == NULL) {
    error = BF_ERROR_MEMORY_ERROR;
    goto cleanup;
  }
  for (BfSize i = 0; i < numBlockCols; ++i)
    factor->numCols[i] = 0;

  /* allocat and initialize each nonzero (diagonal) block */
  factor->block = malloc(numBlocks*sizeof(BfMat));
  if (factor->block == NULL) {
    error = BF_ERROR_MEMORY_ERROR;
    goto cleanup;
  }
  for (BfSize i = 0; i < numBlocks; ++i)
    factor->block[i] = bfGetUninitializedMat();

cleanup:
  if (error) {
    free(factor->rowInd);
    free(factor->colInd);
    free(factor->numRows);
    free(factor->numCols);
    free(factor->block);
  }

  return error;
}

static enum BfError initDiagonalFactor(BfFactor *factor, BfSize numBlocks) {
  enum BfError error = BF_ERROR_NO_ERROR;

  factor->numBlockRows = numBlocks;
  factor->numBlockCols = numBlocks;
  factor->numBlocks = numBlocks;

  /* allocate and initialize row indices for each block */
  factor->rowInd = malloc(numBlocks*sizeof(BfSize));
  if (factor->rowInd == NULL) {
    error = BF_ERROR_MEMORY_ERROR;
    goto cleanup;
  }
  for (BfSize i = 0; i < numBlocks; ++i)
    factor->rowInd[i] = i;

  /* allocate and initialize column indices for each block */
  factor->colInd = malloc(numBlocks*sizeof(BfSize));
  if (factor->colInd == NULL) {
    error = BF_ERROR_MEMORY_ERROR;
    goto cleanup;
  }
  for (BfSize i = 0; i < numBlocks; ++i)
    factor->colInd[i] = i;

  /* allocate and initialize the number of rows corresponding to each
   * row of blocks */
  factor->numRows = malloc(numBlocks*sizeof(BfSize));
  if (factor->numRows == NULL) {
    error = BF_ERROR_MEMORY_ERROR;
    goto cleanup;
  }
  for (BfSize i = 0; i < numBlocks; ++i)
    factor->numRows[i] = BF_SIZE_BAD_VALUE;

  /* allocate and initialize the number of columns corresponding to
   * column of blocks */
  factor->numCols = malloc(numBlocks*sizeof(BfSize));
  if (factor->numCols == NULL) {
    error = BF_ERROR_MEMORY_ERROR;
    goto cleanup;
  }
  for (BfSize i = 0; i < numBlocks; ++i)
    factor->numCols[i] = BF_SIZE_BAD_VALUE;

  /* allocat and initialize each nonzero (diagonal) block */
  factor->block = malloc(numBlocks*sizeof(BfMat));
  if (factor->block == NULL) {
    error = BF_ERROR_MEMORY_ERROR;
    goto cleanup;
  }
  for (BfSize i = 0; i < numBlocks; ++i)
    factor->block[i] = bfGetUninitializedMat();

cleanup:
  if (error) {
    free(factor->rowInd);
    free(factor->colInd);
    free(factor->numRows);
    free(factor->numCols);
    free(factor->block);
  }

  return error;
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

  bfSaveMat(&Z_orig, "Z_or.bin");

  /* compute the kernel matrix mapping charges on the source
   * circle to potentials on the target circle */
  BfMat Z_equiv = bfGetUninitializedMat();
  error = bfGetHelm2KernelMatrix(&Z_equiv, srcPtsEquiv, tgtPts, K);
  if (error) {
    bfFreeMat(&Z_orig);
    bfFreeMat(&Z_equiv);
    return error;
  }

  bfSaveMat(&Z_equiv, "Z_eq.bin");

  /* set the "shift matrix" to Z_equiv\Z_orig */
  *Z_shift = bfGetUninitializedMat();
  error = bfMatLstSq(&Z_equiv, &Z_orig, Z_shift);
  if (error) {
    bfFreeMat(&Z_orig);
    bfFreeMat(&Z_equiv);
    bfFreeMat(Z_shift);
    return error;
  }

  bfSaveMat(Z_shift, "Z_shift.bin");

  return error;
}

static enum BfError
makeFirstFactor(BfFactor *factor, BfReal K,
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
   * circle surrounding the source box */
  for (BfSize i = 0; i < numBlocks; ++i) {
    /* get the current source node and its bounding circle */
    BfQuadtreeNode const *srcNode;
    bfPtrArrayGet(srcLevelNodes, i, (BfPtr *)&srcNode);
    BfCircle2 srcCirc = bfGetQuadtreeNodeBoundingCircle(srcNode);

    /* get the original points in the current source node box */
    BfPoints2 srcPts;
    bfGetQuadtreeNodePoints(srcNode, &srcPts);

    /* get the rank estimate for the current pair of source and
`     * target bounding circles */
    BfSize p;
    bfHelm2RankEstForTwoCircles(srcCirc, tgtCirc, K, 1, 1e-15, &p);

    /* sample points on the source and target circles */
    BfPoints2 tgtCircPts = bfSamplePointsOnCircle2(&tgtCirc, p);
    BfPoints2 srcCircPts = bfSamplePointsOnCircle2(&srcCirc, p);

    /* compute the shift matrix and store it in the current block */
    error = getShiftMat(&srcPts, &srcCircPts, &tgtCircPts, K, &factor->block[i]);

    /* update factor's number of rows for this block */
    if (factor->numRows[i] == BF_SIZE_BAD_VALUE)
      factor->numRows[i] = factor->block[i].numRows;
    else
      assert(factor->numRows[i] == factor->block[i].numRows);

    /* update factor's number of columns for this block */
    if (factor->numCols[i] == BF_SIZE_BAD_VALUE)
      factor->numCols[i] = factor->block[i].numCols;
    else
      assert(factor->numCols[i] == factor->block[i].numCols);

    /* free all temporaries */
    bfFreePoints2(&srcPts);
    bfFreePoints2(&tgtCircPts);
    bfFreePoints2(&srcCircPts);

    if (error)
      break;
  }

  // TODO: if (error) { free blocks }

  return error;
}

/* This function computes an inner factor in a butterfly
 * factorization.
 *
 * ... explain how this works ...
 */
static
enum BfError
makeFactor(BfSize factorIndex,
           BfFactor *factor, BfFactor const *prevFactor, BfReal K,
           BfPtrArray const *srcLevelNodes, BfPtrArray const *tgtLevelNodes)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  char filename[128];
  sprintf(filename, "ij%lu.txt", factorIndex);

  FILE *fp = fopen(filename, "w");

  printf("makeFactor()\n");

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

  printf("* numBlockRows: %lu\n", numBlockRows);
  printf("* numBlockCols: %lu\n", numBlockCols);
  printf("* numBlocks: %lu\n", numBlocks);

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
          BfSize numSrcChildPts = prevFactor->numRows[j];

          /* set number of columns for current block */
          assert(j < numBlockCols);
          if (numSrcChildPts > factor->numCols[j])
            factor->numCols[j] = numSrcChildPts;

          /* get the bounding circle for the current source child */
          BfCircle2 srcChildCirc = bfGetQuadtreeNodeBoundingCircle(srcChild[q1]);

          /* a priori rank estimate for the original circles */
          BfSize rankOr;
          bfHelm2RankEstForTwoCircles(srcChildCirc,tgtCirc,K,1,1e-15,&rankOr);

          /* a priori rank estimate for the new circles */
          BfSize rankEq;
          bfHelm2RankEstForTwoCircles(srcCirc,tgtChildCirc,K,1,1e-15,&rankEq);

          /* use the larger of the two rank estimates... not sure if
           * this is totally necessary, probably being a little
           * paranoid... they should be nearly the same, but might
           * differ a little due to rounding */
          BfSize rank = rankOr > rankEq ? rankOr : rankEq;

          /* update number of rows for current block */
          assert(i < factor->numBlockRows);
          if (rank > factor->numRows[i])
            factor->numRows[i] = rank;

          ++q;
          qmax = q > qmax ? q : qmax;
        }

        ++i;
      }
      ++p;
    }
    dj += qmax;
  }

  // foreach tgtNode:
  //     foreach tgtChildNode:
  //         foreach srcParentNode:
  //             foreach srcChildNode:
  //                 set corresponding block to shift matrix

  BfSize blockIndex = 0;

  for (BfSize p0 = 0, p = 0, i = 0, dj = 0; p0 < bfPtrArraySize(tgtLevelNodes); ++p0) {
    printf("next target node\n");

    BfQuadtreeNode const *tgtNode;
    bfPtrArrayGet(tgtLevelNodes, p0, (BfPtr *)&tgtNode);

    BfQuadtreeNode const *tgtChild[4];
    BfSize numTgtChildren = getChildren(tgtNode, tgtChild);

    BfSize qmax = 0;

    for (BfSize p1 = 0; p1 < numTgtChildren; ++p1) {
      printf("  next target child\n");

      BfCircle2 tgtChildCirc = bfGetQuadtreeNodeBoundingCircle(tgtChild[p1]);

      for (BfSize q0 = 0, q = 0; q0 < bfPtrArraySize(srcLevelNodes); ++q0) {
        printf("    next source node\n");

        BfQuadtreeNode const *srcNode;
        bfPtrArrayGet(srcLevelNodes, q0, (BfPtr *)&srcNode);

        BfCircle2 srcCirc = bfGetQuadtreeNodeBoundingCircle(srcNode);

        BfQuadtreeNode const *srcChild[4];
        BfSize numSrcChildren = getChildren(srcNode, srcChild);

        for (BfSize q1 = 0; q1 < numSrcChildren; ++q1) {
          BfSize j = q + dj;

          printf("    * (i, j) = (%lu, %lu)\n", i, j);
          fprintf(fp, "%lu %lu\n", i, j);

          BfCircle2 srcChildCirc = bfGetQuadtreeNodeBoundingCircle(srcChild[q1]);

          /* sample points on each of the circles */
          BfPoints2 srcChildPts = bfSamplePointsOnCircle2(&srcChildCirc, factor->numCols[j]);
          BfPoints2 srcPts = bfSamplePointsOnCircle2(&srcCirc, factor->numRows[i]);
          BfPoints2 tgtChildPts = bfSamplePointsOnCircle2(&tgtChildCirc, factor->numRows[i]);

          /* compute the shift matrix for this configuration of circles */
          BfMat *Z_shift = &factor->block[blockIndex];
          error = getShiftMat(&srcChildPts, &srcPts, &tgtChildPts, K, Z_shift);
          assert(!error);

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

      printf("  qmax = %lu\n", qmax);

      ++p;
    }

    dj += qmax;
    printf("dj = %lu\n", dj);
  }

  fclose(fp);

  return error;
}

static enum BfError
makeLastFactor(BfFactor *factor, BfFactor const *prevFactor, BfReal K,
               BfPtrArray const *srcLevelNodes, BfPtrArray const *tgtLevelNodes)
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
    /* get the proxy points on the source circle for the current
     * source block */
    BfPoints2 srcCircPts = bfSamplePointsOnCircle2(&srcCirc, prevFactor->numRows[i]);

    /* get the current target node */
    BfQuadtreeNode const *tgtNode;
    bfPtrArrayGet(tgtLevelNodes, i, (BfPtr *)&tgtNode);

    /* get the current set of target points */
    BfPoints2 tgtPts;
    bfGetQuadtreeNodePoints(tgtNode, &tgtPts);

    error = bfGetHelm2KernelMatrix(&factor->block[i], &srcCircPts, &tgtPts, K);
    assert(!error);

    assert(factor->numRows[i] == BF_SIZE_BAD_VALUE);
    factor->numRows[i] = factor->block[i].numRows;

    assert(factor->numCols[i] == BF_SIZE_BAD_VALUE);
    factor->numCols[i] = factor->block[i].numCols;

    bfFreePoints2(&srcCircPts);
    bfFreePoints2(&tgtPts);

    if (error)
      break;
  }

  // TODO: if (error) { free stuff }

  return error;
}

enum BfError
bfMakeFac(BfQuadtreeNode const *srcNode, BfQuadtreeNode const *tgtNode,
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

  // skip a few levels on the source side
  for (BfSize i = 0; i < 5; ++i)
    bfQuadtreeLevelIterNext(&srcLevelIter);

  BfSize current_src_depth;
  bfQuadtreeLevelIterCurrentDepth(&srcLevelIter, &current_src_depth);
  printf("starting source depth: %lu\n", current_src_depth);

  BfSize current_tgt_depth;
  bfQuadtreeLevelIterCurrentDepth(&tgtLevelIter, &current_tgt_depth);
  printf("starting target depth: %lu\n", current_tgt_depth);

  assert(current_tgt_depth <= current_src_depth);

  /* get number of factors in the butterfly factorization */
  *numFactors = current_src_depth - current_tgt_depth + 2;
  printf("computing butterfly factorization with %lu factors\n", *numFactors);

  /* allocate space for the butterfly factors
   *
   * when multiplying, the factors will be multiplied in the order
   * that they're stored here (e.g., bf_factors[0] is the first factor
   * that will be multiplied) */
  *factors = malloc(*numFactors*sizeof(BfFactor));

  BfFactor *factor = *factors;

  makeFirstFactor(&factor[0], K, &srcLevelIter.level_nodes,
                  &tgtLevelIter.level_nodes);

  for (BfSize i = 1; i < *numFactors - 1; ++i) {
    /* go up a level on the source tree */
    bfQuadtreeLevelIterNext(&srcLevelIter);

    /* make the next factor */
    makeFactor(i, &factor[i], &factor[i - 1], K,
               &srcLevelIter.level_nodes, &tgtLevelIter.level_nodes);

    /* go down a level on the target tree */
    bfQuadtreeLevelIterNext(&tgtLevelIter);
  }

  makeLastFactor(&factor[*numFactors - 1], &factor[*numFactors - 2], K,
                 &srcLevelIter.level_nodes, &tgtLevelIter.level_nodes);

  return BF_ERROR_NO_ERROR;
}
