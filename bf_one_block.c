#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "fac.h"
#include "helm2.h"
#include "mat.h"
#include "rand.h"
#include "quadtree.h"

int main(int argc, char const *argv[]) {
  if (argc != 2) {
    printf("usage: %s <points.bin>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  enum BfError error = BF_ERROR_NO_ERROR;

  printf("reading points from %s...\n", argv[1]);

  BfPoints2 points;
  error = bfReadPoints2FromFile(argv[1], &points);
  assert(!error);

  printf("building quadtree...\n");

  BfQuadtree tree;
  bfInitQuadtreeFromPoints(&tree, &points);

  BfReal K = 3000;

  bfSeed(1234);

  /* Get source and target nodes from quadtree and check that their
   * indices are OK */

  size_t src_depth = 3, srcNodeIndex = 20;
  size_t tgt_depth = 3, tgtNodeIndex = 16*3 + 4*3 + 2;

  BfQuadtreeNode *srcNode;
  error = bfGetQuadtreeNode(&tree, src_depth, srcNodeIndex, &srcNode);
  assert(!error);

  BfQuadtreeNode *tgtNode;
  error = bfGetQuadtreeNode(&tree, tgt_depth, tgtNodeIndex, &tgtNode);
  assert(!error);

  /* Compute the groundtruth subblock of the kernel matrix induced by
   * the source and target nodes */

  BfPoints2 tgtPts;
  error = bfGetQuadtreeNodePoints(tgtNode, &tgtPts);
  assert(!error);

  BfPoints2 srcPts;
  error = bfGetQuadtreeNodePoints(srcNode, &srcPts);
  assert(!error);

  BfMat Z_gt = bfGetUninitializedMat();
  bfGetHelm2KernelMatrix(&Z_gt, &srcPts, &tgtPts, K);

  BfSize num_bytes;
  bfMatNumBytes(&Z_gt, &num_bytes);

  printf("computed groundtruth subblock of kernel matrix:\n");
  printf("- rows: %lu\n", Z_gt.numRows);
  printf("- columns: %lu\n", Z_gt.numCols);
  printf("- size: %1.2f MB\n", ((double)num_bytes)/(1024*1024));

  bfSaveMat(&Z_gt, "Z_gt.bin");
  printf("wrote Z_gt.bin\n");

  /* compute a butterfly factorization of the selected source and
   * target nodes */

  BfSize numFactors;
  BfFactor *factors;
  bfMakeFac(srcNode, tgtNode, K, &numFactors, &factors);

  /* cleanup */

  bfFreePoints2(&points);
}
