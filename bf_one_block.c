#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>

#include "fac.h"
#include "helm2.h"
#include "mat.h"
#include "rand.h"
#include "quadtree.h"

enum BfError func(BfQuadtree *tree, BfQuadtreeNode *node, void *arg) {
  (void)tree;
  (void)arg;

  BfPoints2 points;
  bfGetQuadtreeNodePoints(tree, node, &points);

  assert(bfBbox2ContainsPoints(&node->bbox, &points));

  bfFreePoints2(&points);

  return BF_ERROR_NO_ERROR;
}

int main(int argc, char const *argv[]) {
  if (argc != 2) {
    printf("usage: %s <points.bin>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  enum BfError error = BF_ERROR_NO_ERROR;

  BfPoints2 points;
  error = bfReadPoints2FromFile(argv[1], &points);
  assert(!error);
  printf("read points from %s\n", argv[1]);

  BfQuadtree tree;
  bfInitQuadtreeFromPoints(&tree, &points);
  puts("built quadtree");

  bfMapQuadtree(&tree, BF_TREE_TRAVERSAL_LR_LEVEL_ORDER, func, NULL);

  BfReal K = 3000;

  bfSeed(1234);

  /* Get source and target nodes from quadtree and check that their
   * indices are OK */

  BfSize src_depth = 3, srcNodeIndex = 20;
  BfSize tgt_depth = 3, tgtNodeIndex = 16*3 + 4*3 + 2;

  BfQuadtreeNode *srcNode;
  error = bfGetQuadtreeNode(&tree, src_depth, srcNodeIndex, &srcNode);
  assert(!error);

  BfQuadtreeNode *tgtNode;
  error = bfGetQuadtreeNode(&tree, tgt_depth, tgtNodeIndex, &tgtNode);
  assert(!error);

  /* Compute the groundtruth subblock of the kernel matrix induced by
   * the source and target nodes */

  BfPoints2 tgtPts;
  error = bfGetQuadtreeNodePoints(&tree, tgtNode, &tgtPts);
  assert(!error);

  bfSavePoints2(&tgtPts, "tgtPts.bin");
  puts("wrote target points to tgtPts.bin");

  BfPoints2 srcPts;
  error = bfGetQuadtreeNodePoints(&tree, srcNode, &srcPts);
  assert(!error);

  bfSavePoints2(&srcPts, "srcPts.bin");
  puts("wrote source points to srcPts.bin");

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
  BfFactor *factor;
  bfMakeFac(&tree, srcNode, tgtNode, K, &numFactors, &factor);
  printf("computed kernel matrix's butterfly factorization\n");

  /* write factors to disk */

  char cmd[1024];
  sprintf(cmd, "for i in $(seq 0 %lu); do rm -rf ./factor$i;"
          "mkdir ./factor$i; done", numFactors - 1);
  system(cmd);

  char cwd[1024];
  assert(getcwd(cwd, 1024) != NULL);
  printf("cwd: %s\n", cwd);

  for (BfSize i = 0; i < numFactors; ++i) {
    char path[1024];
    sprintf(path, "factor%lu", i);

    chdir(path);

    FILE *fp = fopen("info.txt", "w");
    fprintf(fp, "numBlockRows %lu\n", factor[i].numBlockRows);
    fprintf(fp, "numBlockCols %lu\n", factor[i].numBlockCols);
    fprintf(fp, "numBlocks %lu\n", factor[i].numBlocks);
    fclose(fp);

    fp = fopen("rowInd.bin", "w");
    fwrite(factor[i].rowInd, sizeof(BfSize), factor[i].numBlocks, fp);
    fclose(fp);

    fp = fopen("colInd.bin", "w");
    fwrite(factor[i].colInd, sizeof(BfSize), factor[i].numBlocks, fp);
    fclose(fp);

    fp = fopen("rowOffset.bin", "w");
    fwrite(factor[i].rowOffset, sizeof(BfSize), factor[i].numBlockRows + 1, fp);
    fclose(fp);

    fp = fopen("colOffset.bin", "w");
    fwrite(factor[i].colOffset, sizeof(BfSize), factor[i].numBlockCols + 1, fp);
    fclose(fp);

    for (BfSize j = 0; j < factor[i].numBlocks; ++j) {
      char filename[1024];

      sprintf(filename, "block%lu.bin", j);
      bfSaveMat(&factor[i].block[j], filename);

#ifdef BF_DEBUG
      sprintf(filename, "srcPtsOrig%lu.bin", j);
      bfSavePoints2(&factor[i].srcPtsOrig[j], filename);

      if (i != numFactors - 1) {
        sprintf(filename, "srcPtsEquiv%lu.bin", j);
        bfSavePoints2(&factor[i].srcPtsEquiv[j], filename);
      }

      sprintf(filename, "tgtPts%lu.bin", j);
      bfSavePoints2(&factor[i].tgtPts[j], filename);
#endif
    }

    chdir(cwd);
  }

  /* test multiplication */

  BfComplex q[3];
  bfRandn(6, (BfReal *)q);
  BfMat Q = bfGetUninitializedMat();
  bfInitEmptyMat(&Q, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE, srcPts.size, 1);
  for (BfSize i = 0; i < srcPts.size; ++i) {
    BfComplex *ptr;
    bfGetMatEltPtr(&Q, i, 0, (BfPtr *)&ptr);
    *ptr = srcPts.data[i][0]*q[0] + srcPts.data[i][1]*q[1] + q[2];
  }
  puts("set up test problem");

  BfMat Phi_gt = bfGetUninitializedMat();
  bfMatMul(&Z_gt, &Q, &Phi_gt);
  puts("did MVP with groundtruth kernel matrix");

  BfMat *Phi = malloc(numFactors*sizeof(BfMat));
  for (BfSize i = 0; i < numFactors; ++i)
    Phi[i] = bfGetUninitializedMat();
  bfMulFac(&factor[0], &Q, &Phi[0]);
  for (BfSize i = 1; i < numFactors; ++i)
    bfMulFac(&factor[i], &Phi[i - 1], &Phi[i]);
  puts("did MVPs with each butterfly factor");

  /* write simulation info to a text file */

  FILE *fp = fopen("info.txt", "w");
  fprintf(fp, "numSrcPts %lu\n", srcPts.size);
  fprintf(fp, "numTgtPts %lu\n", tgtPts.size);
  fprintf(fp, "numFactors %lu\n", numFactors);
  fprintf(fp, "K %g\n", K);
  fclose(fp);
  puts("wrote simulation information to info.txt");

  /* cleanup */

  for (BfSize i = 0; i < numFactors; ++i)
    bfFreeMat(&Phi[i]);
  free(Phi);
  bfFreeMat(&Phi_gt);
  bfFreeMat(&Q);
  bfFreeFac(numFactors, &factor);
  bfFreeMat(&Z_gt);
  bfFreePoints2(&srcPts);
  bfFreePoints2(&tgtPts);
  bfFreeQuadtree(&tree);
  bfFreePoints2(&points);
}
