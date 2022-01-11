#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>

#include "error_macros.h"
#include "fac.h"
#include "helm2.h"
#include "mat_dense_complex.h"
#include "rand.h"
#include "quadtree.h"

int main(int argc, char const *argv[]) {
  if (argc != 2) {
    printf("usage: %s <points.bin>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  BEGIN_ERROR_HANDLING();

  BfSize numFactors = BF_SIZE_BAD_VALUE;
  BfFactor *factor = NULL;

  BfPoints2 points;
  bfReadPoints2FromFile(argv[1], &points);
  HANDLE_ERROR();
  printf("read points from %s\n", argv[1]);

  BfQuadtree tree;
  bfInitQuadtreeFromPoints(&tree, &points);
  HANDLE_ERROR();
  puts("built quadtree");

  BfReal K = 3000;

  bfSeed(1234);

  /* Get source and target nodes from quadtree and check that their
   * indices are OK */

  BfSize src_depth = 3, srcNodeIndex = 20;
  BfSize tgt_depth = 3, tgtNodeIndex = 16*3 + 4*3 + 2;

  BfQuadtreeNode *srcNode;
  bfGetQuadtreeNode(&tree, src_depth, srcNodeIndex, &srcNode);
  HANDLE_ERROR();

  BfQuadtreeNode *tgtNode;
  bfGetQuadtreeNode(&tree, tgt_depth, tgtNodeIndex, &tgtNode);
  HANDLE_ERROR();

  /* Compute the groundtruth subblock of the kernel matrix induced by
   * the source and target nodes */

  BfPoints2 tgtPts;
  bfGetQuadtreeNodePoints(&tree, tgtNode, &tgtPts);
  HANDLE_ERROR();

  bfSavePoints2(&tgtPts, "tgtPts.bin");
  HANDLE_ERROR();
  puts("wrote target points to tgtPts.bin");

  BfPoints2 srcPts;
  bfGetQuadtreeNodePoints(&tree, srcNode, &srcPts);
  HANDLE_ERROR();

  bfSavePoints2(&srcPts, "srcPts.bin");
  HANDLE_ERROR();
  puts("wrote source points to srcPts.bin");

  BfMatDenseComplex *Z_gt = bfGetHelm2KernelMatrix(&srcPts, &tgtPts, K);
  HANDLE_ERROR();

  printf("computed groundtruth subblock of kernel matrix\n");

  bfMatDenseComplexSave(Z_gt, "Z_gt.bin");
  HANDLE_ERROR();
  printf("wrote Z_gt.bin\n");

  /* compute a butterfly factorization of the selected source and
   * target nodes */

  bfMakeFac(&tree, srcNode, tgtNode, K, &numFactors, &factor);
  HANDLE_ERROR();
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
    fprintf(fp, "numBlockRows %lu\n", factor[i].mat->super.numRows);
    fprintf(fp, "numBlockCols %lu\n", factor[i].mat->super.numCols);
    fprintf(fp, "numBlocks %lu\n", factor[i].mat->numBlocks);
    fclose(fp);

    fp = fopen("rowInd.bin", "w");
    fwrite(factor[i].mat->rowInd, sizeof(BfSize), factor[i].mat->numBlocks, fp);
    fclose(fp);

    fp = fopen("colInd.bin", "w");
    fwrite(factor[i].mat->colInd, sizeof(BfSize), factor[i].mat->numBlocks, fp);
    fclose(fp);

    fp = fopen("rowOffset.bin", "w");
    fwrite(factor[i].mat->rowOffset, sizeof(BfSize), factor[i].mat->super.numRows + 1, fp);
    fclose(fp);

    fp = fopen("colOffset.bin", "w");
    fwrite(factor[i].mat->colOffset, sizeof(BfSize), factor[i].mat->super.numCols + 1, fp);
    fclose(fp);

    for (BfSize j = 0; j < factor[i].mat->numBlocks; ++j) {
      char filename[1024];

      sprintf(filename, "block%lu.bin", j);
      bfMatSave(factor[i].mat->block[j], filename);

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

  BfMatDenseComplex *Q = bfMatDenseComplexNew();
  bfMatDenseComplexInit(Q, srcPts.size, 1);

  for (BfSize i = 0; i < srcPts.size; ++i)
    Q->data[i] = srcPts.data[i][0]*q[0] + srcPts.data[i][1]*q[1] + q[2];

  puts("set up test problem");

  /* write simulation info to a text file */

  FILE *fp = fopen("info.txt", "w");
  fprintf(fp, "numSrcPts %lu\n", srcPts.size);
  fprintf(fp, "numTgtPts %lu\n", tgtPts.size);
  fprintf(fp, "numFactors %lu\n", numFactors);
  fprintf(fp, "K %g\n", K);
  fclose(fp);
  puts("wrote simulation information to info.txt");

  END_ERROR_HANDLING() {
    puts("error!");
  }

  bfMatDenseComplexDeinitAndDelete(&Q);
  bfFreeFac(numFactors, &factor);
  bfMatDenseComplexDeinitAndDelete(&Z_gt);
  bfFreePoints2(&srcPts);
  bfFreePoints2(&tgtPts);
  bfFreeQuadtree(&tree);
  bfFreePoints2(&points);
}
