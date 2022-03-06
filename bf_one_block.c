#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>

#include "error_macros.h"
#include "fac.h"
#include "helm2.h"
#include "mat_block.h"
#include "mat_block_coo.h"
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
  BfMatProduct *factorization = NULL;
  BfMat *factor = NULL;

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

  factorization = bfFacMakeSingleLevelHelm2(&tree, srcNode, tgtNode, K);
  HANDLE_ERROR();
  printf("computed kernel matrix's butterfly factorization\n");

  numFactors = bfMatProductNumFactors(factorization);

  /* write factors to disk */

  char cmd[1024];
  sprintf(cmd, "for i in $(seq 0 %lu); do rm -rf ./factor$i;"
          "mkdir ./factor$i; done", numFactors - 1);
  system(cmd);

  char cwd[1024];
  assert(getcwd(cwd, 1024) != NULL);
  printf("cwd: %s\n", cwd);

  for (BfSize i = 0; i < numFactors; ++i) {
    factor = bfMatProductGetFactor(factorization, i);

    // TODO: do this cast with a function w/ error handling
    BfMatBlock *matBlock = (BfMatBlock *)factor;

    BfMatType matType = bfMatGetType(factor);
    BfSize numBlockRows = bfMatGetNumRows(factor);
    BfSize numBlockCols = bfMatGetNumCols(factor);
    BfSize numBlocks = bfMatBlockNumBlocks(matBlock);

    char path[1024];
    sprintf(path, "factor%lu", i);

    chdir(path);

    FILE *fp = fopen("info.txt", "w");
    fprintf(fp, "numBlockRows %lu\n", numBlockRows);
    fprintf(fp, "numBlockCols %lu\n", numBlockCols);
    fprintf(fp, "numBlocks %lu\n", numBlocks);
    fclose(fp);

    fp = fopen("rowInd.bin", "w");
    if (matType == BF_MAT_TYPE_BLOCK_COO)
      fwrite(((BfMatBlockCoo *)factor)->rowInd, sizeof(BfSize), numBlocks, fp);
    else if (matType == BF_MAT_TYPE_BLOCK_DIAG)
      for (BfSize j = 0; j < numBlocks; ++j)
        fwrite(&j, sizeof(BfSize), 1, fp);
    else
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);
    fclose(fp);

    fp = fopen("colInd.bin", "w");
    if (matType == BF_MAT_TYPE_BLOCK_COO)
      fwrite(((BfMatBlockCoo *)factor)->colInd, sizeof(BfSize), numBlocks, fp);
    else if (matType == BF_MAT_TYPE_BLOCK_DIAG)
      for (BfSize j = 0; j < numBlocks; ++j)
        fwrite(&j, sizeof(BfSize), 1, fp);
    else
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);
    fclose(fp);

    fp = fopen("rowOffset.bin", "w");
    fwrite(matBlock->rowOffset, sizeof(BfSize), numBlockRows + 1, fp);
    fclose(fp);

    fp = fopen("colOffset.bin", "w");
    fwrite(matBlock->colOffset, sizeof(BfSize), numBlockCols + 1, fp);
    fclose(fp);

    for (BfSize j = 0; j < numBlocks; ++j) {
      char filename[1024];

      sprintf(filename, "block%lu.bin", j);
      bfMatSave(matBlock->block[j], filename);

#ifdef BF_DEBUG
      BfPoints2 *auxPts = (BfPoints2 *)matBlock->block[j]->aux;
      BfPoints2 *srcChildPts = &auxPts[0];
      BfPoints2 *srcPts = &auxPts[1];
      BfPoints2 *tgtChildPts = &auxPts[2];

      sprintf(filename, "srcPtsOrig%lu.bin", j);
      bfSavePoints2(srcChildPts, filename);

      if (i != numFactors - 1) {
        sprintf(filename, "srcPtsEquiv%lu.bin", j);
        bfSavePoints2(srcPts, filename);
      }

      sprintf(filename, "tgtPts%lu.bin", j);
      bfSavePoints2(tgtChildPts, filename);
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
  // bfFreeFac(numFactors, &factor);
  bfMatDenseComplexDeinitAndDelete(&Z_gt);
  bfFreePoints2(&srcPts);
  bfFreePoints2(&tgtPts);
  bfFreeQuadtree(&tree);
  bfFreePoints2(&points);
}
