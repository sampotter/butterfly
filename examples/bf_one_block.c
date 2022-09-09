#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>

#include <openblas/cblas.h>

#include <bf/error_macros.h>
#include <bf/fac.h>
#include <bf/helm2.h>
#include <bf/mat_block.h>
#include <bf/mat_block_coo.h>
#include <bf/mat_dense_complex.h>
#include <bf/points.h>
#include <bf/quadtree.h>
#include <bf/rand.h>

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

  bfQuadtreeSaveBoxesToTextFile(&tree, "quadtree.txt");
  puts("saved quadtree boxes to quadtree.txt");

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

  bfMatDenseComplexSave(bfMatDenseComplexToMat(Z_gt), "Z_gt.bin");
  HANDLE_ERROR();
  printf("wrote Z_gt.bin\n");

  /* compute a butterfly factorization of the selected source and
   * target nodes */

  BfQuadtreeLevelIter srcLevelIter, tgtLevelIter;
  numFactors = bfFacHelm2Prepare(
    srcNode, tgtNode, K, &srcLevelIter, &tgtLevelIter);
  if (numFactors == 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  factorization = bfFacHelm2Make(
    &tree, K, &srcLevelIter, &tgtLevelIter, numFactors);
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
    factor = bfMatProductGetFactor(factorization, i);

    BfMatType matType = bfMatGetType(factor);

    BfMatBlock *matBlock = bfMatToMatBlock(factor);
    BfSize numRowBlocks = bfMatBlockGetNumRowBlocks(matBlock);
    BfSize numColBlocks = bfMatBlockGetNumColBlocks(matBlock);
    BfSize numBlocks = bfMatBlockNumBlocks(matBlock);

    char path[1024];
    sprintf(path, "factor%lu", i);

    chdir(path);

    FILE *fp = fopen("info.txt", "w");
    fprintf(fp, "numBlockRows %lu\n", numRowBlocks);
    fprintf(fp, "numBlockCols %lu\n", numColBlocks);
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
    fwrite(matBlock->rowOffset, sizeof(BfSize), numRowBlocks + 1, fp);
    fclose(fp);

    fp = fopen("colOffset.bin", "w");
    fwrite(matBlock->colOffset, sizeof(BfSize), numColBlocks + 1, fp);
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

      if (i > 0) {
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
  bfComplexRandn(3, q);

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

  bfMatDenseComplexDeinitAndDealloc(&Q);
  bfMatDenseComplexDeinitAndDealloc(&Z_gt);
  bfFreePoints2(&srcPts);
  bfFreePoints2(&tgtPts);
  bfFreeQuadtree(&tree);
  bfFreePoints2(&points);
}
