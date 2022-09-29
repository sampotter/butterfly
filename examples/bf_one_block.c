#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
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

void print_usage_and_exit(char const *argv0) {
  printf("usage: %s <K> <src_depth> <src_index> <tgt_depth> <tgt_index> "
         "<points.bin> [layerpot] [normals.bin]\n"
         "\n"
         "where:\n"
         "  layerpot is one of: S, Sp, D\n"
         "\n"
         "NOTE: double-layer potential (D) not yet implemented\n", argv0);
  exit(EXIT_FAILURE);
}

int main(int argc, char const *argv[]) {
  if (argc < 7)
    print_usage_and_exit(argv[0]);

  char const *layerPotStr = argc >= 8 ? argv[7] : "S";
  if (!strcmp(layerPotStr, "Sp") && argc != 9) {
    printf("usage: %s <K> <src_depth> <src_index> <tgt_depth> <tgt_index> "
           "<points.bin> Sp <normals.bin>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  BEGIN_ERROR_HANDLING();

  BfLayerPotential layerPot = BF_LAYER_POTENTIAL_UNKNOWN;
  if (!strcmp(layerPotStr, "S")) {
    layerPot = BF_LAYER_POTENTIAL_SINGLE;
    printf("using single-layer potential\n");
  } else if (!strcmp(layerPotStr, "Sp")) {
    layerPot = BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_SINGLE;
    printf("using PV of normal derivative of single-layer potential\n");
  } else if (!strcmp(layerPotStr, "D")) {
    printf("ERROR: double layer not yet implemented\n");
    exit(EXIT_FAILURE);
  } else {
    print_usage_and_exit(argv[0]);
  }

  BfSize numFactors = BF_SIZE_BAD_VALUE;
  BfMatProduct *factorization = NULL;
  BfMat *factor = NULL;

  BfReal K = atoi(argv[1]);

  BfPoints2 points;
  bfReadPoints2FromFile(argv[6], &points);
  HANDLE_ERROR();
  printf("read points from %s\n", argv[2]);

  BfVectors2 tgtNormals, *tgtNormalsPtr = NULL;
  if (layerPot != BF_LAYER_POTENTIAL_SINGLE) {
    bfReadVectors2FromFile(argv[8], &tgtNormals);
    HANDLE_ERROR();
    printf("read unit normals from %s\n", argv[8]);
    tgtNormalsPtr = &tgtNormals;
  }

  BfQuadtree tree;
  bfInitQuadtreeFromPoints(&tree, &points, tgtNormalsPtr);
  HANDLE_ERROR();
  puts("built quadtree");

  bfQuadtreeSaveBoxesToTextFile(&tree, "quadtree.txt");
  puts("saved quadtree boxes to quadtree.txt");

  bfSeed(1234);

  /* Get source and target nodes from quadtree and check that their
   * indices are OK */

  BfSize src_depth = atoi(argv[2]);
  BfSize srcNodeIndex = atoi(argv[3]);
  BfSize tgt_depth = atoi(argv[4]);
  BfSize tgtNodeIndex = atoi(argv[5]);

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

  if (layerPot != BF_LAYER_POTENTIAL_SINGLE) {
    BfVectors2 tgtNormals_;

    bfQuadtreeNodeGetUnitNormals(&tree, tgtNode, &tgtNormals_);
    HANDLE_ERROR();

    bfSaveVectors2(&tgtNormals_, "tgtNormals.bin");
    HANDLE_ERROR();

    puts("wrote target unit normals to tgtNormals.bin");
  }

  BfMat *Z_gt = bfGetHelm2KernelMatrix(&srcPts, &tgtPts, tgtNormalsPtr, K, layerPot);
  HANDLE_ERROR();

  printf("computed groundtruth subblock of kernel matrix\n");

  bfMatDenseComplexSave(Z_gt, "Z_gt.bin");
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
    &tree, K, layerPot, &srcLevelIter, &tgtLevelIter, numFactors);
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

    BfType type = bfMatGetType(factor);

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
    if (type == BF_TYPE_MAT_BLOCK_COO)
      fwrite(((BfMatBlockCoo *)factor)->rowInd, sizeof(BfSize), numBlocks, fp);
    else if (type == BF_TYPE_MAT_BLOCK_DIAG)
      for (BfSize j = 0; j < numBlocks; ++j)
        fwrite(&j, sizeof(BfSize), 1, fp);
    else
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);
    fclose(fp);

    fp = fopen("colInd.bin", "w");
    if (type == BF_TYPE_MAT_BLOCK_COO)
      fwrite(((BfMatBlockCoo *)factor)->colInd, sizeof(BfSize), numBlocks, fp);
    else if (type == BF_TYPE_MAT_BLOCK_DIAG)
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
      BfFacAux *facAux = matBlock->block[j]->aux;

      BfPoints2 *srcChildPts = &facAux->srcPts[0];
      BfPoints2 *srcPts = &facAux->srcPts[1];
      BfPoints2 *tgtChildPts = &facAux->tgtPts;
      BfVectors2 *tgtNormals = &facAux->tgtNormals;

      sprintf(filename, "srcPtsOrig%lu.bin", j);
      bfSavePoints2(srcChildPts, filename);

      if (i > 0) {
        sprintf(filename, "srcPtsEquiv%lu.bin", j);
        bfSavePoints2(srcPts, filename);
      }

      sprintf(filename, "tgtPts%lu.bin", j);
      bfSavePoints2(tgtChildPts, filename);

      if (layerPot != BF_LAYER_POTENTIAL_SINGLE) {
        sprintf(filename, "tgtNormals%lu.bin", j);
        bfSaveVectors2(tgtNormals, filename);
      }
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

  BfMat const *Z_mat = bfMatProductToMat(factorization);
  BfMat const *Q_mat = bfMatDenseComplexToMat(Q);

  BfMat *V_gt = bfMatMul(Z_gt, Q_mat);
  bfMatSave(V_gt, "V_gt.bin");

  BfMat *V = bfMatMul(Z_mat, Q_mat);
  bfMatSave(V, "V.bin");

  /* write simulation info to a text file */

  FILE *fp = fopen("info.txt", "w");
  fprintf(fp, "layerPot %s\n", layerPotStr);
  fprintf(fp, "numSrcPts %lu\n", srcPts.size);
  fprintf(fp, "numTgtPts %lu\n", tgtPts.size);
  fprintf(fp, "numTgtPts %lu\n", tgtPts.size);
  fprintf(fp, "numFactors %lu\n", numFactors);
  fprintf(fp, "K %g\n", K);
  fclose(fp);
  puts("wrote simulation information to info.txt");

  END_ERROR_HANDLING() {
    puts("error!");
  }

  bfMatDenseComplexDeinitAndDealloc(&Q);
  bfMatDelete(&Z_gt);
  bfFreePoints2(&srcPts);
  bfFreePoints2(&tgtPts);
  bfFreeQuadtree(&tree);
  bfFreePoints2(&points);
}
