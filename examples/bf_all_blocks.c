#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <bf/error_macros.h>
#include <bf/fac.h>
#include <bf/helm2.h>
#include <bf/layer_pot.h>
#include <bf/mat_block_coo.h>
#include <bf/mat_block_dense.h>
#include <bf/mat_block_diag.h>
#include <bf/mat_dense_complex.h>
#include <bf/points.h>
#include <bf/quadtree.h>
#include <bf/rand.h>
#include <bf/util.h>
#include <bf/vec_real.h>

/* Set PERM_DENSE equal to 1 to permute A_dense instead of permuting
 * the righthand side. This is also required if we want to compare the
 * kernel matrices themselves. */
#define PERM_DENSE 1

/* Right multiply the butterflied kernel matrices with the identity
 * matrix to do a comparison between groundtruth dense kernel matrix
 * and the butterflied version. */
#define COMPARE_KERNEL_MATRICES 1

void print_usage_and_exit(char const *argv0) {
  printf("usage: %s <K> <points.bin> <blocks.txt> [layerpot] [normals.bin]\n"
         "\n"
         "where:\n"
         "  layerpot is one of: S, Sp, D\n"
         "\n"
         "NOTE: double-layer potential (D) not yet implemented\n", argv0);
  exit(EXIT_FAILURE);
}

int main(int argc, char const *argv[]) {
  if (argc < 4)
    print_usage_and_exit(argv[0]);

  char const *K_str = argv[1];
  char const *points_path_str = argv[2];
  char const *blocks_path_str = argv[3];
  char const *layerPotStr = argc >= 5 ? argv[4] : "S";

  if (!strcmp(layerPotStr, "Sp") && argc != 6) {
    printf("usage: %s <K> <points.bin> <blocks.txt> %s <normals.bin>\n",
           argv[0], layerPotStr);
    exit(EXIT_FAILURE);
  }

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

  char const *normals_path_str = !strcmp(layerPotStr, "Sp") ? argv[5] : NULL;

  BEGIN_ERROR_HANDLING();

  BfReal K = atoi(K_str);

  BfMat *pointsMat = bfMatFromFile(points_path_str, -1, 2, BF_DTYPE_REAL);
  HANDLE_ERROR();
  printf("read points from %s [%0.2fs]\n", points_path_str, bfToc());

  BfPoints2 const *points = bfPoints2ConstViewFromMat(pointsMat);

  BfSize numPoints = points->size;

  BfVectors2 normals, *normalsPtr = NULL;
  if (layerPot != BF_LAYER_POTENTIAL_SINGLE) {
    bfReadVectors2FromFile(normals_path_str, &normals);
    HANDLE_ERROR();
    printf("read unit normals from %s\n", normals_path_str);
    normalsPtr = &normals;
  }

  bfToc();
  BfQuadtree quadtree;
  bfQuadtreeInit(&quadtree, points, normalsPtr);
  HANDLE_ERROR();
  printf("built quadtree [%0.2fs]\n", bfToc());

  BfTree *tree = bfQuadtreeToTree(&quadtree);

  BfPerm const *perm = bfTreeGetPermConst(tree);
  BfPerm revPerm = bfPermGetReversePerm(perm);

#if PERM_DENSE
  BfMat *pointsPermMat = bfMatCopy(pointsMat);
  bfMatPermuteRows(pointsPermMat, &revPerm);

  BfPoints2 const *pointsPerm = bfPoints2ConstViewFromMat(pointsPermMat);

  BfMat *normalsPermMat = NULL;
  BfVectors2 const *normalsPermPtr = NULL;
  if (layerPot != BF_LAYER_POTENTIAL_SINGLE) {
    normalsPermMat = bfMatFromFile(normals_path_str, -1, 2, BF_DTYPE_REAL);
    HANDLE_ERROR();
    bfMatPermuteRows(normalsPermMat, &revPerm);
    normalsPermPtr = bfVectors2ConstViewFromMat(normalsPermMat);
  }

  bfToc();
  BfMat *A_dense = bfGetHelm2KernelMatrix(
    pointsPerm, pointsPerm, NULL, normalsPermPtr, K, layerPot, NULL, NULL);
  printf("computed dense kernel matrix [%0.2fs]\n", bfToc());
#else
  bfToc();
  BfMat *A_dense = bfGetHelm2KernelMatrix(points, points, normalsPtr, K, layerPot);
  printf("computed dense kernel matrix [%0.2fs]\n", bfToc());
#endif

  BfMat *A_BF = bfFacHelm2MakeMultilevel(&quadtree, K, layerPot, NULL, NULL);
  printf("assembled HODBF matrix [%0.2fs]\n", bfToc());

  BfMat *x; {
    BfMatDenseComplex *_ = bfMatDenseComplexNew();
    bfMatDenseComplexInit(_, numPoints, 1);
    bfSeed(0);
    bfComplexRandn(numPoints, _->data);
    x = bfMatDenseComplexToMat(_);
  }
  printf("set up random righthand side [%0.2fs]\n", bfToc());

  BfMat *y_dense = bfMatMul(A_dense, x);
  printf("multiplied with dense kernel matrix [%0.2fs]\n", bfToc());

#if PERM_DENSE
  BfMat *y_BF = bfMatMul(A_BF, x);
#else
  BfMat *x_perm = bfMatCopy(x);
  bfMatPermuteRows(x_perm, &revPerm);
  BfMat *y_BF = bfMatMul(A_BF, x_perm);
  bfMatPermuteRows(y_BF, &tree.perm);
#endif
  printf("multiplied with HODBF matrix [%0.2fs]\n", bfToc());

  BfVec *err2 = bfMatColDists(y_dense, y_BF);
  BfVec *y_dense_norm = bfMatColNorms(y_dense);
  BfReal rel_err_l2 =
    *bfVecToVecReal(err2)->data / *bfVecToVecReal(y_dense_norm)->data;
  printf("relative l2 error: %g [%0.2fs]\n", rel_err_l2, bfToc());

  FILE *fp = fopen(blocks_path_str, "w");
  bfPrintBlocks(A_BF, 2, fp);
  fclose(fp);
  printf("wrote blocks to %s [%0.2fs]\n", blocks_path_str, bfToc());

  bfMatSave(A_dense, "A_dense.bin");
  printf("wrote dense kernel matrix to A_dense.bin [%0.2fs]\n", bfToc());

#if COMPARE_KERNEL_MATRICES
#if !PERM_DENSE
#error "Must have PERM_DENSE == 1 if COMPARE_KERNEL_MATRICES == 1"
#endif
  BfMat *A_BF_dense = bfMatEmptyLike(A_dense, numPoints, numPoints);
  for (BfSize j = 0; j < numPoints; ++j) {
    BfMat *ej; {
      BfMatDenseComplex *_ = bfMatDenseComplexNew();
      bfMatDenseComplexInit(_, numPoints, 1);
      for (BfSize i = 0; i < numPoints; ++i)
        *(_->data + i*_->rowStride) = i == j ? 1 : 0;
      ej = bfMatDenseComplexToMat(_);
    }
    BfMat *aj = bfMatMul(A_BF, ej);
    BfVec *ajVecView = bfMatGetColView(aj, 0);
    bfMatSetCol(A_BF_dense, j, ajVecView);
    bfVecDelete(&ajVecView);
    bfMatDelete(&aj);
  }
  printf("sampled dense version of BF'd kernel matrix [%0.2fs]\n", bfToc());

  bfMatSave(A_BF_dense, "A_BF_dense.bin");
  printf("saved dense version of BF'd kernel matrix to A_BF_dense.bin [%0.2fs]\n", bfToc());
#endif

  END_ERROR_HANDLING() {}

  bfMatDelete(&y_BF);
  bfMatDelete(&y_dense);
  bfMatDelete(&x);
  bfMatDelete(&A_BF);
  bfMatDelete(&A_dense);
  bfQuadtreeDeinit(&quadtree);
#if PERM_DENSE
  bfMatDelete(&pointsPermMat);
#endif
  bfMatDelete(&pointsMat);
}
