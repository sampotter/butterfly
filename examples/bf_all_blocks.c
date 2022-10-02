#include <assert.h>
#include <stdio.h>

#include <bf/error_macros.h>
#include <bf/fac.h>
#include <bf/helm2.h>
#include <bf/layer_pot.h>
#include <bf/mat_block_coo.h>
#include <bf/mat_block_dense.h>
#include <bf/mat_block_diag.h>
#include <bf/mat_dense_complex.h>
#include <bf/mat_util.h>
#include <bf/points.h>
#include <bf/quadtree.h>
#include <bf/rand.h>
#include <bf/util.h>
#include <bf/vec_real.h>

int main(int argc, char const *argv[]) {
  if (argc != 4) {
    printf("usage: %s <K> <points.bin> <blocks.txt>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  BEGIN_ERROR_HANDLING();

  char const *K_str = argv[1];
  char const *points_path_str = argv[2];
  char const *blocks_path_str = argv[3];

  BfReal K = atoi(K_str);

  BfMat *pointsMat = bfMatFromFile(points_path_str, -1, 2, BF_DTYPE_REAL);
  HANDLE_ERROR();
  printf("read points from %s [%0.2fs]\n", points_path_str, bfToc());

  BfPoints2 const *points = bfPoints2ConstViewFromMat(pointsMat);

  BfSize numPoints = points->size;

  bfToc();
  BfQuadtree tree;
  bfInitQuadtreeFromPoints(&tree, points, NULL);
  HANDLE_ERROR();
  printf("built quadtree [%0.2fs]\n", bfToc());

  BfPerm revPerm = bfPermGetReversePerm(&tree.perm);

  BfMat *pointsPermMat = bfMatCopy(pointsMat);
  bfMatPermuteRows(pointsPermMat, &revPerm);

  BfPoints2 const *pointsPerm = bfPoints2ConstViewFromMat(pointsPermMat);

  bfToc();

  BfMat *A_dense = bfGetHelm2KernelMatrix(
    pointsPerm, pointsPerm, NULL, K, BF_LAYER_POTENTIAL_SINGLE);
  printf("computed dense kernel matrix [%0.2fs]\n", bfToc());

  BfMat *A_BF = bfFacHelm2MakeMultilevel(&tree, K, BF_LAYER_POTENTIAL_SINGLE);
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

  BfMat *y_BF = bfMatMul(A_BF, x);
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

  END_ERROR_HANDLING() {}

  bfMatDelete(&y_BF);
  bfMatDelete(&y_dense);
  bfMatDelete(&x);
  bfMatDelete(&A_BF);
  bfMatDelete(&A_dense);
  bfFreeQuadtree(&tree);
  bfMatDelete(&pointsPermMat);
  bfMatDelete(&pointsMat);
}
