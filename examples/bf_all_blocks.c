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
#include <bf/points.h>
#include <bf/quadtree.h>
#include <bf/rand.h>
#include <bf/util.h>

int main(int argc, char const *argv[]) {
  if (argc != 4) {
    printf("usage: %s <K> <points.bin> <blocks.txt>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  BEGIN_ERROR_HANDLING();

  char const *K_str = argv[1];
  char const *points_path_str = argv[2];
  char const *blocks_path_str = argv[3];

  bfToc();

  BfPoints2 points;
  bfReadPoints2FromFile(points_path_str, &points);
  HANDLE_ERROR();
  printf("read points from %s [%0.2fs]\n", points_path_str, bfToc());

  BfQuadtree tree;
  bfInitQuadtreeFromPoints(&tree, &points);
  HANDLE_ERROR();
  printf("built quadtree [%0.2fs]\n", bfToc());

  BfReal K = atoi(K_str);

  BfMatDenseComplex *A_true = bfGetHelm2KernelMatrix(&points, &points, K);
  printf("computed dense kernel matrix [%0.2fs]\n", bfToc());

  assert(bfMatDenseComplexIsFinite(A_true));

  BfMat *A = bfFacHelm2MakeMultilevel(&tree, K, BF_LAYER_POTENTIAL_SINGLE);
  printf("assembled HODBF matrix [%0.2fs]\n", bfToc());

  BfMatDenseComplex *x = bfMatDenseComplexNew();
  bfMatDenseComplexInit(x, points.size, 1);
  bfSeed(0);
  bfComplexRandn(points.size, x->data);
  printf("set up random RHS [%0.2fs]\n", bfToc());

  assert(bfMatDenseComplexIsFinite(x));

  BfMat *y_true = bfMatMul(bfMatDenseComplexToMat(A_true), bfMatDenseComplexToMat(x));
  printf("multiplied with dense kernel matrix [%0.2fs]\n", bfToc());

  BfMat *y = bfMatMul(A, bfMatDenseComplexToMat(x));
  printf("multiplied with HODBF matrix [%0.2fs]\n", bfToc());

  FILE *fp = fopen(blocks_path_str, "w");
  bfPrintBlocks(A, 2, fp);
  fclose(fp);
  printf("wrote blocks to %s [%0.2fs]\n", blocks_path_str, bfToc());

  END_ERROR_HANDLING() {}

  bfMatDelete(&y);
  bfMatDelete(&y_true);
  bfMatDenseComplexDeinitAndDealloc(&x);
  bfMatDelete(&A);
  bfMatDenseComplexDeinitAndDealloc(&A_true);
  bfFreeQuadtree(&tree);
  bfFreePoints2(&points);
}
