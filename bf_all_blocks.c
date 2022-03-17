#include <stdio.h>

#include "error_macros.h"
#include "fac.h"
#include "mat_block_dense.h"
#include "mat_dense_complex.h"
#include "quadtree.h"
#include "rand.h"
#include "util.h"

int main(int argc, char const *argv[]) {
  if (argc != 2) {
    printf("usage: %s <points.bin>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  BEGIN_ERROR_HANDLING();

  bfToc();

  BfPoints2 points;
  bfReadPoints2FromFile(argv[1], &points);
  HANDLE_ERROR();
  printf("read points from %s [%0.2fs]\n", argv[1], bfToc());

  BfQuadtree tree;
  bfInitQuadtreeFromPoints(&tree, &points);
  HANDLE_ERROR();
  printf("built quadtree [%0.2fs]\n", bfToc());

  BfReal K = 3000;

  BfMatBlockDense *A = bfFacHelm2MakeMultilevel(&tree, K);
  printf("built HODBF matrix [%0.2fs]\n", bfToc());

  BfMatDenseComplex *x = bfMatDenseComplexNew();
  bfMatDenseComplexInit(x, points.size, 1);
  bfComplexRandn(points.size, x->data);

  BfMat *y = bfMatMul(bfMatBlockDenseGetMatPtr(A), bfMatDenseComplexGetMatPtr(x));
  printf("multiplied with random vector [%0.2fs]\n", bfToc());

  END_ERROR_HANDLING() {}

  bfMatDeinitAndDelete(&y);
  bfMatDenseComplexDeinitAndDelete(&x);
  // bfMatBlockDenseDeinitAndDelete(&mat);
  bfFreeQuadtree(&tree);
  bfFreePoints2(&points);
}
