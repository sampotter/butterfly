#include <stdio.h>

#include "error_macros.h"
#include "fac.h"
#include "mat_block_dense.h"
#include "mat_dense_complex.h"
#include "quadtree.h"
#include "rand.h"

int main(int argc, char const *argv[]) {
  if (argc != 2) {
    printf("usage: %s <points.bin>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  BEGIN_ERROR_HANDLING();

  BfPoints2 points;
  bfReadPoints2FromFile(argv[1], &points);
  HANDLE_ERROR();
  printf("read points from %s\n", argv[1]);

  BfQuadtree tree;
  bfInitQuadtreeFromPoints(&tree, &points);
  HANDLE_ERROR();
  printf("built quadtree\n");

  BfReal K = 3000;

  BfMatBlockDense *A = bfFacHelm2MakeMultilevel(&tree, K);
  printf("built HODBF matrix\n");

  BfMatDenseComplex *x = bfMatDenseComplexNew();
  bfMatDenseComplexInit(x, points.size, 1);
  bfComplexRandn(points.size, x->data);

  BfMat *y = bfMatMul(bfMatBlockDenseGetMatPtr(A), bfMatDenseComplexGetMatPtr(x));
  printf("multiplied with random vector\n");

  END_ERROR_HANDLING() {}

  bfMatDeinitAndDelete(&y);
  bfMatDenseComplexDeinitAndDelete(&x);
  // bfMatBlockDenseDeinitAndDelete(&mat);
  bfFreeQuadtree(&tree);
  bfFreePoints2(&points);
}
