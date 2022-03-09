#include <stdio.h>

#include "error_macros.h"
#include "fac.h"
#include "mat_block_dense.h"
#include "quadtree.h"



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
  puts("built quadtree");

  BfReal K = 3000;

  BfMatBlockDense *mat = bfFacHelm2MakeMultilevel(&tree, K);
  (void)mat;

  END_ERROR_HANDLING() {}

  // bfMatBlockDenseDeinitAndDelete(&mat);
  bfFreeQuadtree(&tree);
  bfFreePoints2(&points);
}
