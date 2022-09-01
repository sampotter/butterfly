#include <assert.h>
#include <stdio.h>

#include <openblas/cblas.h>

#include "error_macros.h"
#include "fac.h"
#include "helm2.h"
#include "mat_block_coo.h"
#include "mat_block_dense.h"
#include "mat_block_diag.h"
#include "mat_dense_complex.h"
#include "quadrature.h"
#include "quadtree.h"
#include "rand.h"
#include "util.h"

/* TODO: should make this functionality a part of the library... need
 * some mappers for the hierarchical matrix tree... also need to wrap
 * up the hierarchical matrix as a separate ADT to make the interface
 * clearer */
void printBlocks(BfMat const *mat,FILE *fp,BfSize i0,BfSize j0,BfSize level) {
  BfMatType matType = bfMatGetType(mat);

  if (matType == BF_MAT_TYPE_BLOCK_COO) {
    BfMatBlockCoo const *sub = (BfMatBlockCoo const *)mat;
    BfSize numBlocks = bfMatBlockCooNumBlocks(sub);
    for (BfSize k = 0; k < numBlocks; ++k) {
      BfMat const *block = sub->super.block[k];
      BfSize di = sub->super.rowOffset[sub->rowInd[k]];
      BfSize dj = sub->super.colOffset[sub->colInd[k]];
      printBlocks(block, fp, i0 + di, j0 + dj, level + 1);
    }
  }

  else if (matType == BF_MAT_TYPE_BLOCK_DENSE) {
    BfMatBlockDense const *sub = (BfMatBlockDense const *)mat;
    BfSize numRowBlocks = bfMatBlockGetNumRowBlocks(&sub->super);
    BfSize numColBlocks = bfMatBlockGetNumColBlocks(&sub->super);
    for (BfSize k = 0; k < numRowBlocks; ++k) {
      for (BfSize l = 0; l < numColBlocks; ++l) {
        BfMat const *block = sub->super.block[k*numColBlocks + l];
        BfSize di = sub->super.rowOffset[k];
        BfSize dj = sub->super.colOffset[l];
        printBlocks(block, fp, i0 + di, j0 + dj, level + 1);
      }
    }
  }

  else if (matType == BF_MAT_TYPE_BLOCK_DIAG) {
    BfMatBlockDiag const *sub = (BfMatBlockDiag const *)mat;
    BfSize numBlocks = bfMatBlockDiagNumBlocks(sub);
    for (BfSize k = 0; k < numBlocks; ++k) {
      BfMat const *block = sub->super.block[k];
      BfSize di = sub->super.rowOffset[k];
      BfSize dj = sub->super.colOffset[k];
      printBlocks(block, fp, i0 + di, j0 + dj, level + 1);
    }
  }

  else {
    BfSize i1 = i0 + bfMatGetNumRows(mat);
    BfSize j1 = j0 + bfMatGetNumCols(mat);
    fprintf(fp, "%lu %lu %lu %lu %lu %d\n", level, i0, i1, j0, j1, matType);
  }
}

int main(int argc, char const *argv[]) {
  if (argc != 4) {
    printf("usage: %s <K> <points.bin> <blocks.txt>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  /* Do these comparisons on a single thread for now. */
  openblas_set_num_threads(1);

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

  // TODO: modifications to make:
  // - need Sp matrix
  //   + need to check that it works with a priori rank estimate
  // - A_true -> Nystrom system matrix using Kapur-Rokhlin quadrature
  // - A -> butterfly-factorized Nystrom matrix
  // - "well-separated condition" -> "well-separated" +
  //                                 boxes *don't* contain corrected nodes
  //                                 (they must only belong to dense near-field blocks
  //                                  or diagonal blocks!)
  // - need to write a function that will compute the Nystrom system matrix
  //   + then just use that instead of bfGetHelm2KernelMatrix when we bottom out
  // - bfFacHelm2MakeMultilevel needs to be modified
  // - set up ellipse problem as before
  // - compute LU factorization of A_true and solve
  // - use GMRES to solve using A_true
  // - use GMRES to solve using A
  // - need to make sure I'm permuting and unpermuting before and
  //   after multiplying...

  BfMatDenseComplex *A_true = bfGetHelm2KernelMatrix(&points, &points, K);
  printf("computed dense kernel matrix [%0.2fs]\n", bfToc());

  BfMatBlockDense *A = bfFacHelm2MakeMultilevel(&tree, K);
  printf("assembled HODBF matrix [%0.2fs]\n", bfToc());

  BfMatDenseComplex *x = bfMatDenseComplexNew();
  bfMatDenseComplexInit(x, points.size, 1);
  bfComplexRandn(points.size, x->data);
  printf("set up random RHS [%0.2fs]\n", bfToc());

  BfMat *y_true = bfMatMul(bfMatDenseComplexGetMatPtr(A_true), bfMatDenseComplexGetMatPtr(x));
  printf("multiplied with dense kernel matrix [%0.2fs]\n", bfToc());

  BfMat *y = bfMatMul(bfMatBlockDenseGetMatPtr(A), bfMatDenseComplexGetMatPtr(x));
  printf("multiplied with HODBF matrix [%0.2fs]\n", bfToc());

  FILE *fp = fopen(blocks_path_str, "w");
  printBlocks(bfMatBlockDenseGetMatConstPtr(A), fp, 0, 0, 2);
  fclose(fp);
  printf("wrote blocks to %s [%0.2fs]\n", blocks_path_str, bfToc());

  END_ERROR_HANDLING() {}

  bfMatDeinitAndDelete(&y);
  bfMatDeinitAndDelete(&y_true);
  bfMatDenseComplexDeinitAndDelete(&x);
  // bfMatBlockDenseDeinitAndDelete(&mat);
  bfMatDenseComplexDeinitAndDelete(&A_true);
  bfFreeQuadtree(&tree);
  bfFreePoints2(&points);
}
