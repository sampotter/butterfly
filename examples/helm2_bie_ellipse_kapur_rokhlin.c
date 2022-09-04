#include <assert.h>
#include <stdio.h>

#include <openblas/cblas.h>

#include <bf/error_macros.h>
#include <bf/fac.h>
#include <bf/helm2.h>
#include <bf/mat_block_coo.h>
#include <bf/mat_block_dense.h>
#include <bf/mat_block_diag.h>
#include <bf/mat_dense_complex.h>
#include <bf/quadrature.h>
#include <bf/quadtree.h>
#include <bf/rand.h>
#include <bf/util.h>

void printBlocks(BfMat const *mat,FILE *fp,BfSize i0,BfSize j0,BfSize level) {
  BfMatType matType = bfMatGetType(mat);

  if (matType == BF_MAT_TYPE_BLOCK_COO) {
    BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
    BfMatBlockCoo const *matBlockCoo = bfMatConstToMatBlockCooConst(mat);
    BfSize numBlocks = bfMatBlockCooNumBlocks(matBlock);
    for (BfSize k = 0; k < numBlocks; ++k) {
      BfMat const *block = matBlock->block[k];
      BfSize di = matBlock->rowOffset[matBlockCoo->rowInd[k]];
      BfSize dj = matBlock->colOffset[matBlockCoo->colInd[k]];
      printBlocks(block, fp, i0 + di, j0 + dj, level + 1);
    }
  }

  else if (matType == BF_MAT_TYPE_BLOCK_DENSE) {
    BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
    BfSize numRowBlocks = bfMatBlockGetNumRowBlocks(matBlock);
    BfSize numColBlocks = bfMatBlockGetNumColBlocks(matBlock);
    for (BfSize k = 0; k < numRowBlocks; ++k) {
      for (BfSize l = 0; l < numColBlocks; ++l) {
        BfMat const *block = matBlock->block[k*numColBlocks + l];
        BfSize di = matBlock->rowOffset[k];
        BfSize dj = matBlock->colOffset[l];
        printBlocks(block, fp, i0 + di, j0 + dj, level + 1);
      }
    }
  }

  else if (matType == BF_MAT_TYPE_BLOCK_DIAG) {
    BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
    BfSize numBlocks = bfMatBlockDiagNumBlocks(matBlock);
    for (BfSize k = 0; k < numBlocks; ++k) {
      BfMat const *block = matBlock->block[k];
      BfSize di = matBlock->rowOffset[k];
      BfSize dj = matBlock->colOffset[k];
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

  BfMat *y_true = bfMatMul(bfMatDenseComplexToMat(A_true), bfMatDenseComplexToMat(x));
  printf("multiplied with dense kernel matrix [%0.2fs]\n", bfToc());

  BfMat *y = bfMatMul(bfMatBlockDenseToMat(A), bfMatDenseComplexToMat(x));
  printf("multiplied with HODBF matrix [%0.2fs]\n", bfToc());

  FILE *fp = fopen(blocks_path_str, "w");
  printBlocks(bfMatBlockDenseConstToMatConst(A), fp, 0, 0, 2);
  fclose(fp);
  printf("wrote blocks to %s [%0.2fs]\n", blocks_path_str, bfToc());

  END_ERROR_HANDLING() {}

  bfMatDelete(&y);
  bfMatDelete(&y_true);
  bfMatDenseComplexDeinitAndDealloc(&x);
  // bfMatBlockDenseDeinitAndDealloc(&mat);
  bfMatDenseComplexDeinitAndDealloc(&A_true);
  bfFreeQuadtree(&tree);
  bfFreePoints2(&points);
}
