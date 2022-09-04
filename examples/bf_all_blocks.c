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
#include <bf/quadtree.h>
#include <bf/rand.h>
#include <bf/util.h>

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

static BfMatBlockCoo *unravel(BfMatBlockDense *A) {
  BfMatBlockCoo *coo = bfMatBlockCooNew();

  for (size_t i = 0; i < bfMatBlockGetNumRowBlocks(&A->super); ++i) {
    for (size_t j = 0; j < bfMatBlockGetNumColBlocks(&A->super); ++j) {
      BfMat const *block = bfMatBlockDenseGetBlock(A, i, j);
      BfMatType matType = bfMatGetType(block);
      printf("%lu, %lu: %d\n", i, j, matType);
    }
  }

  return coo;
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

  BfMatDenseComplex *A_true = bfGetHelm2KernelMatrix(&points, &points, K);
  printf("computed dense kernel matrix [%0.2fs]\n", bfToc());

  assert(bfMatDenseComplexIsFinite(A_true));

  BfMatBlockDense *A = bfFacHelm2MakeMultilevel(&tree, K);
  printf("assembled HODBF matrix [%0.2fs]\n", bfToc());

  BfMatBlockCoo *A_unravelled = unravel(A);
  (void)A_unravelled;

  BfMatDenseComplex *x = bfMatDenseComplexNew();
  bfMatDenseComplexInit(x, points.size, 1);
  bfSeed(0);
  bfComplexRandn(points.size, x->data);
  printf("set up random RHS [%0.2fs]\n", bfToc());

  assert(bfMatDenseComplexIsFinite(x));

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
  bfMatBlockDenseDeinitAndDealloc(&A);
  bfMatDenseComplexDeinitAndDealloc(&A_true);
  bfFreeQuadtree(&tree);
  bfFreePoints2(&points);
}
