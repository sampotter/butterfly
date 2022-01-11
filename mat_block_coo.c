#include "mat_block_coo.h"

#include <stdlib.h>

BfMatBlockCoo *bfMatBlockCooNew() {
  return malloc(sizeof(BfMatBlockCoo));
}

void bfMatBlockCooDeinit(BfMatBlockCoo *mat) {
  free(mat->rowInd);
  free(mat->colInd);
  free(mat->rowOffset);
  free(mat->colOffset);

  for (BfSize i = 0; i < mat->numBlocks; ++i)
    bfMatDeinitAndDelete(&mat->block[i]);
  free(mat->block);
}

void bfMatBlockCooDelete(BfMatBlockCoo **mat) {
  free(*mat);
  *mat = NULL;
}

BfMat *bfMatBlockCooGetMatPtr(BfMatBlockCoo *mat) {
  return &mat->super;
}

// static BfSize getNumRows(BfFactor const *factor) {
//   BfMatBlockCoo *mat = factor->mat;
//   return mat->rowOffset[mat->super.numBlockRows];
// }

// static BfSize getNumCols(BfFactor const *factor) {
//   BfMatBlockCoo *mat = factor->mat;
//   return mat->colOffset[mat->super.numBlockCols];
// }

// void bfMulFac(BfFactor const *factor, BfMat const *X, BfMat *Y) {
//   BEGIN_ERROR_HANDLING();

//   BfSize p = getNumCols(factor);
//   if (X->numRows != p)
//     RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

//   BfSize m = getNumRows(factor), n = X->numCols;

//   *Y = bfGetUninitializedMat();
//   bfMatZeros(Y, BF_DTYPE_COMPLEX, m, n);
//   HANDLE_ERROR();

//   BfMat tmp;

//   for (BfSize k = 0; k < factor->numBlocks; ++k) {
//     BfSize i0 = factor->rowOffset[factor->rowInd[k]];
//     BfSize i1 = factor->rowOffset[factor->rowInd[k] + 1];

//     BfSize j0 = factor->colOffset[factor->colInd[k]];
//     BfSize j1 = factor->colOffset[factor->colInd[k] + 1];

//     BfMat Xrows = bfGetMatRowRange(X, j0, j1);

//     tmp = bfGetUninitializedMat();
//     bfMatMul(&factor->block[k], &Xrows, &tmp);
//     HANDLE_ERROR();

//     BfMat Yrows = bfGetMatRowRange(Y, i0, i1);
//     bfMatAddInplace(&Yrows, &tmp);
//     HANDLE_ERROR();

//     bfMatFree(&tmp);
//   }

//   END_ERROR_HANDLING() {
//     bfMatFree(&tmp);
//     bfMatFree(Y);
//   }
// }
