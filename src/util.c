#include <bf/util.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_block_coo.h>
#include <bf/mat_block_diag.h>

#include <time.h>

BfReal bfToc() {
  static clock_t t1 = 0;
  clock_t t0 = t1;
  t1 = clock();
  return ((double)t1 - (double)t0)/CLOCKS_PER_SEC;
}

void bfSizeSetConstant(BfSize numSizes, BfSize *size, BfSize value) {
  for (BfSize i = 0; i < numSizes; ++i)
    size[i] = value;
}

void bfSizeRunningSum(BfSize numSizes, BfSize *size) {
  for (BfSize i = 1; i < numSizes; ++i)
    size[i] += size[i - 1];
}

static void printBlocksRec(BfMat const *mat,
                           BfSize level, BfSize i0, BfSize j0,
                           FILE *fp) {
  BfType type = bfMatGetType(mat);

  if (type == BF_TYPE_MAT_BLOCK_COO) {
    BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
    BfMatBlockCoo const *matBlockCoo = bfMatConstToMatBlockCooConst(mat);
    BfSize numBlocks = bfMatBlockCooNumBlocks(matBlock);
    for (BfSize k = 0; k < numBlocks; ++k) {
      BfMat const *block = matBlock->block[k];
      BfSize di = matBlock->rowOffset[matBlockCoo->rowInd[k]];
      BfSize dj = matBlock->colOffset[matBlockCoo->colInd[k]];
      printBlocksRec(block, level + 1, i0 + di, j0 + dj, fp);
    }
  }

  else if (type == BF_TYPE_MAT_BLOCK_DENSE) {
    BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
    BfSize numRowBlocks = bfMatBlockGetNumRowBlocks(matBlock);
    BfSize numColBlocks = bfMatBlockGetNumColBlocks(matBlock);
    for (BfSize k = 0; k < numRowBlocks; ++k) {
      for (BfSize l = 0; l < numColBlocks; ++l) {
        BfMat const *block = matBlock->block[k*numColBlocks + l];
        BfSize di = matBlock->rowOffset[k];
        BfSize dj = matBlock->colOffset[l];
        printBlocksRec(block, level + 1, i0 + di, j0 + dj, fp);
      }
    }
  }

  else if (type == BF_TYPE_MAT_BLOCK_DIAG) {
    BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
    BfSize numBlocks = bfMatBlockDiagNumBlocks(matBlock);
    for (BfSize k = 0; k < numBlocks; ++k) {
      BfMat const *block = matBlock->block[k];
      BfSize di = matBlock->rowOffset[k];
      BfSize dj = matBlock->colOffset[k];
      printBlocksRec(block, level + 1, i0 + di, j0 + dj, fp);
    }
  }

  else {
    BfSize i1 = i0 + bfMatGetNumRows(mat);
    BfSize j1 = j0 + bfMatGetNumCols(mat);
    fprintf(fp, "%lu %lu %lu %lu %lu %d\n", level, i0, i1, j0, j1, type);
  }
}

void bfPrintBlocks(BfMat const *mat, BfSize level, FILE *fp) {
  printBlocksRec(mat, level, 0, 0, fp);
}

BfSize bfGetFileSizeInBytes(char const *path) {
  BEGIN_ERROR_HANDLING();

  BfSize numBytes = BF_SIZE_BAD_VALUE;

  FILE *fp = fopen(path, "rb");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  fseek(fp, SEEK_END, 0);
  numBytes = ftell(fp);

  END_ERROR_HANDLING() {}

  fclose(fp);

  return numBytes;
}

void bfReadFileToMemory(char const *path, BfSize numBytes, BfByte *ptr) {
  BEGIN_ERROR_HANDLING();

  FILE *fp = fopen(path, "rb");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  fread(ptr, numBytes, sizeof(BfByte), fp);

  END_ERROR_HANDLING() {}

  fclose(fp);
}
