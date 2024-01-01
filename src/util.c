#include <bf/util.h>

#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <bf/assert.h>
#include <bf/const.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_block_coo.h>
#include <bf/mat_block_diag.h>
#include <bf/mem.h>

#include "macros.h"

BfReal bfToc() {
  static clock_t t1 = 0;
  clock_t t0 = t1;
  t1 = clock();
  return ((double)t1 - (double)t0)/CLOCKS_PER_SEC;
}

void bfRealArgsort(BfReal const *values, BfSize n, BfSize *index) {
  BF_ERROR_BEGIN();

  for (BfSize i = 0; i < n; ++i)
    index[i] = i;

  /* TODO: just using selection sort here for now... should upgrade to
   * something better later */
  for (BfSize i = 0; i < n - 1; ++i) {
    BfSize k = i;
    for (BfSize j = i + 1; j < n; ++j)
      if (values[index[j]] < values[index[k]])
        k = j;
    SWAP(index[i], index[k]);
  }

  BfPerm perm = {.index = index, .size = n};

  bfPermReverse(&perm);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfSizeSetConstant(BfSize numSizes, BfSize *size, BfSize value) {
  for (BfSize i = 0; i < numSizes; ++i)
    size[i] = value;
}

void bfSizeRunningSum(BfSize numSizes, BfSize *size) {
  for (BfSize i = 1; i < numSizes; ++i)
    size[i] += size[i - 1];
}

bool bfSizeIsPerm(BfSize n, BfSize const *size) {
  BF_ERROR_BEGIN();

  bool isPerm = true;

  bool *found = bfMemAllocAndZero(n, sizeof(bool));
  HANDLE_ERROR();

  for (BfSize i = 0; i < n; ++i) {
    if (size[i] >= n || found[size[i]]) {
      isPerm = false;
      break;
    }
    found[size[i]] = true;
  }

#if BF_DEBUG
  if (isPerm) for (BfSize i = 0; i < n; ++i) BF_ASSERT(found[i]);
#endif

  BF_ERROR_END() {
    BF_DIE();
  }

  bfMemFree(found);

  return isPerm;
}

static void printBlocksRec(BfMat const *mat,
                           BfSize level, BfSize i0, BfSize j0,
                           FILE *fp) {
  BfType type = bfMatGetType(mat);

  if (type == BF_TYPE_MAT_BLOCK_COO) {
    BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
    BfMatBlockCoo const *matBlockCoo = bfMatConstToMatBlockCooConst(mat);
    BfSize numBlocks = bfMatBlockCooNumBlocks(matBlockCoo);
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
    BfMatBlockDiag const *matBlockDiag = bfMatConstToMatBlockDiagConst(mat);
    BfSize numBlocks = bfMatBlockDiagNumBlocks(matBlockDiag);
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
  BF_ERROR_BEGIN();

  BfSize numBytes = BF_SIZE_BAD_VALUE;

  FILE *fp = fopen(path, "rb");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  fseek(fp, SEEK_END, 0);
  numBytes = ftell(fp);

  BF_ERROR_END() {}

  fclose(fp);

  return numBytes;
}

void bfReadFileToMemory(char const *path, BfSize numBytes, BfByte *ptr) {
  BF_ERROR_BEGIN();

  FILE *fp = fopen(path, "rb");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  fread(ptr, numBytes, sizeof(BfByte), fp);

  BF_ERROR_END() {}

  fclose(fp);
}

void bfSort(BfPtr ptr, BfSize n, BfSize size, BfCompar compar, BfPtr aux) {
  qsort_r(ptr, n, size, compar, aux);
}

static int sgn(BfReal x) {
  if (x > 0) {
    return 1;
  } else if (x < 0) {
    return -1;
  } else {
    return 0;
  }
}

bool bfFindZeroOnInterval(BfReal (*f)(BfReal, void *), BfReal a, BfReal b, void *context, BfReal *t) {
  BfReal c, d, fa, fb, fc, fd, dm, df, ds, dd, tmp;

  fa = f(a, context);
  if (fabs(fa) <= BF_EPS) {
    *t = a;
    return true;
  }

  fb = f(b, context);
  if (fabs(fb) <= BF_EPS) {
    *t = b;
    return true;
  }

  if (sgn(fa) == sgn(fb)) {
    return false;
  }

  c = a;
  fc = fa;
  for (;;) {
    if (fabs(fc) < fabs(fb)) {
      tmp = b; b = c; c = tmp;
      tmp = fb; fb = fc; fc = tmp;
      a = c;
      fa = fc;
    }
    if (fabs(b - c) <= BF_EPS) {
      break;
    }
    dm = (c - b)/2;
    df = fa - fb;
    ds = df == 0 ? dm : -fb*(a - b)/df;
    dd = sgn(ds) != sgn(dm) || fabs(ds) > fabs(dm) ? dm : ds;
    if (fabs(dd) < BF_EPS) {
      dd = BF_EPS*sgn(dm)/2;
    }
    d = b + dd;
    fd = f(d, context);
    if (fd == 0) {
      c = d;
      b = c;
      fc = fd;
      fb = fc;
      break;
    }
    a = b;
    b = d;
    fa = fb;
    fb = fd;
    if (sgn(fb) == sgn(fc)) {
      c = a;
      fc = fa;
    }
  }
  *t = (b + c)/2;
  return true;
}

int bfSignum(BfReal x) {
  if (x > 0) return 1;
  else if (x == 0) return 0;
  else return -1;
}
