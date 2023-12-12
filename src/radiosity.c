#include <bf/radiosity.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_csr_real.h>
#include <bf/real_array.h>
#include <bf/size_array.h>
#include <bf/trimesh.h>

static BfReal getRadiosityKernelValue(BfTrimesh const *trimesh, BfSize srcInd, BfSize tgtInd) {
  (void)trimesh;
  (void)srcInd;
  (void)tgtInd;
  BF_DIE();
}

BfMatCsrReal *bfRadiosityGetViewFactorMatrix(BfTrimesh const *trimesh, BfSizeArray const *rowInds, BfSizeArray const *colInds) {
  BF_ERROR_BEGIN();

  BfSize numRows = bfSizeArrayGetSize(rowInds);
  BfSize numCols = bfSizeArrayGetSize(colInds);

  BfSizeArray *rowptr = bfSizeArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  bfSizeArrayAppend(rowptr, 0);

  BfSizeArray *colind = bfSizeArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  BfRealArray *data = bfRealArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  for (BfSize i = 0; i < bfSizeArrayGetSize(rowInds); ++i) {
    BfSize rowInd = bfSizeArrayGet(rowInds, i);

    BfSizeArray *visibleColInds = bfTrimeshGetVisibility(trimesh, i, colInds);
    HANDLE_ERROR();

    bfSizeArrayExtend(colind, visibleColInds);
    HANDLE_ERROR();

    BfSize numVisibleColInds = bfSizeArrayGetSize(visibleColInds);

    bfSizeArrayAppend(rowptr, numVisibleColInds + bfSizeArrayGetLast(rowptr));
    HANDLE_ERROR();

    for (BfSize j = 0; j < numVisibleColInds; ++j) {
      BfSize colInd = bfSizeArrayGet(visibleColInds, j);
      BfReal value = getRadiosityKernelValue(trimesh, rowInd, colInd);
      bfRealArrayAppend(data, value);
    }

    bfSizeArrayDeinitAndDealloc(&visibleColInds);
  }

  BfMatCsrReal *matCsrReal = bfMatCsrRealNewFromArrays(numRows, numCols, rowptr, colind, data, BF_POLICY_STEAL);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  bfSizeArrayDeinitAndDealloc(&rowptr);
  bfSizeArrayDeinitAndDealloc(&colind);
  bfRealArrayDeinitAndDealloc(&data);

  return matCsrReal;
}
