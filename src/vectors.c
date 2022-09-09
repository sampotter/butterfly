#include <bf/vectors.h>

#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>

BfVectors2 const *bfVectors2ConstViewFromMat(BfMat const *mat) {
  return bfVectors2ConstViewFromMatDenseReal(bfMatConstToMatDenseRealConst(mat));
}

BfVectors2 const *bfVectors2ConstViewFromMatDenseReal(BfMatDenseReal const *matDenseReal) {
  BEGIN_ERROR_HANDLING();

  BfVectors2 *vectors = NULL;

  BfMat const *mat = bfMatDenseRealConstToMatConst(matDenseReal);

  if (bfMatDenseRealGetNumCols(mat) != 2)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (matDenseReal->colStride != 1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  vectors = malloc(sizeof(BfVectors2));
  vectors->data = (BfPoint2 *)matDenseReal->data;
  vectors->size = bfMatDenseRealGetNumRows(mat);

  END_ERROR_HANDLING() {
    free(vectors);
    vectors = NULL;
  }

  return vectors;
}

BfVectors2 bfGetUninitializedVectors2() {
  return (BfVectors2) {.data = NULL, .size = 0};
}

void bfInitEmptyVectors2(BfVectors2 *vectors, BfSize numVectors) {
  if (numVectors == 0) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  vectors->size = numVectors;

  vectors->data = malloc(numVectors*sizeof(BfVector2));
  if (vectors->data == NULL)
    bfSetError(BF_ERROR_MEMORY_ERROR);
}

void bfFreeVectors2(BfVectors2 *vectors) {
  free(vectors->data);
}

void bfGetVectorsByIndex(BfVectors2 const *vectors,
                        BfSize numInds, BfSize const *inds,
                        BfVectors2 *indexedVectors)
{
  BEGIN_ERROR_HANDLING();

  bfInitEmptyVectors2(indexedVectors, numInds);
  HANDLE_ERROR();

  BfVector2 const *vector = (BfVector2 const *)vectors->data;
  BfVector2 *indexedVector = indexedVectors->data;

  for (BfSize i = 0, j; i < numInds; ++i) {
    j = inds[i];
    indexedVector[i][0] = vector[j][0];
    indexedVector[i][1] = vector[j][1];
  }

  END_ERROR_HANDLING() {
    bfFreeVectors2(indexedVectors);
  }
}
