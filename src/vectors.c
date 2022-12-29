#include <bf/vectors.h>

#include <math.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>

void bfVector3Scale(BfVector3 u, BfReal alpha) {
  u[0] *= alpha;
  u[1] *= alpha;
  u[2] *= alpha;
}

BfReal bfVector3Norm(BfVector3 const u) {
  return sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
}

BfReal bfVector3Dot(BfVector3 const u, BfVector3 const v) {
  return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

void bfVector3Cross(BfVector3 const u, BfVector3 const v, BfVector3 w) {
  w[0] = u[1]*v[2] - u[2]*v[1];
  w[1] = u[2]*v[0] - u[0]*v[2];
  w[2] = u[0]*v[1] - u[1]*v[0];
}

BfVectors2 const *bfVectors2ConstViewFromMat(BfMat const *mat) {
  return bfVectors2ConstViewFromMatDenseReal(bfMatConstToMatDenseRealConst(mat));
}

BfVectors2 const *bfVectors2ConstViewFromMatDenseReal(BfMatDenseReal const *matDenseReal) {
  BEGIN_ERROR_HANDLING();

  BfVectors2 *vectors = NULL;

  BfMat const *mat = bfMatDenseRealConstToMatConst(matDenseReal);
  BfMatDense const *matDense = bfMatDenseRealConstToMatDenseConst(matDenseReal);

  if (bfMatDenseRealGetNumCols(mat) != 2)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (bfMatDenseGetColStride(matDense) != 1)
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

void bfReadVectors2FromFile(char const *path, BfVectors2 *vectors) {
  BEGIN_ERROR_HANDLING();

  /* open the file for reading */
  FILE *fp = fopen(path, "r");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* get the size of the file */
  fseek(fp, 0, SEEK_END);
  BfSize size = ftell(fp);
  fseek(fp, 0, SEEK_SET);

  /* make sure the binary file is the right size */
  if (size % sizeof(BfPoint2) != 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* get the number of vectors */
  vectors->size = size/sizeof(BfPoint2);

  /* allocate space for the vectors */
  vectors->data = malloc(size);
  if (vectors->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* read them in */
  fread(vectors->data, sizeof(BfPoint2), vectors->size, fp);
  // TODO: error-handling

  END_ERROR_HANDLING() {
    free(vectors->data);
  }

  fclose(fp);
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

void bfSaveVectors2(BfVectors2 const *vectors, char const *path) {
  BEGIN_ERROR_HANDLING();

  FILE *fp = fopen(path, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  fwrite(vectors->data, vectors->size, sizeof(BfPoint2), fp);
  // TODO: error-handling

  END_ERROR_HANDLING() {}

  fclose(fp);
}
