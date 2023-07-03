#include <bf/vectors.h>

#include <math.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>

void bfVector2Normalize(BfVector2 u) {
  BfReal mag = hypot(u[0], u[1]);
  u[0] /= mag;
  u[1] /= mag;
}

void bfVector2Reject(BfVector2 u, BfVector2 const v) {
  BfReal v_dot_u = v[0]*u[0] + v[1]*u[1];
  u[0] -= v_dot_u*v[0];
  u[1] -= v_dot_u*v[1];
}

void bfVector2Negate(BfVector2 u) {
  u[0] *= -1;
  u[1] *= -1;
}

void bfVector2Rotate(BfVector2 u, BfReal theta) {
  BfReal c = cos(theta);
  BfReal s = sin(theta);
  BfReal tmp = c*u[0] - s*u[1];
  u[1] = s*u[0] + c*u[1];
  u[0] = tmp;
}

void bfVector3Copy(BfVector3 u, BfVector3 const v) {
  u[0] = v[0];
  u[1] = v[1];
  u[2] = v[2];
}

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

BfVectors2 *bfVectors2NewEmpty() {
  BF_ERROR_BEGIN();

  BfVectors2 *vectors = bfMemAlloc(1, sizeof(BfVectors2));
  HANDLE_ERROR();

  bfVectors2InitEmpty(vectors);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return vectors;
}

BfVectors2 *bfVectors2NewWithCapacity(BfSize capacity) {
  BF_ERROR_BEGIN();

  BfVectors2 *vectors = bfMemAlloc(1, sizeof(BfVectors2));
  HANDLE_ERROR();

  bfVectors2InitWithCapacity(vectors, capacity);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return vectors;
}

BfVectors2 *bfVectors2NewFromFile(char const *path) {
  BF_ERROR_BEGIN();

  BfVectors2 *vectors = bfMemAlloc(1, sizeof(BfVectors2));
  HANDLE_ERROR();

  bfVectors2InitFromFile(path, vectors);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return vectors;
}

BfVectors2 const *bfVectors2ConstViewFromMat(BfMat const *mat) {
  return bfVectors2ConstViewFromMatDenseReal(bfMatConstToMatDenseRealConst(mat));
}

BfVectors2 const *bfVectors2ConstViewFromMatDenseReal(BfMatDenseReal const *matDenseReal) {
  BF_ERROR_BEGIN();

  BfVectors2 *vectors = NULL;

  BfMat const *mat = bfMatDenseRealConstToMatConst(matDenseReal);
  BfMatDense const *matDense = bfMatDenseRealConstToMatDenseConst(matDenseReal);

  if (bfMatDenseRealGetNumCols(mat) != 2)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (bfMatDenseGetColStride(matDense) != 1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  vectors = bfMemAlloc(1, sizeof(BfVectors2));
  vectors->data = (BfPoint2 *)matDenseReal->data;
  vectors->size = bfMatDenseRealGetNumRows(mat);
  vectors->capacity = BF_SIZE_BAD_VALUE;
  vectors->isView = true;

  BF_ERROR_END() {
    bfMemFree(vectors);
    vectors = NULL;
  }

  return vectors;
}

BfVectors2 bfGetUninitializedVectors2() {
  return (BfVectors2) {.data = NULL, .size = 0};
}

void bfVectors2Init(BfVectors2 *vectors) {
  BF_ERROR_BEGIN();

  vectors->size = 0;
  vectors->capacity = BF_ARRAY_DEFAULT_CAPACITY;
  vectors->isView = false;

  vectors->data = bfMemAlloc(vectors->capacity, sizeof(BfVector2));
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfVectors2InitEmpty(BfVectors2 *vectors) {
  BF_ERROR_BEGIN();

  vectors->size = 0;
  vectors->capacity = BF_ARRAY_DEFAULT_CAPACITY;
  vectors->isView = false;

  vectors->data = bfMemAlloc(vectors->capacity, sizeof(BfVector2));
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfVectors2InitWithCapacity(BfVectors2 *vectors, BfSize capacity) {
  BF_ERROR_BEGIN();

  vectors->size = 0;
  vectors->capacity = capacity;
  vectors->isView = false;

  vectors->data = bfMemAlloc(vectors->capacity, sizeof(BfVector2));
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfVectors2InitFromFile(char const *path, BfVectors2 *vectors) {
  BF_ERROR_BEGIN();

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
  vectors->data = bfMemAlloc(size, sizeof(char));
  if (vectors->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* read them in */
  fread(vectors->data, sizeof(BfPoint2), vectors->size, fp);
  // TODO: error-handling

  BF_ERROR_END() {
    bfMemFree(vectors->data);
  }

  fclose(fp);
}

void bfFreeVectors2(BfVectors2 *vectors) {
  if (!vectors->isView)
    bfMemFree(vectors->data);
}

void bfVectors2Dealloc(BfVectors2 **vectors) {
  bfMemFree(*vectors);
  *vectors = NULL;
}

void bfVectors2DeinitAndDealloc(BfVectors2 **vectors) {
  bfFreeVectors2(*vectors);
  bfVectors2Dealloc(vectors);
}

void bfGetVectorsByIndex(BfVectors2 const *vectors,
                        BfSize numInds, BfSize const *inds,
                        BfVectors2 *indexedVectors)
{
  BF_ERROR_BEGIN();

  bfVectors2InitEmpty(indexedVectors);
  HANDLE_ERROR();

  BfVector2 const *vector = (BfVector2 const *)vectors->data;

  for (BfSize i = 0; i < numInds; ++i) {
    bfVectors2Append(indexedVectors, vector[inds[i]]);
    HANDLE_ERROR();
  }

  BF_ERROR_END() {
    bfFreeVectors2(indexedVectors);
  }
}

void bfSaveVectors2(BfVectors2 const *vectors, char const *path) {
  BF_ERROR_BEGIN();

  FILE *fp = fopen(path, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  fwrite(vectors->data, vectors->size, sizeof(BfPoint2), fp);
  // TODO: error-handling

  BF_ERROR_END() {}

  fclose(fp);
}

void bfVectors2Append(BfVectors2 *vectors, BfPoint2 const p) {
  BF_ERROR_BEGIN();

  /* Grow the array if we're at capacity */
  if (vectors->size == vectors->capacity) {
    vectors->capacity *= 2;

    BfPoint2 *newData = bfMemRealloc(vectors->data, vectors->capacity, sizeof(BfPoint2));
    if (newData == NULL)
      RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

    vectors->data = newData;
  }

  /* Append new point */
  bfMemCopy(p, 1, sizeof(BfPoint2), &vectors->data[vectors->size++]);

  BF_ERROR_END() {}
}

void bfVectors2Extend(BfVectors2 *vectors, BfVectors2 const *newVectors) {
  /* TODO: bad implementation to start... fix */
  for (BfSize i = 0; i < newVectors->size; ++i) {
    BfPoint2 v;
    bfVectors2Get(newVectors, i, v);
    bfVectors2Append(vectors, v);
  }
}

void bfVectors2Get(BfVectors2 const *vectors, BfSize i, BfVector2 v) {
  BF_ERROR_BEGIN();

  if (i >= vectors->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  bfMemCopy(&vectors->data[i], 1, sizeof(BfVector2), v);

  BF_ERROR_END() {
    BF_DIE();
  }
}

BfSize bfVectors2GetSize(BfVectors2 const *vectors) {
  return vectors->size;
}

void bfVectors2Set(BfVectors2 *vectors, BfSize i, BfVector2 const v) {
  BF_ERROR_BEGIN();

  if (i >= vectors->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  bfMemCopy(v, 1, sizeof(BfVector2), &vectors->data[i]);

  BF_ERROR_END() {
    BF_DIE();
  }
}

BfVectors2 *bfVectors2GetRangeView(BfVectors2 *vectors, BfSize i0, BfSize i1) {
  BF_ERROR_BEGIN();

  if ((i0 >= vectors->size && i1 > vectors->size) || i0 > vectors->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (i0 > i1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfVectors2 *vectorsView = bfMemAlloc(1, sizeof(BfVectors2));
  if (vectorsView == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  vectorsView->size = i1 - i0;
  vectorsView->capacity = BF_SIZE_BAD_VALUE;
  vectorsView->isView = true;
  vectorsView->data = &vectors->data[i0];

  BF_ERROR_END() {
    BF_DIE();
  }

  return vectorsView;
}

void bfVectors3InitEmpty(BfVectors3 *vectors, BfSize numVectors) {
  if (numVectors == 0) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  vectors->size = numVectors;

  vectors->data = bfMemAlloc(numVectors, sizeof(BfVector3));
  if (vectors->data == NULL)
    bfSetError(BF_ERROR_MEMORY_ERROR);
}

void bfVectors3Deinit(BfVectors3 *vectors) {
  (void)vectors;
  BF_DIE();
}

void bfVectors3GetByIndex(BfVectors3 const *vectors, BfSize numInds, BfSize const *inds, BfVectors3 *indexedVectors) {
  BF_ERROR_BEGIN();

  bfVectors3InitEmpty(indexedVectors, numInds);
  HANDLE_ERROR();

  for (BfSize i = 0; i < numInds; ++i)
    bfVector3Copy(indexedVectors->data[i], vectors->data[inds[i]]);

  BF_ERROR_END() {
    bfVectors3Deinit(indexedVectors);
  }
}
