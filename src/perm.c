#include <bf/perm.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_perm.h>
#include <bf/mem.h>

typedef enum {
  PERM_TYPE_BF,
  PERM_TYPE_LAPACK,
} PermType;

BfPerm *bfPermNew() {
  BF_ERROR_BEGIN();

  BfPerm *perm = bfMemAlloc(1, sizeof(BfPerm));
  if (perm == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BF_ERROR_END() {
    perm = NULL;
  }

  return perm;
}

void bfPermInitEmpty(BfPerm *perm, BfSize size) {
  BF_ERROR_BEGIN();

  perm->index = bfMemAlloc(size, sizeof(BfSize));
  if (perm->index == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0; i < size; ++i)
    perm->index[i] = BF_SIZE_BAD_VALUE;

  perm->size = size;

  BF_ERROR_END() {
    bfPermDeinit(perm);
  }
}

void bfPermDeinit(BfPerm *perm) {
  bfMemFree(perm->index);
  perm->index = NULL;
  perm->size = BF_SIZE_BAD_VALUE;
}

void bfPermDealloc(BfPerm **perm) {
  bfMemFree(*perm);
  *perm = NULL;
}

void bfPermDelete(BfPerm **perm) {
  bfPermDeinit(*perm);
  bfPermDealloc(perm);
}

BfPerm *bfPermCopy(BfPerm const *perm) {
  BF_ERROR_BEGIN();

  BfPerm *permCopy = bfPermNew();
  HANDLE_ERROR();

  bfPermInitEmpty(permCopy, perm->size);
  HANDLE_ERROR();

  permCopy->index = bfMemAlloc(perm->size, sizeof(BfSize));
  if (permCopy->index == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  bfMemCopy(perm->index, perm->size, sizeof(BfSize), permCopy->index);

  BF_ERROR_END() {
    BF_ASSERT(false);
  }

  return permCopy;
}

BfPerm bfPermIdentity(BfSize size) {
  BF_ERROR_BEGIN();

  BfPerm perm = {.index = bfMemAlloc(size, sizeof(BfSize)), .size = size};

  if (perm.index == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0; i < size; ++i)
    perm.index[i] = i;

  BF_ERROR_END() {
    bfMemFree(perm.index);
    perm.index = NULL;
    perm.size = BF_SIZE_BAD_VALUE;
  }

  return perm;
}

BfPerm bfPermGetReversePerm(BfPerm const *perm) {
  BF_ERROR_BEGIN();

  BfPerm revPerm = {
    .index = bfMemAlloc(perm->size, sizeof(BfSize)),
    .size = perm->size
  };

  if (revPerm.index == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0; i < revPerm.size; ++i)
    revPerm.index[perm->index[i]] = i;

  BF_ERROR_END() {
    bfMemFree(revPerm.index);
    revPerm.index = NULL;
    revPerm.size = BF_SIZE_BAD_VALUE;
  }

  return revPerm;
}

BfSize bfPermGetSize(BfPerm const *perm) {
  return perm->size;
}
