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

BfPerm *bfPermGetView(BfPerm *perm) {
  BF_ERROR_BEGIN();

  BfPerm *permView = bfPermNew();
  HANDLE_ERROR();

  permView->index = perm->index;
  permView->size = perm->size;
  permView->isView = true;

  BF_ERROR_END() {
    BF_DIE();
  }

  return permView;
}

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

  perm->isView = false;

  BF_ERROR_END() {
    bfPermDeinit(perm);
  }
}

void bfPermDeinit(BfPerm *perm) {
  if (!perm->isView)
    bfMemFree(perm->index);

  perm->index = NULL;
  perm->size = BF_SIZE_BAD_VALUE;
}

void bfPermDealloc(BfPerm **perm) {
  bfMemFree(*perm);
  *perm = NULL;
}

void bfPermDeinitAndDealloc(BfPerm **perm) {
  bfPermDeinit(*perm);
  bfPermDealloc(perm);
}

BfPerm *bfPermCopy(BfPerm const *perm) {
  BF_ERROR_BEGIN();

  BfPerm *permCopy = bfPermNew();
  HANDLE_ERROR();

  permCopy->size = perm->size;

  permCopy->index = bfMemAlloc(perm->size, sizeof(BfSize));
  if (permCopy->index == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  bfMemCopy(perm->index, perm->size, sizeof(BfSize), permCopy->index);

  permCopy->isView = false;

  BF_ERROR_END() {
    BF_DIE();
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

  perm.isView = false;

  BF_ERROR_END() {
    BF_DIE();
  }

  return perm;
}

BfPerm *bfPermGetReversePerm(BfPerm const *perm) {
  BF_ERROR_BEGIN();

  BfPerm *revPerm = bfPermNew();
  HANDLE_ERROR();

  bfPermInitEmpty(revPerm, perm->size);
  HANDLE_ERROR();

  for (BfSize i = 0; i < revPerm->size; ++i)
    revPerm->index[perm->index[i]] = i;

  BF_ERROR_END() {
    BF_DIE();
  }

  return revPerm;
}

BfSize bfPermGetSize(BfPerm const *perm) {
  return perm->size;
}
