#include <bf/perm.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_perm.h>
#include <bf/mem.h>

#include "macros.h"

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

BfPerm *bfPermGetRangeView(BfPerm *perm, BfSize i0, BfSize i1) {
  BF_ERROR_BEGIN();

  if (i0 > i1 || i1 > perm->size)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  BfPerm *permRangeView = bfMemAlloc(1, sizeof(BfPerm));
  HANDLE_ERROR();

  permRangeView->index = &perm->index[i0];
  permRangeView->size = i1 - i0;
  permRangeView->isView = true;

  BF_ERROR_END() {
    BF_DIE();
  }

  return permRangeView;
}

BfPerm *bfPermNew(void) {
  BF_ERROR_BEGIN();

  BfPerm *perm = bfMemAlloc(1, sizeof(BfPerm));
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return perm;
}

BfPerm *bfPermNewEmpty(BfSize size) {
  BF_ERROR_BEGIN();

  BfPerm *perm = bfPermNew();
  HANDLE_ERROR();

  bfPermInitEmpty(perm, size);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return perm;
}

BfPerm *bfPermNewIdentity(BfSize size) {
  BF_ERROR_BEGIN();

  BfPerm *perm = bfPermNew();
  HANDLE_ERROR();

  bfPermInitIdentity(perm, size);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return perm;
}

void bfPermInitEmpty(BfPerm *perm, BfSize size) {
  BF_ERROR_BEGIN();

  perm->index = bfMemAlloc(size, sizeof(BfSize));
  HANDLE_ERROR();

  for (BfSize i = 0; i < size; ++i)
    perm->index[i] = BF_SIZE_BAD_VALUE;

  perm->size = size;

  perm->isView = false;

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfPermInitIdentity(BfPerm *perm, BfSize size) {
  BF_ERROR_BEGIN();

  perm->index = bfMemAlloc(size, sizeof(BfSize));
  HANDLE_ERROR();

  for (BfSize i = 0; i < size; ++i)
    perm->index[i] = i;

  perm->size = size;

  perm->isView = false;

  BF_ERROR_END() {
    BF_DIE();
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

  for (BfSize i = 0; i < perm->size; ++i)
    if (perm->index[i] >= perm->size)
      RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfPerm *revPerm = bfPermNew();
  HANDLE_ERROR();

  bfPermInitEmpty(revPerm, perm->size);
  HANDLE_ERROR();

  for (BfSize i = 0; i < revPerm->size; ++i) {
    /* Make sure we don't set an entry twice... */
    if (BF_SIZE_OK(revPerm->index[perm->index[i]]))
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);
    revPerm->index[perm->index[i]] = i;
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return revPerm;
}

BfSize bfPermGetSize(BfPerm const *perm) {
  return perm->size;
}

BfSize bfPermGetIndex(BfPerm const *perm, BfSize i) {
  if (i >= perm->size) {
    bfSetError(BF_ERROR_OUT_OF_RANGE);
    return BF_SIZE_BAD_VALUE;
  } else {
    return perm->index[i];
  }
}

void bfPermReverse(BfPerm *perm) {
  BF_ERROR_BEGIN();

  BfSize *tmp = bfMemAllocCopy(perm->index, perm->size, sizeof(BfSize));
  HANDLE_ERROR();

  for (BfSize i = 0; i < perm->size; ++i)
    perm->index[tmp[i]] = i;

  BF_ERROR_END() {
    BF_DIE();
  }

  bfMemFree(tmp);
}
