#include <bf/perm.h>

#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>

BfPerm *bfPermNew() {
  BEGIN_ERROR_HANDLING();

  BfPerm *perm = malloc(sizeof(BfPerm));
  if (perm == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {
    perm = NULL;
  }

  return perm;
}

void bfPermInitEmpty(BfPerm *perm, BfSize size) {
  BEGIN_ERROR_HANDLING();

  perm->index = malloc(size*sizeof(BfSize));
  if (perm->index == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0; i < size; ++i)
    perm->index[i] = BF_SIZE_BAD_VALUE;

  perm->size = size;

  END_ERROR_HANDLING() {
    bfPermDeinit(perm);
  }
}

void bfPermDeinit(BfPerm *perm) {
  free(perm->index);
  perm->index = NULL;
  perm->size = BF_SIZE_BAD_VALUE;
}

void bfPermDealloc(BfPerm **perm) {
  free(*perm);
  *perm = NULL;
}

void bfPermDelete(BfPerm **perm) {
  bfPermDeinit(*perm);
  bfPermDealloc(perm);
}

BfPerm bfPermIdentity(BfSize size) {
  BEGIN_ERROR_HANDLING();

  BfPerm perm = {.index = malloc(size*sizeof(BfSize)), .size = size};

  if (perm.index == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0; i < size; ++i)
    perm.index[i] = i;

  END_ERROR_HANDLING() {
    free(perm.index);
    perm.index = NULL;
    perm.size = BF_SIZE_BAD_VALUE;
  }

  return perm;
}

BfPerm bfPermGetReversePerm(BfPerm const *perm) {
  BEGIN_ERROR_HANDLING();

  BfPerm revPerm = {
    .index = malloc(perm->size*sizeof(BfSize)),
    .size = perm->size
  };

  if (revPerm.index == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0; i < revPerm.size; ++i)
    revPerm.index[perm->index[i]] = i;

  END_ERROR_HANDLING() {
    free(revPerm.index);
    revPerm.index = NULL;
    revPerm.size = BF_SIZE_BAD_VALUE;
  }

  return revPerm;
}
