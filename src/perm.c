#include <bf/perm.h>

#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>

void bfPermDeinit(BfPerm *perm) {
  free(perm->index);
  perm->index = NULL;

  perm->size = BF_SIZE_BAD_VALUE;
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
