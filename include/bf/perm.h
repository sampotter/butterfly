#pragma once

#include <bf/def.h>

typedef struct BfPerm {
  BfSize *index;
  BfSize size;
} BfPerm;

void bfPermDeinit(BfPerm *perm);
BfPerm bfPermIdentity(BfSize size);
BfPerm bfPermGetReversePerm(BfPerm const *perm);
