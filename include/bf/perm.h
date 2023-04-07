#pragma once

#include "def.h"
#include "types.h"

typedef struct BfPerm {
  BfSize *index;
  BfSize size;
} BfPerm;

BfPerm *bfPermNew();
void bfPermInitEmpty(BfPerm *perm, BfSize size);
void bfPermDeinit(BfPerm *perm);
void bfPermDealloc(BfPerm **perm);
void bfPermDelete(BfPerm **perm);
BfPerm bfPermIdentity(BfSize size);
BfPerm bfPermGetReversePerm(BfPerm const *perm);
