#pragma once

#include "def.h"
#include "types.h"

typedef struct BfPerm {
  BfSize *index;
  BfSize size;
  bool isView;
} BfPerm;

BfPerm *bfPermGetView(BfPerm *perm);
BfPerm *bfPermNew(void);
void bfPermInitEmpty(BfPerm *perm, BfSize size);
void bfPermDeinit(BfPerm *perm);
void bfPermDealloc(BfPerm **perm);
void bfPermDeinitAndDealloc(BfPerm **perm);
BfPerm *bfPermCopy(BfPerm const *perm);
BfPerm bfPermIdentity(BfSize size);
BfPerm *bfPermGetReversePerm(BfPerm const *perm);
BfSize bfPermGetSize(BfPerm const *perm);
BfSize bfPermGetNumBytes(BfPerm const *perm);
