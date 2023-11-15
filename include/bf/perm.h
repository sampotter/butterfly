#pragma once

#include "def.h"
#include "types.h"

typedef struct BfPerm {
  BfSize *index;
  BfSize size;
  bool isView;
} BfPerm;

BfPerm *bfPermGetView(BfPerm *perm);
BfPerm *bfPermGetRangeView(BfPerm *perm, BfSize i0, BfSize i1);
BfPerm *bfPermNew(void);
BfPerm *bfPermNewEmpty(BfSize size);
BfPerm *bfPermNewIdentity(BfSize size);
void bfPermInitEmpty(BfPerm *perm, BfSize size);
void bfPermInitIdentity(BfPerm *perm, BfSize size);
void bfPermDeinit(BfPerm *perm);
void bfPermDealloc(BfPerm **perm);
void bfPermDeinitAndDealloc(BfPerm **perm);
BfPerm *bfPermCopy(BfPerm const *perm);
BfPerm bfPermIdentity(BfSize size);
BfPerm *bfPermGetReversePerm(BfPerm const *perm);
BfSize bfPermGetSize(BfPerm const *perm);
BfSize bfPermGetNumBytes(BfPerm const *perm);
BfSize bfPermGetIndex(BfPerm const *perm, BfSize i);
void bfPermReverse(BfPerm *perm);
