#pragma once

#include "mat.h"

/** Interface: Mat */

void bfMatIdentityDelete(BfMatIdentity **matIdentity);
BfType bfMatIdentityGetType(BfMatIdentity const *matIdentity);
BfSize bfMatIdentityGetNumRows(BfMatIdentity const *matIdentity);
BfSize bfMatIdentityGetNumCols(BfMatIdentity const *matIdentity);
BfMat *bfMatIdentityGetRowRangeCopy(BfMatIdentity const *matIdentity, BfSize i0, BfSize i1);

/** Upcasting: MatIdentity -> Mat */

BfMat *bfMatIdentityToMat(BfMatIdentity *matIdentity);
BfMat const *bfMatIdentityConstToMatConst(BfMatIdentity const *matIdentity);

/** Downcasting: Mat -> MatIdentity */

BfMatIdentity *bfMatToMatIdentity(BfMat *mat);

/** Implementation: MatIdentity */

struct BfMatIdentity {
  BfMat super;
};

BfMatIdentity *bfMatIdentityNew();
void bfMatIdentityInit(BfMatIdentity *matIdentity, BfSize n);
void bfMatIdentityDeinit(BfMatIdentity *matIdentity);
void bfMatIdentityDealloc(BfMatIdentity **matIdentity);
void bfMatIdentityDeinitAndDealloc(BfMatIdentity **matIdentity);
