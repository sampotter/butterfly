#pragma once

#include "mat.h"

/** Interface: Mat */

void bfMatIdentityDelete(BfMatIdentity **matIdentity);
BfType bfMatIdentityGetType(BfMatIdentity const *matIdentity);

/** Implementation: MatIdentity */

struct BfMatIdentity {
  BfMat super;
};

BfMat *bfMatIdentityToMat(BfMatIdentity *matIdentity);

BfMatIdentity *bfMatIdentityNew();
void bfMatIdentityInit(BfMatIdentity *matIdentity, BfSize n);
void bfMatIdentityDeinit(BfMatIdentity *matIdentity);
void bfMatIdentityDealloc(BfMatIdentity **matIdentity);
void bfMatIdentityDeinitAndDealloc(BfMatIdentity **matIdentity);
