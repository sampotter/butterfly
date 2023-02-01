#include <bf/mat_identity.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .Delete = (__typeof__(&bfMatDelete))bfMatIdentityDelete,
  .GetType = (__typeof__(&bfMatGetType))bfMatIdentityGetType,
};

void bfMatIdentityDelete(BfMatIdentity **matIdentity) {
  bfMatIdentityDeinitAndDealloc(matIdentity);
}

BfType bfMatIdentityGetType(BfMatIdentity const *matIdentity) {
  (void)matIdentity;
  return BF_TYPE_MAT_IDENTITY;
}

/** Upcasting: MatIdentity -> Mat */

BfMat *bfMatIdentityToMat(BfMatIdentity *matIdentity) {
  return &matIdentity->super;
}

/** Implementation: MatIdentity */

BfMatIdentity *bfMatIdentityNew() {
  BEGIN_ERROR_HANDLING();

  BfMatIdentity *matIdentity = malloc(sizeof(BfMatIdentity));
  if (matIdentity == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return matIdentity;
}

void bfMatIdentityInit(BfMatIdentity *matIdentity, BfSize n) {
  BEGIN_ERROR_HANDLING();

  bfMatInit(&matIdentity->super, &MAT_VTABLE, n, n);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    bfMatIdentityDeinit(matIdentity);
  }
}

void bfMatIdentityDeinit(BfMatIdentity *matIdentity) {
  (void)matIdentity;
  assert(false);
}

void bfMatIdentityDealloc(BfMatIdentity **matIdentity) {
  free(*matIdentity);
  *matIdentity = NULL;
}

void bfMatIdentityDeinitAndDealloc(BfMatIdentity **matIdentity) {
  bfMatIdentityDeinit(*matIdentity);
  bfMatIdentityDealloc(matIdentity);
}
