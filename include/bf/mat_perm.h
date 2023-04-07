#pragma once

#include "mat.h"

/** Interface: Mat */

BfType bfMatPermGetType(BfMatPerm const *matPerm);
BfSize bfMatPermGetNumRows(BfMatPerm const *matPerm);
BfSize bfMatPermGetNumCols(BfMatPerm const *matPerm);
BfMat *bfMatPermSolve(BfMatPerm const *matPerm, BfMat const *mat);
BfMat *bfMatPermGetInverse(BfMatPerm const *matPerm);

/** Upcasting: MatPerm -> Mat */

BfMat *bfMatPermToMat(BfMatPerm *matPerm);

/** Downcasting: Mat -> MatPerm */

/** Implementation: MatPerm */

typedef struct BfMatPermImpl BfMatPermImpl;

struct BfMatPerm {
  BfMat super;
  BfMatPermImpl *impl;
};

BfMatPerm *bfMatPermNew();
BfMatPerm *bfMatPermNewFromPerm(BfPerm const *perm);
BfMatPerm *bfMatPermNewViewFromLapackPivots(BfSize size, BfLapackInt *ipiv);
void bfMatPermInitFromPerm(BfMatPerm *matPerm, BfPerm const *perm);
void bfMatPermInitViewFromLapackPivots(BfMatPerm *matPerm, BfSize size, BfLapackInt *ipiv);
void bfMatPermDeinit(BfMatPerm *matPerm);
void bfMatPermDealloc(BfMatPerm **matPerm);
void bfMatPermDelete(BfMatPerm **matPerm);
