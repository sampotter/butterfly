#pragma once

#include "types.h"

/** Interface: */

void bfCholDelete(BfChol **chol);
BfMat *bfCholSolve(BfChol const *chol, BfMat const *B);
BfVec *bfCholSolveVec(BfChol const *chol, BfVec const *b);

typedef struct BfCholVtable {
  __typeof__(&bfCholDelete) Delete;
  __typeof__(&bfCholSolve) Solve;
  __typeof__(&bfCholSolveVec) SolveVec;
} BfCholVtable;

/** Implementation: Chol */

struct BfChol {
  BfCholVtable *vtbl;
};

void bfCholInit(BfChol *chol, BfCholVtable *vtbl);
void bfCholDeinit(BfChol *chol);
