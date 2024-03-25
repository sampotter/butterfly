#include <bf/chol.h>

/** Interface: Chol */

void bfCholDelete(BfChol **chol) {
  (*chol)->vtbl->Delete(chol);
}

BfMat *bfCholSolve(BfChol const *chol, BfMat const *B) {
  return chol->vtbl->Solve(chol, B);
}

BfVec *bfCholSolveVec(BfChol const *chol, BfVec const *b) {
  return chol->vtbl->SolveVec(chol, b);
}

/** Implementation: Chol */

void bfCholInit(BfChol *chol, BfCholVtable *vtbl) {
  chol->vtbl = vtbl;
}

void bfCholDeinit(BfChol *chol) {
  (void)chol;
}
