#include <bf/lu.h>

/** Interface: Lu */

BfMat *bfLuSolve(BfLu const *lu, BfMat const *B) {
  return lu->vtbl->Solve(lu, B);
}

BfMat *bfLuSolveLower(BfLu const *lu, BfMat const *B, bool permute) {
  return lu->vtbl->SolveLower(lu, B, permute);
}

BfMat *bfLuSolveUpper(BfLu const *lu, BfMat const *B, bool permute) {
  return lu->vtbl->SolveUpper(lu, B, permute);
}

BfMat *bfLuScale(BfLu const *lu, BfMat const *B) {
  return lu->vtbl->Scale(lu, B);
}

BfVec *bfLuSolveVec(BfLu const *lu, BfVec const *b) {
  return lu->vtbl->SolveVec(lu, b);
}

BfVec *bfLuSolveLowerVec(BfLu const *lu, BfVec const *b, bool permute) {
  return lu->vtbl->SolveLowerVec(lu, b, permute);
}

BfVec *bfLuSolveUpperVec(BfLu const *lu, BfVec const *b, bool permute) {
  return lu->vtbl->SolveUpperVec(lu, b, permute);
}

BfVec *bfLuScaleVec(BfLu const *lu, BfVec const *b) {
  return lu->vtbl->ScaleVec(lu, b);
}

BfMat *bfLuGetMatView(BfLu *lu) {
  return lu->vtbl->GetMatView(lu);
}

/** Implementation: Lu */

void bfLuInit(BfLu *lu, BfLuVtable *vtbl) {
  lu->vtbl = vtbl;
}

void bfLuDeinit(BfLu *lu) {}
