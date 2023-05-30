#pragma once

#include "types.h"

/** Interface: */

void bfLuDelete(BfLu **lu);
BfMat *bfLuSolve(BfLu const *lu, BfMat const *B);
BfMat *bfLuSolveLower(BfLu const *lu, BfMat const *B, bool permute);
BfMat *bfLuSolveUpper(BfLu const *lu, BfMat const *B, bool permute);
BfMat *bfLuScale(BfLu const *lu, BfMat const *B);
BfVec *bfLuSolveVec(BfLu const *lu, BfVec const *b);
BfVec *bfLuSolveLowerVec(BfLu const *lu, BfVec const *b, bool permute);
BfVec *bfLuSolveUpperVec(BfLu const *lu, BfVec const *b, bool permute);
BfVec *bfLuScaleVec(BfLu const *lu, BfVec const *b);
BfMat *bfLuGetMatView(BfLu *lu);

typedef struct BfLuVtable {
  __typeof__(&bfLuDelete) Delete;
  __typeof__(&bfLuSolve) Solve;
  __typeof__(&bfLuSolveLower) SolveLower;
  __typeof__(&bfLuSolveUpper) SolveUpper;
  __typeof__(&bfLuScale) Scale;
  __typeof__(&bfLuSolveVec) SolveVec;
  __typeof__(&bfLuSolveLowerVec) SolveLowerVec;
  __typeof__(&bfLuSolveUpperVec) SolveUpperVec;
  __typeof__(&bfLuScaleVec) ScaleVec;
  __typeof__(&bfLuGetMatView) GetMatView;
} BfLuVtable;

/** Implementation: Lu */

struct BfLu {
  BfLuVtable *vtbl;
};

void bfLuInit(BfLu *lu, BfLuVtable *vtbl);
void bfLuDeinit(BfLu *lu);
