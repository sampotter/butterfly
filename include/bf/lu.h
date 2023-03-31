#pragma once

#include "types.h"

BfLu *bfLuNew();
void bfLuInit(BfLu *lu, BfMat const *mat);
void bfLuDeinit(BfLu *lu);
void bfLuDealloc(BfLu **lu);
void bfLuDeinitAndDealloc(BfLu **lu);
BfVec *bfLuSolveVec(BfLu const *lu, BfVec const *b);
BfVec *bfLuSolveLowerVec(BfLu const *lu, BfVec const *b, bool permute);
BfVec *bfLuSolveUpperVec(BfLu const *lu, BfVec const *b, bool permute);
BfVec *bfLuScaleVec(BfLu const *lu, BfVec const *b);
void bfLuDump(BfLu const *lu);
