#pragma once

#include "lu.h"

/** Interface: Lu */

void bfLuDenseComplexDelete(BfLuDenseComplex **luDenseComplex);
BfMat *bfLuDenseComplexSolve(BfLuDenseComplex const *luDenseComplex, BfMat const *B);
BfMat *bfLuDenseComplexSolveLower(BfLuDenseComplex const *luDenseComplex, BfMat const *B, bool permute);
BfMat *bfLuDenseComplexSolveUpper(BfLuDenseComplex const *luDenseComplex, BfMat const *B, bool permute);
BfMat *bfLuDenseComplexScale(BfLuDenseComplex const *luDenseComplex, BfMat const *B);
BfVec *bfLuDenseComplexSolveVec(BfLuDenseComplex const *luDenseComplex, BfVec const *b);
BfVec *bfLuDenseComplexSolveLowerVec(BfLuDenseComplex const *luDenseComplex, BfVec const *b, bool permute);
BfVec *bfLuDenseComplexSolveUpperVec(BfLuDenseComplex const *luDenseComplex, BfVec const *b, bool permute);
BfVec *bfLuDenseComplexScaleVec(BfLuDenseComplex const *luDenseComplex, BfVec const *b);
BfMat *bfLuDenseComplexGetMatView(BfLuDenseComplex *luDenseComplex);

/** Upcasting: LuDenseComplex -> Lu */

BfLu *bfLuDenseComplexToLu(BfLuDenseComplex *luDenseComplex);

/** Downcasting: Lu -> LuDenseComplex */

/** Implementation: LuDenseComplex */

typedef struct BfLuDenseComplexImpl BfLuDenseComplexImpl;

struct BfLuDenseComplex {
  BfLu super;
  BfLuDenseComplexImpl *impl;
};

BfLuDenseComplex *bfLuDenseComplexNew();
void bfLuDenseComplexInit(BfLuDenseComplex *luDenseComplex, BfMat const *mat);
void bfLuDenseComplexDeinit(BfLuDenseComplex *luDenseComplex);
void bfLuDenseComplexDealloc(BfLuDenseComplex **luDenseComplex);
void bfLuDenseComplexDeinitAndDealloc(BfLuDenseComplex **luDenseComplex);
void bfLuDenseComplexDump(BfLuDenseComplex const *luDenseComplex);
