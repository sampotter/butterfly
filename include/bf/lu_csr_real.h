#pragma once

#include "def.h"
#include "lu.h"

/** Interface: Lu */

BfMat *bfLuCsrRealSolve(BfLuCsrReal const *luCsrReal, BfMat const *B);
BfMat *bfLuCsrRealSolveLower(BfLuCsrReal const *luCsrReal, BfMat const *B, bool permute);
BfMat *bfLuCsrRealSolveUpper(BfLuCsrReal const *luCsrReal, BfMat const *B, bool permute);
BfMat *bfLuCsrRealScale(BfLuCsrReal const *luCsrReal, BfMat const *B);
BfVec *bfLuCsrRealSolveVec(BfLuCsrReal const *luCsrReal, BfVec const *b);
BfVec *bfLuCsrRealSolveLowerVec(BfLuCsrReal const *luCsrReal, BfVec const *b, bool permute);
BfVec *bfLuCsrRealSolveUpperVec(BfLuCsrReal const *luCsrReal, BfVec const *b, bool permute);
BfVec *bfLuCsrRealScaleVec(BfLuCsrReal const *luCsrReal, BfVec const *b);

/** Upcasting: LuCsrReal -> Lu */

/** Downcasting: Lu -> LuCsrReal */

/** Implementation: LuCsrReal */

struct BfLuCsrReal {
  BfLu super;

  /** Private data used internally by UMFPACK: */

  int *rowptr;
  int *colind;
  BfReal *data;
  void *symbolic;
  void *numeric;
};

BfLuCsrReal *bfLuCsrRealNew(void);
void bfLuCsrRealInit(BfLuCsrReal *luCsrReal, BfMat const *mat);
void bfLuCsrRealDeinit(BfLuCsrReal *luCsrReal);
void bfLuCsrRealDealloc(BfLuCsrReal **luCsrReal);
void bfLuCsrRealDeinitAndDealloc(BfLuCsrReal **luCsrReal);
void bfLuCsrRealDump(BfLuCsrReal const *luCsrReal);
