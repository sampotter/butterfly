#pragma once

#include "mat.h"

typedef BfMat *(*MatMulFunc)(BfMat const *, void *);

/** Interface: Mat */

BfSize bfMatFuncGetNumRows(BfMatFunc const *matFunc);
BfSize bfMatFuncGetNumCols(BfMatFunc const *matFunc);
BfMat *bfMatFuncMul(BfMatFunc const *matFunc, BfMat const *otherMat);

/** Upcasting: MatFunc -> Mat */

BfMat *bfMatFuncToMat(BfMatFunc *matFunc);
BfMat const *bfMatFuncConstToMatConst(BfMatFunc const *matFunc);

/** Implementation: MatFunc */

struct BfMatFunc {
  BfMat super;
  MatMulFunc matMul;
  void *aux;
};

BfMatFunc *bfMatFuncNew();
void bfMatFuncInit(BfMatFunc *matFunc, BfSize numRows, BfSize numCols, MatMulFunc matMul, void *);
void bfMatFuncDeinit(BfMatFunc *matFunc);
void bfMatFuncDealloc(BfMatFunc **matFunc);
void bfMatFuncDeinitAndDealloc(BfMatFunc **matFunc);
