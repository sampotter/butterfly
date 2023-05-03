#pragma once

#include "vec.h"

/** Interface: Vec */

BfVec *bfVecRealCopy(BfVecReal const *vec);
void bfVecRealDelete(BfVec **mat);
BfType bfVecRealGetType(BfVec const *vec);
BfPtr bfVecRealGetEltPtr(BfVec *vec, BfSize i);
BfVec *bfVecRealGetSubvecCopy(BfVec const *vec, BfSize i0, BfSize i1);
BfVecReal *bfVecRealGetSubvecView(BfVecReal *vecReal, BfSize i0, BfSize i1);
BfVecReal const *bfVecRealGetSubvecViewConst(BfVecReal const *vecReal, BfSize i0, BfSize i1);
void bfVecRealSetRange(BfVecReal *vecReal, BfSize i0, BfSize i1, BfVec const *other);
void bfVecRealPrint(BfVec const *vec, FILE *fp);
BfReal bfVecRealDistMax(BfVecReal const *vecReal, BfVec const *otherVec);
BfReal bfVecRealNormMax(BfVec const *vec);
void bfVecRealAddInplace(BfVecReal *vecReal, BfVec const *otherVec);
void bfVecRealRecipInplace(BfVec *vec);
void bfVecRealPermute(BfVec *vec, BfPerm const *perm);
BfVec *bfVecRealConcat(BfVec const *vec, BfVec const *otherVec);
void bfVecRealSave(BfVecReal const *vecReal, char const *path);
void bfVecRealDaxpy(BfVecReal *vecReal, BfReal scale, BfVecReal const *otherVecReal);
void bfVecRealDscal(BfVecReal *vecReal, BfReal scale);

struct BfVecReal {
  BfVec super;
  BfSize stride;
  BfReal *data;
};

BfVec *bfVecRealToVec(BfVecReal *vecReal);
BfVec const *bfVecRealConstToVecConst(BfVecReal const *vecReal);

BfVecReal *bfVecToVecReal(BfVec *vec);
BfVecReal const *bfVecConstToVecRealConst(BfVec const *vec);

BfVecReal *bfVecRealNew();
BfVecReal *bfVecRealNewEmpty(BfSize n);
BfVecReal *bfVecRealNewWithValue(BfSize n, BfReal value);
BfVecReal *bfVecRealNewRandn(BfSize n);
BfVecReal *bfVecRealNewStdBasis(BfSize n, BfSize i);
BfVecReal *bfVecRealFromFile(char const *path, BfSize size);
void bfVecRealInit(BfVecReal *vecReal, BfSize size);
void bfVecRealInitFrom(BfVecReal *vecReal, BfSize size, BfSize stride, BfReal const *data);
void bfVecRealInitView(BfVecReal *vecReal, BfSize size, BfSize stride, BfReal *data);
void bfVecRealDeinit(BfVecReal *vecReal);
void bfVecRealDealloc(BfVecReal **vecReal);
void bfVecRealDeinitAndDealloc(BfVecReal **vecReal);
void bfVecRealDump(BfVecReal const *vecReal, char const *path);
