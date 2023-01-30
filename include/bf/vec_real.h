#pragma once

#include "vec.h"

/** Interface: Vec */

BfVec *bfVecRealCopy(BfVec const *vec);
void bfVecRealDelete(BfVec **mat);
BfType bfVecRealGetType(BfVec const *vec);
BfPtr bfVecRealGetEltPtr(BfVec *vec, BfSize i);
BfVec *bfVecRealGetSubvecCopy(BfVec const *vec, BfSize i0, BfSize i1);
void bfVecRealPrint(BfVec const *vec, FILE *fp);
BfReal bfVecRealNormMax(BfVec const *vec);
void bfVecRealRecipInplace(BfVec *vec);
void bfVecRealPermute(BfVec *vec, BfPerm const *perm);
BfVec *bfVecRealConcat(BfVec const *vec, BfVec const *otherVec);

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
BfVecReal *bfVecRealFromFile(char const *path, BfSize size);
void bfVecRealInit(BfVecReal *vecReal, BfSize size);
void bfVecRealInitFrom(BfVecReal *vecReal, BfSize size, BfSize stride, BfReal const *data);
void bfVecRealInitView(BfVecReal *vecReal, BfSize size, BfSize stride, BfReal *data);
void bfVecRealDeinit(BfVecReal *vecReal);
void bfVecRealDealloc(BfVecReal **vecReal);
void bfVecRealDeinitAndDealloc(BfVecReal **vecReal);
void bfVecRealDump(BfVecReal const *vecReal, char const *path);
