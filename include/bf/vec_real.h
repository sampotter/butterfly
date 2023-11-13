#pragma once

#include "vec.h"

/** Interface: Vec */

BfVecReal *bfVecRealCopy(BfVecReal const *vec);
void bfVecRealDelete(BfVecReal **mat);
BfType bfVecRealGetType(BfVec const *vec);
BfReal *bfVecRealGetEltPtr(BfVecReal *vecReal, BfSize i);
BfVecReal *bfVecRealGetSubvecCopy(BfVecReal const *vec, BfSize i0, BfSize i1);
BfVecReal *bfVecRealGetSubvecView(BfVecReal *vecReal, BfSize i0, BfSize i1);
BfVecReal const *bfVecRealGetSubvecViewConst(BfVecReal const *vecReal, BfSize i0, BfSize i1);
void bfVecRealSetRange(BfVecReal *vecReal, BfSize i0, BfSize i1, BfVec const *other);
void bfVecRealSetMask(BfVecReal *vecReal, bool const *mask, BfVec const *otherVec);
void bfVecRealPrint(BfVec const *vec, FILE *fp);
BfReal bfVecRealDistMax(BfVecReal const *vecReal, BfVec const *otherVec);
BfReal bfVecRealNormMax(BfVecReal const *vecReal);
void bfVecRealAddInplace(BfVecReal *vecReal, BfVec const *otherVec);
void bfVecRealRecipInplace(BfVec *vec);
void bfVecRealPermute(BfVecReal *vecReal, BfPerm const *perm);
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
BfVecReal *bfVecRealNewFromRealArray(BfRealArray *realArray, BfPolicy policy);
BfVecReal *bfVecRealNewFromCsv(char const *path);
void bfVecRealInit(BfVecReal *vecReal, BfSize size);
void bfVecRealInitView(BfVecReal *vecReal, BfSize size, BfSize stride, BfReal *data);
void bfVecRealInitFromPtr(BfVecReal *vecReal, BfSize size, BfSize stride, BfReal *data, BfPolicy policy);
void bfVecRealInitFromRealArray(BfVecReal *vecReal, BfRealArray *realArray, BfPolicy policy);
void bfVecRealInitFromCsv(BfVecReal *vecReal, char const *path);
void bfVecRealDeinit(BfVecReal *vecReal);
void bfVecRealDealloc(BfVecReal **vecReal);
void bfVecRealDeinitAndDealloc(BfVecReal **vecReal);
BfReal bfVecRealGetElt(BfVecReal const *vecReal, BfSize i);
void bfVecRealGetValues(BfVecReal const *vecReal, BfSize n, BfSize const *inds, BfReal *values);
BfRealArray *bfVecRealGetArrayView(BfVecReal *vecReal);
void bfVecRealDump(BfVecReal const *vecReal, char const *path);
BfPerm *bfVecRealArgsort(BfVecReal const *vecReal);
BfSize bfVecRealGetSize(BfVecReal const *vecReal);
