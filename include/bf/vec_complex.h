#pragma once

#include "vec.h"

/** Interface: Vec */

BfVec *bfVecComplexCopy(BfVec const *vec);
void bfVecComplexDelete(BfVec **mat);
BfType bfVecComplexGetType(BfVec const *vec);
BfVec *bfVecComplexGetSubvecView(BfVec *vec, BfSize i0, BfSize i1);
void bfVecComplexPrint(BfVec const *vec, FILE *fp);
BfReal bfVecComplexNormMax(BfVec const *vec);
void bfVecComplexAddInplace(BfVec *vec, BfVec const *otherVec);
void bfVecComplexMulInplace(BfVec *vec, BfMat const *mat);
void bfVecComplexSolveInplace(BfVec *vec, BfMat const *mat);
BfMat *bfVecComplexGetGivensRotation(BfVec const *vec, BfSize srcInd, BfSize elimInd);
BfVec *bfVecComplexConcat(BfVec const *vec, BfVec const *otherVec);
void bfVecComplexSave(BfVecComplex const *vecComplex, char const *path);

struct BfVecComplex {
  BfVec super;
  BfSize stride;
  BfComplex *data;
};

BfVec *bfVecComplexToVec(BfVecComplex *vecComplex);

BfVecComplex *bfVecToVecComplex(BfVec *vec);
BfVecComplex const *bfVecConstToVecComplexConst(BfVec const *vec);

BfVecComplex *bfVecComplexNew();
BfVecComplex *bfVecComplexFromFile(char const *path, BfSize size);
void bfVecComplexInit(BfVecComplex *vecComplex, BfSize size);
void bfVecComplexInitView(BfVecComplex *vecComplex, BfSize size, BfSize stride, BfComplex *data);
void bfVecComplexDeinit(BfVecComplex *vecComplex);
void bfVecComplexDealloc(BfVecComplex **vecComplex);
void bfVecComplexDeinitAndDealloc(BfVecComplex **vecComplex);
