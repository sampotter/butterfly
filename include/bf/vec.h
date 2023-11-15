#pragma once

#include "def.h"
#include "dtype.h"
#include "perm.h"
#include "types.h"

#include <stdio.h>

typedef enum BfVecProps {
  BF_VEC_PROPS_NONE = 0,
  BF_VEC_PROPS_VIEW = (1 << 0)
} BfVecProps;

/** Interface: Vec */

BfVec *bfVecCopy(BfVec const *);
void bfVecDelete(BfVec **);
BfType bfVecGetType(BfVec const *);
bool bfVecInstanceOf(BfVec const *, BfType);
BfPtr bfVecGetEltPtr(BfVec *, BfSize);
BfVec *bfVecGetSubvecCopy(BfVec const *, BfSize, BfSize);
BfVec *bfVecGetSubvecView(BfVec *, BfSize, BfSize);
BfVec const *bfVecGetSubvecViewConst(BfVec const *, BfSize, BfSize);
void bfVecSetRange(BfVec *, BfSize, BfSize, BfVec const *);
void bfVecSetMask(BfVec *, bool const *, BfVec const *);
void bfVecPrint(BfVec const *, FILE *);
BfReal bfVecDist(BfVec const *, BfVec const *);
BfReal bfVecDistMax(BfVec const *, BfVec const *);
BfReal bfVecNormMax(BfVec const *);
void bfVecScaleByReal(BfVec *, BfReal);
void bfVecAddInplace(BfVec *, BfVec const *);
void bfVecMulInplace(BfVec *, BfMat const *);
void bfVecSolveInplace(BfVec *, BfMat const *);
void bfVecRecipInplace(BfVec *);
BfMat *bfVecGetGivensRotation(BfVec const *, BfSize, BfSize);
void bfVecPermute(BfVec *, BfPerm const *);
BfVec *bfVecConcat(BfVec const *, BfVec const *);
void bfVecSave(BfVec const *, char const *);
void bfVecDaxpy(BfVec *, BfReal, BfVec const *);
void bfVecDscal(BfVec *, BfReal);
BfSize bfVecGetSize(BfVec const *);

typedef struct BfVecVtable {
  __typeof__(&bfVecCopy) Copy;
  __typeof__(&bfVecDelete) Delete;
  __typeof__(&bfVecGetType) GetType;
  __typeof__(&bfVecInstanceOf) InstanceOf;
  __typeof__(&bfVecGetEltPtr) GetEltPtr;
  __typeof__(&bfVecGetSubvecCopy) GetSubvecCopy;
  __typeof__(&bfVecGetSubvecView) GetSubvecView;
  __typeof__(&bfVecGetSubvecViewConst) GetSubvecViewConst;
  __typeof__(&bfVecSetRange) SetRange;
  __typeof__(&bfVecSetMask) SetMask;
  __typeof__(&bfVecPrint) Print;
  __typeof__(&bfVecDist) Dist;
  __typeof__(&bfVecDistMax) DistMax;
  __typeof__(&bfVecNormMax) NormMax;
  __typeof__(&bfVecScaleByReal) ScaleByReal;
  __typeof__(&bfVecAddInplace) AddInplace;
  __typeof__(&bfVecMulInplace) MulInplace;
  __typeof__(&bfVecSolveInplace) SolveInplace;
  __typeof__(&bfVecRecipInplace) RecipInplace;
  __typeof__(&bfVecGetGivensRotation) GetGivensRotation;
  __typeof__(&bfVecPermute) Permute;
  __typeof__(&bfVecConcat) Concat;
  __typeof__(&bfVecSave) Save;
  __typeof__(&bfVecDaxpy) Daxpy;
  __typeof__(&bfVecDscal) Dscal;
} BfVecVtable;

/** Implementation: Vec */

struct BfVec {
  BfVecVtable *vtbl;
  BfVecProps props;
  BfSize size;
};

void bfVecInit(BfVec *vec, BfVecVtable *vtbl, BfSize size);
void bfVecDeinit(BfVec *vec);
BfVec *bfVecFromFile(char const *path, BfSize size, BfDtype dtype);
bool bfVecInstanceOf(BfVec const *vec, BfType type);
BfSize bfVecGetSize(BfVec const *vec);
