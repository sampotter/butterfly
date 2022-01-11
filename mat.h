#pragma once

#include "def.h"

typedef enum BfMatTypes {
  BF_MAT_TYPE_MAT,
  BF_MAT_TYPE_DENSE_COMPLEX,
  BF_MAT_TYPE_DIAG_REAL,
  BF_MAT_TYPE_COUNT
} BfMatType;

typedef struct BfMat BfMat;
typedef struct BfMatBlockCoo BfMatBlockCoo;
typedef struct BfMatDenseComplex BfMatDenseComplex;
typedef struct BfMatDiagReal BfMatDiagReal;

typedef struct BfMatVtable {
  void (*deinit)(BfMat *);
  void (*delete)(BfMat **);
  void (*deinitAndDelete)(BfMat **);
  BfMatType (*getType)(BfMat const *);
  BfSize (*numBytes)(BfMat const *);
  void (*save)(BfMat const *, char const *);
  BfMat *(*mul)(BfMat const *, BfMat const *);
  BfMat *(*lstSq)(BfMat const *, BfMat const *);
} BfMatVtable;

typedef enum BfMatProps {
  BF_MAT_PROPS_NONE = 0,
  BF_MAT_PROPS_VIEW = 1 << 0,
  BF_MAT_PROPS_TRANS = 1 << 1,
  BF_MAT_PROPS_CONJ = 1 << 2,
  BF_MAT_PROPS_ORTHO = 1 << 3
} BfMatProps;

struct BfMat {
  BfMatVtable *vtbl;
  BfMatProps props;
  BfSize numRows;
  BfSize numCols;
};

void bfMatInit(BfMat *mat, BfMatVtable *vtbl, BfSize numRows, BfSize numCols);
void bfMatDeinit(BfMat *mat);
void bfMatDelete(BfMat **mat);
void bfMatDeinitAndDelete(BfMat **mat);
BfMatType bfMatGetType(BfMat const *mat);
BfSize bfMatNumBytes(BfMat const *mat);
void bfMatSave(BfMat const *mat, char const *path);
bool bfMatIsTransposed(BfMat const *mat);
BfSize bfMatGetNumRows(BfMat const *mat);
BfSize bfMatGetNumCols(BfMat const *mat);
BfMat *bfMatConjTrans(BfMat *mat);
BfMat *bfMatGetRowRange(BfMat *mat, BfSize i0, BfSize i1);
BfMat *bfMatGetColRange(BfMat *mat, BfSize j0, BfSize j1);
BfMat *bfMatMul(BfMat const *lhs, BfMat const *rhs);
BfMat *bfMatLstSq(BfMat const *lhs, BfMat const *rhs);
