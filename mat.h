#pragma once

#include "def.h"
#include "dtype.h"
#include "vec.h"

enum BfMatProps {
  BF_MAT_PROP_NONE          = 0,
  BF_MAT_PROP_VIEW          = (1 << 0),
  BF_MAT_PROP_DIAGONAL      = (1 << 1),
  BF_MAT_PROP_TRANS         = (1 << 2),
  BF_MAT_PROP_CONJ_TRANS    = (1 << 3),
  BF_MAT_PROP_UNITARY       = (1 << 4),
  BF_MAT_PROP_SEMI_UNITARY  = (1 << 5)
};

typedef struct BfMat {
  enum BfDtypes dtype;
  enum BfMatProps props;
  BfSize numRows, numCols;
  BfSize rowStride, colStride;
  BfPtr data;
} BfMat;

BfMat bfGetUninitializedMat();

bool bfMatIsInitialized(BfMat const *A);

BfSize bfMatSize(BfMat const *A);
void bfMatNumBytes(BfMat const *A, BfSize *nbytes);

BfSize bfMatNumRows(BfMat const *A);
BfSize bfMatNumCols(BfMat const *A);

bool bfMatIsAligned(BfMat const *A);
BfSize bfMatRowStride(BfMat const *A);
BfSize bfMatColStride(BfMat const *A);

void bfFreeMat(BfMat *A);

void
bfInitEmptyMat(BfMat *A, enum BfDtypes dtype, enum BfMatProps props,
               BfSize numRows, BfSize numCols);

void
bfMatZeros(BfMat *A, enum BfDtypes dtype, BfSize numRows, BfSize numCols);

void
bfSaveMat(BfMat const *A, char const *path);

void
bfFillMatRandn(BfMat *A);

void
bfGetMatElt(BfMat const *A, BfSize i, BfSize j, BfPtr ptr);

void
bfGetMatEltPtr(BfMat const *A, BfSize i, BfSize j, BfPtr *ptr);

BfVec bfGetMatRow(BfMat const *A, BfSize i);

void
bfGetRowPtr(BfMat const *A, BfSize i, BfPtr *ptr);

void
bfSetMatRow(BfMat *A, BfSize i, void const *data);

void
bfCopyMatRow(BfMat const *A, BfSize i, BfMat *B, BfSize j);

BfMat bfGetMatRowRange(BfMat const *A, BfSize i0, BfSize i1);

BfMat bfGetMatColRange(BfMat const *A, BfSize j0, BfSize j1);

BfMat bfGetMatContSubblock(BfMat const *A, BfSize i0, BfSize i1,
                           BfSize j0, BfSize j1);

bool bfMatIsTransposed(BfMat const *A);

BfMat bfConjTrans(BfMat const *A);

void bfMatAddInplace(BfMat *A, BfMat const *B);

void
bfMatMul(BfMat const *A, BfMat const *B, BfMat *C);

void
bfMatSolve(BfMat const *A, BfMat const *B, BfMat *C);

void
bfInitEmptySvdMats(BfMat const *A, BfMat *U, BfMat *S, BfMat *Vt);

void
bfComputeMatSvd(BfMat const *A, BfMat *U, BfMat *S, BfMat *Vt);

void
bfMatLstSq(BfMat const *A, BfMat const *B, BfMat *C);
