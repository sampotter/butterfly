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
  BF_MAT_PROP_SEMI_UNITARY  = (1 << 5),
  BF_MAT_PROP_SPARSE_CSR    = (1 << 6)
};

typedef struct BfMat {
  enum BfDtypes dtype;
  enum BfMatProps props;
  BfSize numRows, numCols;
  BfSize rowStride, colStride;
  BfPtr data;
} BfMat;

BfSize bfMatSize(BfMat const *A);
enum BfError bfMatNumBytes(BfMat const *A, BfSize *nbytes);

BfSize bfMatNumRows(BfMat const *A);
BfSize bfMatNumCols(BfMat const *A);

bool bfMatIsAligned(BfMat const *A);
BfSize bfMatRowStride(BfMat const *A);
BfSize bfMatColStride(BfMat const *A);

enum BfError bfFreeMat(BfMat *A);

enum BfError
bfInitEmptyMat(BfMat *A, enum BfDtypes dtype, enum BfMatProps props,
               BfSize numRows, BfSize numCols);

enum BfError
bfSaveMat(BfMat const *A, char const *path);

enum BfError
bfFillMatRandn(BfMat *A);

enum BfError
bfGetMatElt(BfMat const *A, BfSize i, BfSize j, BfPtr ptr);

enum BfError
bfGetMatEltPtr(BfMat const *A, BfSize i, BfSize j, BfPtr *ptr);

BfVec bfGetMatRow(BfMat const *A, BfSize i);

enum BfError
bfGetRowPtr(BfMat const *A, BfSize i, BfPtr *ptr);

enum BfError
bfSetMatRow(BfMat *A, BfSize i, void const *data);

enum BfError
bfCopyMatRow(BfMat const *A, BfSize i, BfMat *B, BfSize j);

BfMat bfGetMatColRange(BfMat const *A, BfSize j0, BfSize j1);

BfMat bfGetMatContSubblock(BfMat const *A, BfSize i0, BfSize i1,
                           BfSize j0, BfSize j1);

bool bfMatIsTransposed(BfMat const *A);

BfMat bfConjTrans(BfMat const *A);

enum BfError
bfMatMul(BfMat const *A, BfMat const *B, BfMat *C);

enum BfError
bfMatSolve(BfMat const *A, BfMat const *B, BfMat *C);

enum BfError
bfInitEmptySvdMats(BfMat const *A, BfMat *U, BfMat *S, BfMat *Vt);

enum BfError
bfComputeMatSvd(BfMat const *A, BfMat *U, BfMat *S, BfMat *Vt);

enum BfError
bfComputePinv(BfMat const *A, BfReal atol, BfReal rtol, BfMat *pinv);

enum BfError
bfMatLstSq(BfMat const *A, BfMat const *B, BfMat *C);
