#pragma once

#include <stdio.h>

#include "backends.h"
#include "def.h"
#include "dtype.h"
#include "perm.h"
#include "types.h"
#include "vec.h"

typedef enum BfMatProps {
  BF_MAT_PROPS_NONE = 0,
  BF_MAT_PROPS_VIEW = 1 << 0,
  BF_MAT_PROPS_TRANS = 1 << 1,
  BF_MAT_PROPS_CONJ = 1 << 2,
  BF_MAT_PROPS_ORTHO = 1 << 3,
  BF_MAT_PROPS_LOWER_TRI = 1 << 4,
  BF_MAT_PROPS_UPPER_TRI = 1 << 5,
  BF_MAT_PROPS_UNIT = 1 << 6
} BfMatProps;

static BfMatProps const BF_MAT_PROPS_TRI = BF_MAT_PROPS_LOWER_TRI | BF_MAT_PROPS_UPPER_TRI;

/** Interface: Mat */

BfMat *bfMatGetView(BfMat *);
BfMat *bfMatCopy(BfMat const *);
BfMat *bfMatSteal(BfMat *);
BfVec *bfMatGetRowCopy(BfMat const *, BfSize);
BfVec *bfMatGetRowView(BfMat *, BfSize);
BfVec *bfMatGetColView(BfMat *, BfSize);
BfVec *bfMatGetColRangeView(BfMat *, BfSize, BfSize, BfSize);
void bfMatDelete(BfMat **);
BfMat *bfMatEmptyLike(BfMat const *, BfSize, BfSize);
BfMat *bfMatZerosLike(BfMat const *, BfSize, BfSize);
BfType bfMatGetType(BfMat const *);
BfSize bfMatNumBytes(BfMat const *);
void bfMatSave(BfMat const *, char const *);
void bfMatPrint(BfMat const *, FILE *);
BfSize bfMatGetNumRows(BfMat const *);
BfSize bfMatGetNumCols(BfMat const *);
void bfMatSetRow(BfMat *, BfSize, BfVec const *);
void bfMatSetCol(BfMat *, BfSize, BfVec const *);
void bfMatSetColRange(BfMat *, BfSize, BfSize, BfSize, BfVec const *);
BfMat *bfMatGetRowRange(BfMat *, BfSize, BfSize);
BfMat *bfMatGetColRange(BfMat *, BfSize, BfSize);
BfMat *bfMatGetRowRangeCopy(BfMat const *, BfSize, BfSize);
BfMat *bfMatGetColRangeCopy (BfMat const *, BfSize, BfSize);
void bfMatSetRowRange(BfMat *, BfSize, BfSize, BfMat const *);
void bfMatPermuteRows(BfMat *, BfPerm const *);
void bfMatPermuteCols(BfMat *, BfPerm const *);
BfVec *bfMatRowDists(BfMat const *, BfMat const *);
BfVec *bfMatColDists(BfMat const *, BfMat const *);
BfVec *bfMatColDots(BfMat const *, BfMat const *);
BfVec *bfMatColNorms(BfMat const *);
void bfMatScale(BfMat *, BfComplex);
void bfMatScaleRows(BfMat *, BfVec const *);
void bfMatScaleCols(BfMat *, BfVec const *);
BfVec *bfMatSumCols(BfMat const *);
BfMat *bfMatAdd(BfMat const *, BfMat const *);
void bfMatAddInplace(BfMat *, BfMat const *);
void bfMatAddDiag(BfMat *, BfMat const *);
BfMat *bfMatSub(BfMat const *, BfMat const *);
void bfMatSubInplace(BfMat *, BfMat const *);
BfMat *bfMatMul(BfMat const *, BfMat const *);
BfVec *bfMatMulVec(BfMat const *, BfVec const *);
void bfMatMulInplace(BfMat *, BfMat const *);
BfVec *bfMatRmulVec(BfMat const *, BfVec const *);
BfMat *bfMatSolve(BfMat const *, BfMat const *);
BfMat *bfMatSolveLU(BfMat const *, BfMat const *);
BfMat *bfMatLstSq(BfMat const *, BfMat const *);
bool bfMatIsUpperTri(BfMat const *);
BfVec *bfMatForwardSolveVec(BfMat const *, BfVec const *);
BfVec *bfMatBackwardSolveVec(BfMat const *, BfVec const *);
bool bfMatIsZero(BfMat const *);
void bfMatNegate(BfMat *);
BfMat *bfMatToType(BfMat const *, BfType);
BfMat *bfMatCholesky(BfMat const *);
BfSizeArray *bfMatGetNonzeroColumnRanges(BfMat const *);
void bfMatPrintBlocksDeep(BfMat const *, FILE *, BfSize, BfSize, BfSize);
BfMat *bfMatGetBlockView(BfMat *mat, BfSize, BfSize, BfSize, BfSize);
BfLu *bfMatGetLu(BfMat const *mat);
BfMat *bfMatGetInverse(BfMat const *mat);
void bfMatDivideCols(BfMat *, BfVec const *);

typedef struct BfMatVtable {
  __typeof__(&bfMatGetView) GetView;
  __typeof__(&bfMatCopy) Copy;
  __typeof__(&bfMatSteal) Steal;
  __typeof__(&bfMatGetRowCopy) GetRowCopy;
  __typeof__(&bfMatGetRowView) GetRowView;
  __typeof__(&bfMatGetColView) GetColView;
  __typeof__(&bfMatGetColRangeView) GetColRangeView;
  __typeof__(&bfMatDelete) Delete;
  __typeof__(&bfMatEmptyLike) EmptyLike;
  __typeof__(&bfMatZerosLike) ZerosLike;
  __typeof__(&bfMatGetType) GetType;
  __typeof__(&bfMatNumBytes) NumBytes;
  __typeof__(&bfMatSave) Save;
  __typeof__(&bfMatPrint) Print;
  __typeof__(&bfMatGetNumRows) GetNumRows;
  __typeof__(&bfMatGetNumCols) GetNumCols;
  __typeof__(&bfMatSetRow) SetRow;
  __typeof__(&bfMatSetCol) SetCol;
  __typeof__(&bfMatSetColRange) SetColRange;
  __typeof__(&bfMatGetRowRange) GetRowRange;
  __typeof__(&bfMatGetColRange) GetColRange;
  __typeof__(&bfMatGetRowRangeCopy) GetRowRangeCopy;
  __typeof__(&bfMatGetColRangeCopy ) GetColRangeCopy ;
  __typeof__(&bfMatSetRowRange) SetRowRange;
  __typeof__(&bfMatPermuteRows) PermuteRows;
  __typeof__(&bfMatPermuteCols) PermuteCols;
  __typeof__(&bfMatRowDists) RowDists;
  __typeof__(&bfMatColDists) ColDists;
  __typeof__(&bfMatColDots) ColDots;
  __typeof__(&bfMatColNorms) ColNorms;
  __typeof__(&bfMatScale) Scale;
  __typeof__(&bfMatScaleRows) ScaleRows;
  __typeof__(&bfMatScaleCols) ScaleCols;
  __typeof__(&bfMatSumCols) SumCols;
  __typeof__(&bfMatAdd) Add;
  __typeof__(&bfMatAddInplace) AddInplace;
  __typeof__(&bfMatAddDiag) AddDiag;
  __typeof__(&bfMatSub) Sub;
  __typeof__(&bfMatSubInplace) SubInplace;
  __typeof__(&bfMatMul) Mul;
  __typeof__(&bfMatMulVec) MulVec;
  __typeof__(&bfMatMulInplace) MulInplace;
  __typeof__(&bfMatRmulVec) RmulVec;
  __typeof__(&bfMatSolve) Solve;
  __typeof__(&bfMatSolveLU) SolveLU;
  __typeof__(&bfMatLstSq) LstSq;
  __typeof__(&bfMatIsUpperTri) IsUpperTri;
  __typeof__(&bfMatForwardSolveVec) ForwardSolveVec;
  __typeof__(&bfMatBackwardSolveVec) BackwardSolveVec;
  __typeof__(&bfMatIsZero) IsZero;
  __typeof__(&bfMatNegate) Negate;
  __typeof__(&bfMatToType) ToType;
  __typeof__(&bfMatCholesky) Cholesky;
  __typeof__(&bfMatGetNonzeroColumnRanges) GetNonzeroColumnRanges;
  __typeof__(&bfMatPrintBlocksDeep) PrintBlocksDeep;
  __typeof__(&bfMatGetBlockView) GetBlockView;
  __typeof__(&bfMatGetLu) GetLu;
  __typeof__(&bfMatGetInverse) GetInverse;
  __typeof__(&bfMatDivideCols) DivideCols;
} BfMatVtable;

/** Implementation: Mat */

struct BfMat {
  BfMatVtable *vtbl;
  BfMatProps props;
  BfSize numRows;
  BfSize numCols;
#if BF_DEBUG
  void *aux; /* pointer to extra user-defined data for purposes of debugging */
#endif
};

BfMat *bfMatGet(BfMat *mat, BfPolicy policy);
void bfMatInvalidate(BfMat *mat);
void bfMatInit(BfMat *mat, BfMatVtable *vtbl, BfSize numRows, BfSize numCols);
void bfMatDeinit(BfMat *mat);
BfMat *bfMatFromFile(char const *path, BfSize numRows, BfSize numCols, BfDtype dtype);
bool bfMatInstanceOf(BfMat const *mat, BfType type);
bool bfMatIsTransposed(BfMat const *mat);
BfMat *bfMatTrans(BfMat *mat);
BfMat *bfMatConjTrans(BfMat *mat);
bool bfMatIsBlock(BfMat const *mat);
bool bfMatIsView(BfMat const *mat);
