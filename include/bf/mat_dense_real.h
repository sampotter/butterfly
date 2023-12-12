#pragma once

#include "mat_dense.h"

/** Interface: Mat */

BfMat *bfMatDenseRealCopy(BfMatDenseReal const *matDenseReal);
BfMat *bfMatDenseRealGetView(BfMat *mat);
BfMat *bfMatDenseRealSteal(BfMatDenseReal *matDenseReal);
BfVecReal *bfMatDenseRealGetRowView(BfMatDenseReal *matDenseReal, BfSize i);
BfVec *bfMatDenseRealGetColView(BfMat *mat, BfSize j);
void bfMatDenseRealDelete(BfMatDenseReal **matDenseReal);
BfMat *bfMatDenseRealZerosLike(BfMatDenseReal const *matDenseReal, BfSize numRows, BfSize numCols);
BfType bfMatDenseRealGetType(BfMat const *mat);
BfSize bfMatDenseRealNumBytes(BfMatDenseReal const *);
void bfMatDenseRealSave(BfMat const *mat, char const *path);
void bfMatDenseRealDump(BfMatDenseReal const *matDenseReal, FILE *fp);
void bfMatDenseRealPrint(BfMat const *mat, FILE *fp);
BfSize bfMatDenseRealGetNumRows(BfMatDenseReal const *matDenseReal);
BfSize bfMatDenseRealGetNumCols(BfMatDenseReal const *matDenseReal);
void bfMatDenseRealSetRow(BfMat *mat, BfSize i, BfVec const *row);
void bfMatDenseRealSetCol(BfMatDenseReal *matDenseReal, BfSize j, BfVec const *col);
BfMat *bfMatDenseRealGetRowRange(BfMat *mat, BfSize i0, BfSize i1);
BfMat *bfMatDenseRealGetRowRangeCopy(BfMatDenseReal const *matDenseReal, BfSize i0, BfSize i1);
BfMat *bfMatDenseRealGetColRangeCopy(BfMatDenseReal const *matDenseReal, BfSize j0, BfSize j1);
void bfMatDenseRealPermuteRows(BfMatDenseReal *matDenseReal, BfPerm const *perm);
void bfMatDenseRealPermuteCols(BfMatDenseReal *matDenseReal, BfPerm const *perm);
void bfMatDenseRealScaleRows(BfMatDenseReal *matDenseReal, BfVec const *vec);
BfMat *bfMatDenseRealMul(BfMatDenseReal const *matDenseReal, BfMat const *otherMat);
BfVec *bfMatDenseRealMulVec(BfMatDenseReal const *matDenseReal, BfVec const *vec);
void bfMatDenseRealMulInplace(BfMatDenseReal const *op1, BfMat const *op2, BfMat *res);
BfVec *bfMatDenseRealRmulVec(BfMatDenseReal const *matDenseReal, BfVec const *vec);
void bfMatDenseRealPrintBlocksDeep(BfMatDenseReal const *matDenseReal, FILE *fp, BfSize i0, BfSize j0, BfSize depth);
BfReal bfMatDenseRealNormMax(BfMatDenseReal const *matDenseReal);
BfReal bfMatDenseRealDistMax(BfMatDenseReal const *matDenseReal, BfMatDenseReal const *otherMatDenseReal);

/** Upcasting: MatDenseReal -> Mat */

BfMat *bfMatDenseRealToMat(BfMatDenseReal *matDenseReal);
BfMat const *bfMatDenseRealConstToMatConst(BfMatDenseReal const *matDenseReal);

/** Upcasting: MatDenseReal -> MatDense */

BfMatDense *bfMatDenseRealToMatDense(BfMatDenseReal *matDenseReal);
BfMatDense const *bfMatDenseRealConstToMatDenseConst(BfMatDenseReal const *matDenseReal);

/** Downcasting: Mat -> MatDenseReal */

BfMatDenseReal *bfMatToMatDenseReal(BfMat *mat);
BfMatDenseReal const *bfMatConstToMatDenseRealConst(BfMat const *mat);

/** Implementation: MatDenseReal */

struct BfMatDenseReal {
  BfMatDense super;
  BfReal *data;
};

BfMatDenseReal *bfMatDenseRealNew(void);
BfMatDenseReal *bfMatDenseRealNewViewFromPtr(BfSize numRows, BfSize numCols, BfReal *data, BfSize rowStride, BfSize colStride);
BfMatDenseReal const *bfMatDenseRealNewViewConstFromPtr(BfSize numRows, BfSize numCols, BfReal const *data, BfSize rowStride, BfSize colStride);
BfMatDenseReal *bfMatDenseRealNewWithValue(BfSize numRows, BfSize numCols, BfReal value);
BfMatDenseReal *bfMatDenseRealNewFromMatrix(BfMat const *mat);
BfMatDenseReal *bfMatDenseRealNewZeros(BfSize numRows, BfSize numCols);
BfMatDenseReal *bfMatDenseRealFromFile(char const *path, BfSize numRows, BfSize numCols);
BfMatDenseReal *bfMatDenseRealNewFromRealArray(BfRealArray *realArray, BfSize numRows, BfSize numCols, BfPolicy policy);
BfMatDenseReal *bfMatDenseRealNewFromCsv(char const *path);
void bfMatDenseRealInit(BfMatDenseReal *mat, BfSize numRows, BfSize numCols);
void bfMatDenseRealInitViewFromPtr(BfMatDenseReal *matDenseReal, BfSize numRows, BfSize numCols, BfReal *data, BfSize rowStride, BfSize colStride);
void bfMatDenseRealInitViewConstFromPtr(BfMatDenseReal const *matDenseReal, BfSize numRows, BfSize numCols, BfReal const *data, BfSize rowStride, BfSize colStride);
void bfMatDenseRealInitCopy(BfMatDenseReal *mat, BfMatDenseReal const *otherMat);
void bfMatDenseRealInitWithValue(BfMatDenseReal *mat, BfSize numRows, BfSize numCols, BfReal fillValue);
void bfMatDenseRealInitZeros(BfMatDenseReal *mat, BfSize numRows, BfSize numCols);
void bfMatDenseRealInitFromPtr(BfMatDenseReal *matDenseReal, BfSize numRows, BfSize numCols, BfReal *data, BfPolicy policy);
void bfMatDenseRealInitFromRealArray(BfMatDenseReal *matDenseReal, BfRealArray *realArray, BfSize numRows, BfSize numCols, BfPolicy policy);
void bfMatDenseRealDeinit(BfMatDenseReal *mat);
void bfMatDenseRealDealloc(BfMatDenseReal **mat);
void bfMatDenseRealDeinitAndDealloc(BfMatDenseReal **mat);
void bfMatDenseRealSvd(BfMatDenseReal const *mat, BfMatDenseReal **UPtr,
                       BfMatDiagReal **SPtr, BfMatDenseReal **VTPtr);
void bfMatDenseRealSetBlock(BfMatDenseReal *matDenseReal, BfSize i0, BfSize i1, BfSize j0, BfSize j1, BfMatDenseReal const *block);
