#pragma once

#include "mat.h"

/** Interface: Mat */

BfMat *bfMatDenseComplexGetView(BfMat *mat);
BfMat *bfMatDenseComplexCopy(BfMat const *mat);
BfMat *bfMatDenseComplexSteal(BfMatDenseComplex *matDenseComplex);
BfVec *bfMatDenseComplexGetRowCopy(BfMat const *mat, BfSize i);
BfVec *bfMatDenseComplexGetRowView(BfMat *mat, BfSize i);
BfVec *bfMatDenseComplexGetColView(BfMat *mat, BfSize j);
BfVec *bfMatDenseComplexGetColRangeView(BfMat *mat, BfSize i0, BfSize i1, BfSize j);
void bfMatDenseComplexDelete(BfMatDenseComplex **matDenseComplex);
BfMat *bfMatDenseComplexEmptyLike(BfMat const *mat, BfSize numRows, BfSize numCols);
BfMat *bfMatDenseComplexZerosLike(BfMat const *mat, BfSize numRows, BfSize numCols);
BfType bfMatDenseComplexGetType(BfMat const *mat);
BfSize bfMatDenseComplexNumBytes(BfMatDenseComplex const *matDenseComplex);
void bfMatDenseComplexSave(BfMat const *mat, char const *path);
void bfMatDenseComplexPrint(BfMat const *mat, FILE *fp);
BfSize bfMatDenseComplexGetNumRows(BfMatDenseComplex const *matDenseComplex);
BfSize bfMatDenseComplexGetNumCols(BfMatDenseComplex const *matDenseComplex);
void bfMatDenseComplexSetRow(BfMat *mat, BfSize i, BfVec const *rowVec);
void bfMatDenseComplexSetCol(BfMat *mat, BfSize j, BfVec const *vec);
void bfMatDenseComplexSetColRange(BfMat *mat, BfSize j, BfSize i0, BfSize i1, BfVec const *vec);
BfMat *bfMatDenseComplexGetRowRange(BfMat *mat, BfSize i0, BfSize i1);
BfMat *bfMatDenseComplexGetRowRangeCopy(BfMatDenseComplex const *matDenseComplex, BfSize i0, BfSize i1);
BfMat const *bfMatDenseComplexGetRowRangeConst(BfMatDenseComplex const *matDenseComplex, BfSize i0, BfSize i1);
BfMat *bfMatDenseComplexGetColRange(BfMatDenseComplex *matDenseComplex, BfSize j0, BfSize j1);
BfMat const *bfMatDenseComplexGetColRangeConst(BfMatDenseComplex const *matDenseComplex, BfSize j0, BfSize j1);
void bfMatDenseComplexSetRowRange(BfMat *mat, BfSize i0, BfSize i1, BfMat const *rows);
void bfMatDenseComplexPermuteRows(BfMat *mat, BfPerm const *perm);
BfVec *bfMatDenseComplexColDists(BfMat const *mat, BfMat const *otherMat);
BfVec *bfMatDenseComplexColDots(BfMat const *mat, BfMat const *otherMat);
BfVec *bfMatDenseComplexColNorms(BfMat const *mat);
void bfMatDenseComplexScale(BfMatDenseComplex *matDenseComplex, BfComplex scalar);
void bfMatDenseComplexScaleCols(BfMat *mat, BfVec const *vec);
void bfMatDenseComplexAddInplace(BfMat *mat, BfMat const *otherMat);
BfMat *bfMatDenseComplexSub(BfMat const *mat, BfMat const *otherMat);
void bfMatDenseComplexSubInplace(BfMat *mat, BfMat const *otherMat);
BfMat *bfMatDenseComplexMul(BfMat const *op1, BfMat const *op2);
BfVec *bfMatDenseComplexMulVec(BfMat const *mat, BfVec const *vec);
BfMat *bfMatDenseComplexRmul(BfMatDenseComplex const *matDenseComplex, BfMat const *mat);
BfMat *bfMatDenseComplexSolve(BfMatDenseComplex const *matDenseComplex, BfMat const *otherMat);
BfMat *bfMatDenseComplexSolveLU(BfMat const *A, BfMat const *B);
BfMat *bfMatDenseComplexLstSq(BfMat const *lhs, BfMat const *rhs);
bool bfMatDenseComplexIsUpperTri(BfMat const *mat);
BfVec *bfMatDenseComplexBackwardSolveVec(BfMat const *mat, BfVec const *vec);
void bfMatDenseComplexNegate(BfMat *mat);
void bfMatDenseComplexPrintBlocksDeep(BfMatDenseComplex const *matDenseComplex, FILE *fp, BfSize i0, BfSize j0, BfSize depth);
BfMat *bfMatDenseComplexGetBlockView(BfMatDenseComplex *mat, BfSize i0, BfSize i1, BfSize j0, BfSize j1);
BfLu *bfMatDenseComplexGetLu(BfMatDenseComplex const *mat);
void bfMatDenseComplexDivideCols(BfMatDenseComplex *matDenseComplex, BfVec const *vec);
void bfMatDenseComplexTranspose(BfMatDenseComplex *matDenseComplex);

/** Upcasting: */

BfMat *bfMatDenseComplexToMat(BfMatDenseComplex *matDenseComplex);
BfMat const *bfMatDenseComplexConstToMatConst(BfMatDenseComplex const *matDenseComplex);

/** Downcasting: */

BfMatDenseComplex *bfMatToMatDenseComplex(BfMat *mat);
BfMatDenseComplex const *bfMatConstToMatDenseComplexConst(BfMat const *mat);

/** Implementation: MatDenseComplex */

struct BfMatDenseComplex {
  BfMat super;
  BfSize rowStride;
  BfSize colStride;
  BfComplex *data;

  /* A pointer optionally used to signal that this array is backed by
   * a PyArray. This is the case if `pyArray != NULL`, in which case
   * the lifecycle of this instance will take care of any extra
   * reference counting to make sure `pyArray` isn't deleted while
   * this instance is alive. */
  BfPtr *pyArray;
};

BfMatDenseComplex *bfMatDenseComplexNew(void);
BfMatDenseComplex *bfMatDenseComplexNewViewFromPtr(BfSize numRows, BfSize numCols, BfComplex *data);
BfMatDenseComplex *bfMatDenseComplexNewRandn(BfSize numRows, BfSize numCols);
BfMatDenseComplex *bfMatDenseComplexNewViewFromPyArray(BfPtr *arr);
BfMatDenseComplex *bfMatDenseComplexZeros(BfSize numRows, BfSize numCols);
BfMatDenseComplex *bfMatDenseComplexFromFile(char const *path, BfSize numRows, BfSize numCols);
BfMatDenseComplex *bfMatDenseComplexNewIdentity(BfSize numRows, BfSize numCols);
void bfMatDenseComplexInit(BfMatDenseComplex *mat, BfSize numRows, BfSize numCols);
void bfMatDenseComplexInitViewFromPtr(BfMatDenseComplex *matDenseComplex, BfSize numRows, BfSize numCols, BfComplex *data);
void bfMatDenseComplexInitRandn(BfMatDenseComplex *matDenseComplex, BfSize numRows, BfSize numCols);
void bfMatDenseComplexInitViewFromPyArray(BfMatDenseComplex *matDenseComplex, BfPtr *pyArray);
void bfMatDenseComplexInitIdentity(BfMatDenseComplex *matDenseComplex, BfSize numRows, BfSize numCols);
void bfMatDenseComplexDeinit(BfMatDenseComplex *mat);
void bfMatDenseComplexDealloc(BfMatDenseComplex **mat);
void bfMatDenseComplexDeinitAndDealloc(BfMatDenseComplex **mat);
void bfMatDenseComplexSet(BfMatDenseComplex *dst, BfMatDenseComplex const *src);
BfSize bfMatDenseComplexGetRowStride(BfMatDenseComplex const *mat);
BfSize bfMatDenseComplexGetColStride(BfMatDenseComplex const *mat);
void bfMatDenseComplexSvd(BfMatDenseComplex const *mat, BfMatDenseComplex *U,
                          BfMatDiagReal *S, BfMatDenseComplex *VH);
void bfMatDenseComplexCooComplexAddInplace(BfMatDenseComplex *op1,
                                           BfMatCooComplex const *op2);
void bfMatDenseComplexCooRealAddInplace(BfMatDenseComplex *op1,
                                        BfMatCooReal const *op2);
void bfMatDenseComplexDenseComplexAddInplace(BfMatDenseComplex *op1,
                                             BfMatDenseComplex const *op2);
void bfMatDenseComplexDiagRealAddInplace(BfMatDenseComplex *op1,
                                         BfMatDiagReal const *op2);
BfMatDenseComplex *
bfMatDenseComplexDenseComplexSub(BfMatDenseComplex const *op1,
                                 BfMatDenseComplex const *op2);
BfMatDenseComplex *
bfMatDenseComplexDenseComplexMul(BfMatDenseComplex const *op1,
                                 BfMatDenseComplex const *op2);
BfMatDenseComplex *
bfMatDenseComplexDenseComplexLstSq(BfMatDenseComplex const *lhs,
                                   BfMatDenseComplex const *rhs);
bool bfMatDenseComplexIsFinite(BfMatDenseComplex const *mat);

void bf_zmat_add_diag(BfMatDenseComplex *mat, BfReal value);
BfMatDenseComplex *bfMatDenseComplexDenseComplexSolve(BfMatDenseComplex const *A, BfMatDenseComplex const *B);
void bfMatDenseComplexDenseRealScaleCols(BfMatDenseComplex *mat, BfMatDenseReal const *otherMat);
BfComplex *bfMatDenseComplexGetDataPtr(BfMatDenseComplex *matDenseComplex);
BfComplex const *bfMatDenseComplexGetDataConstPtr(BfMatDenseComplex const *matDenseComplex);
void bfMatDenseComplexSetBlock(BfMatDenseComplex *matDenseComplex, BfSize i0, BfSize j0, BfMat const *mat);
BfMatDenseReal *bfMatDenseComplexGetReView(BfMatDenseComplex *matDenseComplex);
BfMatDenseReal const *bfMatDenseComplexGetReViewConst(BfMatDenseComplex const *matDenseComplex);
BfMatDenseReal *bfMatDenseComplexGetImView(BfMatDenseComplex *matDenseComplex);
BfMatDenseReal const *bfMatDenseComplexGetImViewConst(BfMatDenseComplex const *matDenseComplex);
BfComplex *bfMatDenseComplexGetRowPtr(BfMatDenseComplex *matDenseComplex, BfSize i);
BfComplex const *bfMatDenseComplexGetRowConstPtr(BfMatDenseComplex const *matDenseComplex, BfSize i);
