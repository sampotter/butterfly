#pragma once

#include "mat_dense.h"

/** Interface: Mat */

BfMat *bfMatDenseRealCopy(BfMat const *mat);
BfMat *bfMatDenseRealGetView(BfMat *mat);
BfVec *bfMatDenseRealGetColView(BfMat *mat, BfSize j);
void bfMatDenseRealDelete(BfMat **mat);
BfType bfMatDenseRealGetType(BfMat const *mat);
void bfMatDenseRealSave(BfMat const *mat, char const *path);
void bfMatDenseRealPrint(BfMat const *mat, FILE *fp);
BfSize bfMatDenseRealGetNumRows(BfMat const *mat);
BfSize bfMatDenseRealGetNumCols(BfMat const *mat);
void bfMatDenseRealSetRow(BfMat *mat, BfSize i, BfVec const *row);
void bfMatDenseRealSetCol(BfMat *mat, BfSize j, BfVec const *col);
BfMat *bfMatDenseRealGetColRangeCopy (BfMat const *mat, BfSize j0, BfSize j1);
void bfMatDenseRealPermuteRows(BfMat *mat, BfPerm const *perm);

struct BfMatDenseReal {
  BfMatDense super;
  BfReal *data;
};

BfMat *bfMatDenseRealToMat(BfMatDenseReal *matDenseReal);
BfMat const *bfMatDenseRealConstToMatConst(BfMatDenseReal const *matDenseReal);

BfMatDense *bfMatDenseRealToMatDense(BfMatDenseReal *matDenseReal);
BfMatDense const *bfMatDenseRealConstToMatDenseConst(BfMatDenseReal const *matDenseReal);

BfMatDenseReal *bfMatToMatDenseReal(BfMat *mat);
BfMatDenseReal const *bfMatConstToMatDenseRealConst(BfMat const *mat);

BfMatDenseReal *bfMatDenseRealNew();
BfMatDenseReal *bfMatDenseRealFromFile(char const *path, BfSize numRows, BfSize numCols);
void bfMatDenseRealInit(BfMatDenseReal *mat, BfSize numRows, BfSize numCols);
void bfMatDenseRealInitWithValue(BfMatDenseReal *mat, BfSize numRows,
                                 BfSize numCols, BfReal fillValue);
void bfMatDenseRealDeinit(BfMatDenseReal *mat);
void bfMatDenseRealDealloc(BfMatDenseReal **mat);
void bfMatDenseRealDeinitAndDealloc(BfMatDenseReal **mat);
void bfMatDenseRealSvd(BfMatDenseReal const *mat, BfMatDenseReal *U,
                       BfMatDiagReal *S, BfMatDenseReal *VH);
