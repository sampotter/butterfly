#pragma once

#include "mat.h"

/* Sparse matrix stored in compressed sparse row format. */
struct BfMatCsrReal {
  BfMat super;

  /* Array of length `super.numRows + 1`, where `rowptr[0] == 0` and
   * `rowptr[super.numRows]` gives the number of nonzero entries in
   * the matrix. Other entries gives the offset into `colind` and
   * `data` for the start of each row. */
  BfSize *rowptr;

  /* The column of each entry in `data`. */
  BfSize *colind;

  /* Nonzero (and explicit zero) entries of matrix. */
  BfReal *data;
};

#define INTERFACE BF_INTERFACE_Mat
BF_DECLARE_INTERFACE(MatCsrReal)
#undef INTERFACE

/* Upcasting: */
BfMat *bfMatCsrRealToMat(BfMatCsrReal *matCsrReal);
BfMat const *bfMatCsrRealConstToMatConst(BfMatCsrReal const *matCsrReal);

/* Downcasting: */
BfMatCsrReal *bfMatToMatCsrReal(BfMat *mat);
BfMatCsrReal const *bfMatConstToMatCsrRealConst(BfMat const *mat);

BfMatCsrReal *bfMatCsrRealNew();
void bfMatCsrRealInit(BfMatCsrReal *mat, BfSize numRows, BfSize numCols, BfSize const *rowptr, BfSize const *colind, BfReal const *data);
void bfMatCsrRealDeinit(BfMatCsrReal *mat);
void bfMatCsrRealDealloc(BfMatCsrReal **mat);
void bfMatCsrRealDeinitAndDealloc(BfMatCsrReal **mat);
void bfMatCsrRealDump(BfMatCsrReal const *mat, char const *rowptrPath, char const *colindPath,
                      char const *dataPath);
