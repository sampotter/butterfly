#pragma once

#include <stdio.h>

#include "backends.h"
#include "def.h"
#include "dtype.h"
#include "interface.h"
#include "perm.h"
#include "types.h"
#include "vec.h"

typedef enum BfMatProps {
  BF_MAT_PROPS_NONE = 0,
  BF_MAT_PROPS_VIEW = 1 << 0,
  BF_MAT_PROPS_TRANS = 1 << 1,
  BF_MAT_PROPS_CONJ = 1 << 2,
  BF_MAT_PROPS_ORTHO = 1 << 3
} BfMatProps;

#define BF_INTERFACE_Mat(Type, Subtype, _)                                     \
  _(Type, Subtype, BfMat *, Copy, BfMat const *)                               \
  _(Type, Subtype, BfMat *, GetView, BfMat *)                                  \
  _(Type, Subtype, BfVec *, GetRowCopy, BfMat const *, BfSize)                 \
  _(Type, Subtype, BfVec *, GetRowView, BfMat *, BfSize)                       \
  _(Type, Subtype, BfVec *, GetColView, BfMat *, BfSize)                       \
  _(Type, Subtype, BfVec *, GetColRangeView, BfMat *, BfSize, BfSize, BfSize)  \
  _(Type, Subtype, void, Delete, BfMat **)                                     \
  _(Type, Subtype, BfMat *, EmptyLike, BfMat const *, BfSize, BfSize)          \
  _(Type, Subtype, BfMat *, ZerosLike, BfMat const *, BfSize, BfSize)          \
  _(Type, Subtype, BfType, GetType, BfMat const *)                             \
  _(Type, Subtype, BfSize, NumBytes, BfMat const *)                            \
  _(Type, Subtype, void, Save, BfMat const *, char const *)                    \
  _(Type, Subtype, void, Print, BfMat const *, FILE *)                         \
  _(Type, Subtype, BfSize, GetNumRows, BfMat const *)                          \
  _(Type, Subtype, BfSize, GetNumCols, BfMat const *)                          \
  _(Type, Subtype, void, SetRow, BfMat *, BfSize, BfVec const *)               \
  _(Type, Subtype, void, SetCol, BfMat *, BfSize, BfVec const *)               \
  _(Type, Subtype, void, SetColRange, BfMat *, BfSize, BfSize, BfSize, BfVec const *) \
  _(Type, Subtype, BfMat *, GetRowRange, BfMat *, BfSize, BfSize)              \
  _(Type, Subtype, BfMat *, GetColRange, BfMat *, BfSize, BfSize)              \
  _(Type, Subtype, BfMat *, GetRowRangeCopy, BfMat const *, BfSize, BfSize)    \
  _(Type, Subtype, BfMat *, GetColRangeCopy , BfMat const *, BfSize, BfSize)   \
  _(Type, Subtype, void, SetRowRange, BfMat *, BfSize, BfSize, BfMat const *)  \
  _(Type, Subtype, void, PermuteRows, BfMat *, BfPerm const *)                 \
  _(Type, Subtype, void, PermuteCols, BfMat *, BfPerm const *)                 \
  _(Type, Subtype, BfVec *, RowDists, BfMat const *, BfMat const *)            \
  _(Type, Subtype, BfVec *, ColDists, BfMat const *, BfMat const *)            \
  _(Type, Subtype, BfVec *, ColDots, BfMat const *, BfMat const *)             \
  _(Type, Subtype, BfVec *, ColNorms, BfMat const *)                           \
  _(Type, Subtype, void, ScaleRows, BfMat *, BfVec const *)                    \
  _(Type, Subtype, void, ScaleCols, BfMat *, BfVec const *)                    \
  _(Type, Subtype, BfVec *, SumCols, BfMat const *)                            \
  _(Type, Subtype, void, AddInplace, BfMat *, BfMat const *)                   \
  _(Type, Subtype, void, AddDiag, BfMat *, BfMat const *)                      \
  _(Type, Subtype, BfMat *, Sub, BfMat const *, BfMat const *)                 \
  _(Type, Subtype, void, SubInplace, BfMat *, BfMat const *)                   \
  _(Type, Subtype, BfMat *, Mul, BfMat const *, BfMat const *)                 \
  _(Type, Subtype, BfVec *, MulVec, BfMat const *, BfVec const *)              \
  _(Type, Subtype, void, MulInplace, BfMat *, BfMat const *)                   \
  _(Type, Subtype, BfMat *, SolveLU, BfMat const *, BfMat const *)             \
  _(Type, Subtype, BfMat *, LstSq, BfMat const *, BfMat const *)               \
  _(Type, Subtype, bool, IsUpperTri, BfMat const *)                            \
  _(Type, Subtype, BfVec *, ForwardSolveVec, BfMat const *, BfVec const *)     \
  _(Type, Subtype, BfVec *, BackwardSolveVec, BfMat const *, BfVec const *)    \
  _(Type, Subtype, bool, IsZero, BfMat const *)                                \
  _(Type, Subtype, void, Negate, BfMat *)                                      \
  _(Type, Subtype, BfMat *, ToType, BfMat const *, BfType)                     \
  _(Type, Subtype, BfMat *, Cholesky, BfMat const *)

#define INTERFACE BF_INTERFACE_Mat
BF_DEFINE_VTABLE_STRUCT(Mat)
#undef INTERFACE

struct BfMat {
  BfMatVtable *vtbl;
  BfMatProps props;
  BfSize numRows;
  BfSize numCols;
#if BF_DEBUG
  void *aux; /* pointer to extra user-defined data for purposes of debugging */
#endif
};

#define INTERFACE BF_INTERFACE_Mat
BF_DECLARE_INTERFACE(Mat)
#undef INTERFACE

void bfMatInit(BfMat *mat, BfMatVtable *vtbl, BfSize numRows, BfSize numCols);
void bfMatDeinit(BfMat *mat);
BfMat *bfMatFromFile(char const *path, BfSize numRows, BfSize numCols, BfDtype dtype);
bool bfMatInstanceOf(BfMat const *mat, BfType type);
bool bfMatIsTransposed(BfMat const *mat);
BfMat *bfMatTrans(BfMat *mat);
BfMat *bfMatConjTrans(BfMat *mat);
