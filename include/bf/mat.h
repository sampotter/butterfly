#pragma once

#include "def.h"
#include "interface.h"
#include "mat_types.h"

// TODO: this is a little screwed up. The convention should probably be:
//
//   BfMat<Structure><Type>
//
// with Type = Block (or Blk... so: "BfMatDenseBlk") if we just keep a
// pointer to BfMat
typedef struct BfMat BfMat;
typedef struct BfMatBlock BfMatBlock;
typedef struct BfMatBlockCoo BfMatBlockCoo;
typedef struct BfMatBlockDense BfMatBlockDense;
typedef struct BfMatBlockDiag BfMatBlockDiag;
typedef struct BfMatDenseComplex BfMatDenseComplex;
typedef struct BfMatDiagReal BfMatDiagReal;
typedef struct BfMatProduct BfMatProduct;

typedef enum BfMatProps {
  BF_MAT_PROPS_NONE = 0,
  BF_MAT_PROPS_VIEW = 1 << 0,
  BF_MAT_PROPS_TRANS = 1 << 1,
  BF_MAT_PROPS_CONJ = 1 << 2,
  BF_MAT_PROPS_ORTHO = 1 << 3
} BfMatProps;

#define BF_INTERFACE_Mat(Type, Subtype, _)                                \
  _(Type, Subtype, void,      Delete,      BfMat **)                      \
  _(Type, Subtype, BfMat *,   EmptyLike,   BfMat const *, BfSize, BfSize) \
  _(Type, Subtype, BfMat *,   ZerosLike,   BfMat const *, BfSize, BfSize) \
  _(Type, Subtype, BfMatType, GetType,     BfMat const *)                 \
  _(Type, Subtype, bool,      InstanceOf,  BfMat const *, BfMatType)      \
  _(Type, Subtype, BfSize,    NumBytes,    BfMat const *)                 \
  _(Type, Subtype, void,      Save,        BfMat const *, char const *)   \
  _(Type, Subtype, BfSize,    GetNumRows,  BfMat const *)                 \
  _(Type, Subtype, BfSize,    GetNumCols,  BfMat const *)                 \
  _(Type, Subtype, BfMat *,   GetRowRange, BfMat *, BfSize, BfSize)       \
  _(Type, Subtype, BfMat *,   GetColRange, BfMat *, BfSize, BfSize)       \
  _(Type, Subtype, void,      SetRowRange, BfMat *, BfSize, BfSize, BfMat const *) \
  _(Type, Subtype, void,      AddInplace,  BfMat *, BfMat const *)        \
  _(Type, Subtype, BfMat *,   Mul,         BfMat const *, BfMat const *)  \
  _(Type, Subtype, BfMat *,   LstSq,       BfMat const *, BfMat const *)

#define INTERFACE BF_INTERFACE_Mat
BF_DEFINE_VTABLE_STRUCT(Mat);
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
BF_DECLARE_INTERFACE(Mat);
#undef INTERFACE

void bfMatInit(BfMat *mat, BfMatVtable *vtbl, BfSize numRows, BfSize numCols);
void bfMatDeinit(BfMat *mat);
bool bfMatIsTransposed(BfMat const *mat);
BfMat *bfMatConjTrans(BfMat *mat);
