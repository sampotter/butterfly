#pragma once

#include "def.h"

// TODO: rename: BF_MAT_TYPE_* -> BF_TYPE_MAT_* (so that first type is
// BF_TYPE_MAT, which is consistent with the rest of the typenames)
typedef enum BfMatTypes {
  BF_MAT_TYPE_MAT,
  BF_MAT_TYPE_BLOCK,
  BF_MAT_TYPE_BLOCK_COO,
  BF_MAT_TYPE_BLOCK_DENSE,
  BF_MAT_TYPE_BLOCK_DIAG,
  BF_MAT_TYPE_DENSE_COMPLEX,
  BF_MAT_TYPE_DIAG_REAL,
  BF_MAT_TYPE_PRODUCT,
  BF_MAT_TYPE_COUNT
} BfMatType;

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

typedef struct BfMatVtable {
  BfMat *(*emptyLike)(BfMat const *, BfSize, BfSize);
  BfMat *(*zerosLike)(BfMat const *, BfSize, BfSize);
  void (*deinit)(BfMat *);
  void (*delete)(BfMat **);
  void (*deinitAndDelete)(BfMat **);
  BfMatType (*getType)(BfMat const *);
  BfSize (*numBytes)(BfMat const *);
  void (*save)(BfMat const *, char const *);
  BfSize (*getNumRows)(BfMat const *);
  BfSize (*getNumCols)(BfMat const *);
  BfMat *(*getRowRange)(BfMat *, BfSize, BfSize);
  BfMat *(*getColRange)(BfMat *, BfSize, BfSize);
  void (*setRowRange)(BfMat *, BfSize, BfSize, BfMat const *);
  void (*addInplace)(BfMat *, BfMat const *);
  BfMat *(*mul)(BfMat const *, BfMat const *);
  BfMat *(*lstSq)(BfMat const *, BfMat const *);
} BfMatVtable;

#define BF_DEFINE_MAT_VTABLE(Subtype)                                   \
  static BfMatVtable matVtbl = {                                        \
    .emptyLike = (__typeof__(&bfMatEmptyLike))bf##Subtype##EmptyLike,   \
    .zerosLike = (__typeof__(&bfMatZerosLike))bf##Subtype##ZerosLike,   \
    .deinit = (__typeof__(&bfMatDeinit))bf##Subtype##Deinit,            \
    .delete = (__typeof__(&bfMatDelete))bf##Subtype##Delete,            \
    .deinitAndDelete = (__typeof__(&bfMatDeinitAndDelete))bf##Subtype##DeinitAndDelete, \
    .getType = (__typeof__(&bfMatGetType))bf##Subtype##GetType,         \
    .numBytes = (__typeof__(&bfMatNumBytes))bf##Subtype##NumBytes,      \
    .save = (__typeof__(&bfMatSave))bf##Subtype##Save,                  \
    .getNumRows = (__typeof__(&bfMatGetNumRows))bf##Subtype##GetNumRows, \
    .getNumCols = (__typeof__(&bfMatGetNumCols))bf##Subtype##GetNumCols, \
    .getRowRange = (__typeof__(&bfMatGetRowRange))bf##Subtype##GetRowRange, \
    .getColRange = (__typeof__(&bfMatGetColRange))bf##Subtype##GetColRange, \
    .setRowRange = (__typeof__(&bfMatSetRowRange))bf##Subtype##SetRowRange, \
    .addInplace = (__typeof__(&bfMatAddInplace))bf##Subtype##AddInplace, \
    .mul = (__typeof__(&bfMatMul))bf##Subtype##Mul,                     \
    .lstSq = (__typeof__(&bfMatLstSq))bf##Subtype##LstSq                \
  };

#define BF_DECLARE_INTERFACE_MAT(Subtype)                               \
  Bf##Subtype *bf##Subtype##EmptyLike(Bf##Subtype const *, BfSize, BfSize); \
  Bf##Subtype *bf##Subtype##ZerosLike(Bf##Subtype const *, BfSize, BfSize); \
  void bf##Subtype##Deinit(Bf##Subtype *);                              \
  void bf##Subtype##Delete(Bf##Subtype **);                             \
  void bf##Subtype##DeinitAndDelete(Bf##Subtype **);                    \
  BfMatType bf##Subtype##GetType(Bf##Subtype const *);                  \
  BfSize bf##Subtype##NumBytes(Bf##Subtype const *);                    \
  void bf##Subtype##Save(Bf##Subtype const *, char const *);            \
  BfSize bf##Subtype##GetNumRows(Bf##Subtype const *);                  \
  BfSize bf##Subtype##GetNumCols(Bf##Subtype const *);                  \
  Bf##Subtype *bf##Subtype##GetRowRange(Bf##Subtype *, BfSize, BfSize); \
  Bf##Subtype *bf##Subtype##GetColRange(Bf##Subtype *, BfSize, BfSize); \
  void bf##Subtype##SetRowRange(Bf##Subtype *, BfSize, BfSize, BfMat const *); \
  void bf##Subtype##AddInplace(Bf##Subtype *, BfMat const *);           \
  BfMat *bf##Subtype##Mul(Bf##Subtype const *, BfMat const *);          \
  BfMat *bf##Subtype##LstSq(Bf##Subtype const *, BfMat const *);

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
#if BF_DEBUG
  void *aux; /* pointer to extra user-defined data for purposes of debugging */
#endif
};

void bfMatInit(BfMat *mat, BfMatVtable *vtbl, BfSize numRows, BfSize numCols);
bool bfMatIsTransposed(BfMat const *mat);
BfMat *bfMatConjTrans(BfMat *mat);

/* BfMat interface: */
BfMat *bfMatEmptyLike(BfMat const *mat, BfSize numRows, BfSize numCols);
BfMat *bfMatZerosLike(BfMat const *mat, BfSize numRows, BfSize numCols);
void bfMatDeinit(BfMat *mat);
void bfMatDelete(BfMat **mat);
void bfMatDeinitAndDelete(BfMat **mat);
BfMatType bfMatGetType(BfMat const *mat);

BfSize bfMatNumBytes(BfMat const *mat);
void bfMatSave(BfMat const *mat, char const *path);
BfSize bfMatGetNumRows(BfMat const *mat);
BfSize bfMatGetNumCols(BfMat const *mat);
BfMat *bfMatGetRowRange(BfMat *mat, BfSize i0, BfSize i1);
BfMat *bfMatGetColRange(BfMat *mat, BfSize j0, BfSize j1);
void bfMatSetRowRange(BfMat *mat, BfSize i0, BfSize i1, BfMat const *rows);
void bfMatAddInplace(BfMat *lhs, BfMat const *rhs);
BfMat *bfMatMul(BfMat const *lhs, BfMat const *rhs);
BfMat *bfMatLstSq(BfMat const *lhs, BfMat const *rhs);
