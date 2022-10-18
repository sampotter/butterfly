#pragma once

#include <stdio.h>

#include "def.h"
#include "dtype.h"
#include "interface.h"
#include "perm.h"
#include "types.h"

typedef enum BfVecProps {
  BF_VEC_PROPS_NONE = 0,
  BF_VEC_PROPS_VIEW = (1 << 0)
} BfVecProps;

#define BF_INTERFACE_Vec(Type, Subtype, _)                               \
  _(Type, Subtype, BfVec *, Copy, BfVec const *)                         \
  _(Type, Subtype, void, Delete, BfVec **)                               \
  _(Type, Subtype, BfType, GetType, BfVec const *)                       \
  _(Type, Subtype, bool, InstanceOf, BfVec const *, BfType)              \
  _(Type, Subtype, BfPtr, GetEltPtr, BfVec *, BfSize)                    \
  _(Type, Subtype, BfVec *, GetSubvecCopy, BfVec const *, BfSize, BfSize) \
  _(Type, Subtype, BfVec *, GetSubvecView, BfVec *, BfSize, BfSize)      \
  _(Type, Subtype, void, Print, BfVec const *, FILE *)                   \
  _(Type, Subtype, BfReal, Dist, BfVec const *, BfVec const *)           \
  _(Type, Subtype, BfReal, NormMax, BfVec const *)                       \
  _(Type, Subtype, void, ScaleByReal, BfVec *, BfReal)                   \
  _(Type, Subtype, void, AddInplace, BfVec *, BfVec const *)             \
  _(Type, Subtype, void, MulInplace, BfVec *, BfMat const *)             \
  _(Type, Subtype, void, SolveInplace, BfVec *, BfMat const *)           \
  _(Type, Subtype, void, RecipInplace, BfVec *)                          \
  _(Type, Subtype, BfMat *, GetGivensRotation, BfVec const *, BfSize, BfSize) \
  _(Type, Subtype, void, Permute, BfVec *, BfPerm const *)               \
  _(Type, Subtype, BfVec *, Concat, BfVec const *, BfVec const *)

#define INTERFACE BF_INTERFACE_Vec
BF_DEFINE_VTABLE_STRUCT(Vec)
#undef INTERFACE

struct BfVec {
  BfVecVtable *vtbl;
  BfVecProps props;
  BfSize size;
};

#define INTERFACE BF_INTERFACE_Vec
BF_DECLARE_INTERFACE(Vec)
#undef INTERFACE

void bfVecInit(BfVec *vec, BfVecVtable *vtbl, BfSize size);
void bfVecDeinit(BfVec *vec);
BfVec *bfVecFromFile(char const *path, BfSize size, BfDtype dtype);
bool bfVecInstanceOf(BfVec const *vec, BfType type);
