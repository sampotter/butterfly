#pragma once

#include "def.h"
#include "dtype.h"

typedef struct BfMat {
  enum BfDtypes dtype;
  BfSize ndim;
  BfSize *shape;
  void *data;
} BfMat;

BfSize bfMatSize(BfMat const *mat);
enum BfError bfMatNumBytes(BfMat const *mat, BfSize *num_bytes);

enum BfError bfFreeMat(BfMat *mat);

enum BfError
bfMakeEmptyMat(BfMat *mat,
               enum BfDtypes dtype, BfSize ndim, BfSize const *shape);

enum BfError
bfComputeMatSvd(BfMat const *mat, BfMat *U, BfMat *S, BfMat *Vt);
