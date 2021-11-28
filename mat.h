#pragma once

#include "def.h"
#include "dtype.h"
#include "vec.h"

enum BfMatProps {
  BF_MAT_PROP_NONE = 0,
  BF_MAT_PROP_DIAGONAL = (1 << 0)
};

typedef struct BfMat {
  enum BfDtypes dtype;
  enum BfMatProps props;
  BfSize shape[2];
  void *data;
} BfMat;

BfSize bfMatSize(BfMat const *A);
enum BfError bfMatNumBytes(BfMat const *A, BfSize *nbytes);

enum BfError bfFreeMat(BfMat *A);

enum BfError
bfInitEmptyMat(BfMat *A, enum BfDtypes dtype, enum BfMatProps props,
               BfSize const shape[2]);

enum BfError
bfInitEmptySvdMats(BfMat const *A, BfMat *U, BfMat *S, BfMat *Vt);

enum BfError
bfComputeMatSvd(BfMat const *A, BfMat *U, BfMat *S, BfMat *Vt);
