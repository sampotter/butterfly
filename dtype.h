#pragma once

#include "def.h"
#include "error.h"

enum BfDtypes {
  BF_DTYPE_VOID,
  BF_DTYPE_REAL,
  BF_DTYPE_COMPLEX,
  BF_DTYPE_MAT
};

BfSize bfDtypeSize(enum BfDtypes dtype);

bool bfDtypeIsValid(enum BfDtypes dtype);

void bfGetDtypeZero(enum BfDtypes dtype, BfPtr ptr);
