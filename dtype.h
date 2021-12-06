#pragma once

#include "def.h"
#include "error.h"

enum BfDtypes {
  BF_DTYPE_REAL,
  BF_DTYPE_COMPLEX
};

BfSize bfDtypeSize(enum BfDtypes dtype);

enum BfError bfSizeOfDtype(enum BfDtypes dtype, BfSize *size);

bool bfDtypeIsValid(enum BfDtypes dtype);
