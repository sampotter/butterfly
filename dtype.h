#pragma once

#include "def.h"
#include "error.h"

enum BfDtypes {
  BF_DTYPE_REAL,
  BF_DTYPE_COMPLEX
};

enum BfError bfSizeOfDtype(enum BfDtypes dtype, BfSize *size);

bool bfDtypeIsValid(enum BfDtypes dtype);
