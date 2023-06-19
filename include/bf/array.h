#pragma once

#include "def.h"
#include "types.h"

BfArray *bfArrayNewUninitialized(void);
BfArray *bfArrayNewEmpty(BfSize eltSize);
void bfArrayInitEmpty(BfArray *array, BfSize eltSize);
void bfArrayDeinit(BfArray *array);
void bfArrayDealloc(BfArray **array);
void bfArrayDeinitAndDealloc(BfArray **array);
