#pragma once

#include "def.h"

void bfSeed(BfSize seed);

void bfRealUniform(BfSize n, BfReal *x);
void bfRealRandn(BfSize n, BfReal *x);

void bfComplexRandn(BfSize n, BfComplex *x);
