#pragma once

#include "def.h"
#include "geom.h"

void bfSeed(BfSize seed);

BfSize bfSizeUniform1(BfSize low, BfSize high);

BfReal bfRealUniform1();
void bfRealUniform(BfSize n, BfReal *x);

void bfRealRandn(BfSize n, BfReal *x);

void bfComplexRandn(BfSize n, BfComplex *x);

void bfSampleRandomUnitVector2(BfVector2 u);
