#pragma once

#include <stdio.h>

#include <bf/mat.h>

#include "def.h"

BfReal bfToc();

void bfSizeSetConstant(BfSize numSizes, BfSize *size, BfSize value);
void bfSizeRunningSum(BfSize numSizes, BfSize *size);
void bfPrintBlocks(BfMat const *mat, BfSize level, FILE *fp);
