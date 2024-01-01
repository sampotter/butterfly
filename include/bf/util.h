#pragma once

#include "def.h"

#include <stdio.h>

#include <bf/mat.h>

BfReal bfToc(void);

void bfRealArgsort(BfReal const *values, BfSize n, BfSize *perm);
void bfSizeSetConstant(BfSize numSizes, BfSize *size, BfSize value);
void bfSizeRunningSum(BfSize numSizes, BfSize *size);
bool bfSizeIsPerm(BfSize n, BfSize const *size);
void bfPrintBlocks(BfMat const *mat, BfSize level, FILE *fp);
BfSize bfGetFileSizeInBytes(char const *path);
void bfReadFileToMemory(char const *path, BfSize numBytes, BfByte *ptr);
void bfSort(BfPtr ptr, BfSize n, BfSize size, BfCompar compar, BfPtr aux);
bool bfFindZeroOnInterval(BfReal (*f)(BfReal, void *), BfReal a, BfReal b, void *context, BfReal *t);
int bfSignum(BfReal x);
