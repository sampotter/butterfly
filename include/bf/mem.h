#pragma once

#include "def.h"

BfPtr bfMemAlloc(BfSize n, BfSize size);
BfPtr bfMemAllocAndZero(BfSize n, BfSize size);
BfPtr bfMemAllocCopy(BfConstPtr src, BfSize n, BfSize size);
void bfMemFree(BfPtr ptr);
BfPtr bfMemRealloc(BfPtr ptr, BfSize n, BfSize size);
void bfMemCopy(BfConstPtr src, BfSize n, BfSize size, BfPtr dst);
void bfMemMove(BfPtr src, BfSize n, BfSize size, BfPtr dst);
void bfMemZero(BfPtr src, BfSize n, BfSize size);
