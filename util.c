#include "util.h"

void bfSizeSetConstant(BfSize numSizes, BfSize *size, BfSize value) {
  for (BfSize i = 0; i < numSizes; ++i)
    size[i] = value;
}

void bfSizeRunningSum(BfSize numSizes, BfSize *size) {
  for (BfSize i = 1; i < numSizes; ++i)
    size[i] += size[i - 1];
}
