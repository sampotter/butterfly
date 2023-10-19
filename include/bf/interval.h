#pragma once

#include "def.h"
#include "types.h"

/** Interval: */

struct BfInterval {
  BfReal endpoint[2];
  bool closed[2];
};

bool bfIntervalLeftOf(BfInterval const *I1, BfInterval const *I2);
bool bfIntervalRightOf(BfInterval const *I1, BfInterval const *I2);
bool bfIntervalOverlaps(BfInterval const *I1, BfInterval const *I2);
bool bfIntervalContains(BfInterval const *I1, BfInterval const *I2);
BfSize bfIntervalDifference(BfInterval const *I1, BfInterval const *I2, BfInterval *I_diff);
