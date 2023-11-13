#pragma once

#include "def.h"
#include "types.h"

/** Interval: */

struct BfInterval {
  BfReal endpoint[2];
  bool closed[2];
};

void bfIntervalDeinit(BfInterval *interval);
bool bfIntervalLeftOf(BfInterval const *I1, BfInterval const *I2);
bool bfIntervalRightOf(BfInterval const *I1, BfInterval const *I2);
bool bfIntervalOverlaps(BfInterval const *I1, BfInterval const *I2);
bool bfIntervalContainsInterval(BfInterval const *I1, BfInterval const *I2);
bool bfIntervalContainsPoint(BfInterval const *I1, BfReal x);
BfSize bfIntervalDifference(BfInterval const *I1, BfInterval const *I2, BfInterval *I_diff);
BfReal bfIntervalGetMidpoint(BfInterval const* interval);
bool bfIntervalIsEmpty(BfInterval const *interval);
bool bfIntervalEquals(BfInterval const *I1, BfInterval const *I2);
