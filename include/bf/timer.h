#pragma once

#include <bf/def.h>
#include <bf/types.h>

#include <time.h>

struct BfTimer {
  clock_t clock;
};

void bfTimerReset(BfTimer *timer);
BfReal bfTimerGetElapsedTimeInSeconds(BfTimer *timer);
