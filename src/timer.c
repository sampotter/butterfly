#include <bf/timer.h>

void bfTimerReset(BfTimer *timer) {
  timer->clock = clock();
}

BfReal bfTimerGetElapsedTimeInSeconds(BfTimer *timer) {
  clock_t prevClock = timer->clock;
  timer->clock = clock();
  return ((double)timer->clock - (double)prevClock)/CLOCKS_PER_SEC;
}
