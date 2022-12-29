#include <bf/cholmod.h>

#include <bf/error.h>
#include <bf/error_macros.h>

static bool started = false;

static cholmod_common Common;

void bfCholmodStart() {
  if (started)
    bfSetError(BF_ERROR_RUNTIME_ERROR);
  else
    cholmod_start(&Common);
}

void bfCholmodFinish() {
  if (started)
    cholmod_finish(&Common);
  else
    bfSetError(BF_ERROR_RUNTIME_ERROR);
}

bool bfCholmodIsStarted() {
  return started;
}

cholmod_common *bfCholmodGetCommon() {
  return (void *)&Common;
}
