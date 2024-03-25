#pragma once

#include <bf/def.h>
#ifdef BF_SINGLE
#  error "BF_DOUBLE required for Cholmod"
#endif

#include <suitesparse/cholmod.h>

#include <stdbool.h>

void bfCholmodStart(void);
void bfCholmodFinish(void);
bool bfCholmodIsStarted(void);
cholmod_common *bfCholmodGetCommon(void);
