#pragma once

#include <suitesparse/cholmod.h>

#include <stdbool.h>

void bfCholmodStart(void);
void bfCholmodFinish(void);
bool bfCholmodIsStarted(void);
cholmod_common *bfCholmodGetCommon(void);
