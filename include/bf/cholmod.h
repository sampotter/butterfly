#pragma once

#include <suitesparse/cholmod.h>

#include <stdbool.h>

void bfCholmodStart();
void bfCholmodFinish();
bool bfCholmodIsStarted();
cholmod_common *bfCholmodGetCommon();
