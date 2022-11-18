#pragma once

#include "def.h"

typedef struct BfSpmat BfSpmat;

/* Upcasting: */

BfMat *bfSpmatToMat(BfSpmat *spmat);
