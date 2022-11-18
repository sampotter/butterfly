#pragma once

#include "spmat.h"
#include "trimesh.h"

void bfLboGetFemDiscretization(BfTrimesh const *trimesh, BfSpmat **L, BfSpmat **M);
