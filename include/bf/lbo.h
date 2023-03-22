#pragma once

#include "mat_csr_real.h"
#include "trimesh.h"

void bfLboGetFemDiscretization(BfTrimesh const *trimesh, BfMat **L, BfMat **M);
