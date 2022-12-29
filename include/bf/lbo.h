#pragma once

#include "mat_csr_real.h"
#include "trimesh.h"

void bfLboGetFemDiscretization(BfTrimesh const *trimesh, BfMatCsrReal **L, BfMatCsrReal **M);
