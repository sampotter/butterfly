#pragma once

#include "geom.h"
#include "types.h"

BfMatCsrReal *bfRadiosityGetViewFactorMatrix(BfTrimesh const *trimesh, BfSizeArray const *rowInds, BfSizeArray const *colInds);
