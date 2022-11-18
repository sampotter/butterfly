#pragma once

#include "mat.h"

BfReal bfGetEigMax(BfMat const *L, BfMat const *M);
void bfGetEigBand(BfMat const *L, BfMat const *M, BfReal lam0, BfReal lam1,
                  BfMat **Phi, BfMat **Lam);
