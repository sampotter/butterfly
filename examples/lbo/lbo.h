#pragma once

#include <bf/fac_streamer.h>
#include <bf/points.h>

BfPoints1 *convertEigsToFreqs(BfVecReal const *Lam);
void feedFacStreamerNextEigenband(BfFacStreamer *facStreamer, BfPoints1 *freqs, BfMat const *L, BfMat const *M);
