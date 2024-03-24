#pragma once

#include <bf/fac_streamer.h>
#include <bf/points.h>

BfPoints1 *bfLboEigsToFreqs(BfVecReal const *Lam);

typedef struct {
  BfReal eigenbandTime;
  BfReal totalTime;
} BfLboFeedResult;
BfLboFeedResult bfLboFeedFacStreamerNextEigenband(BfFacStreamer *facStreamer, BfPoints1 *freqs, BfMat const *L, BfMat const *M);
