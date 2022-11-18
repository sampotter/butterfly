#pragma once

#include "tree.h"

typedef enum BfStreamModes { BF_STREAM_MODE_POST_ORDER } BfStreamMode;

typedef struct BfFacStreamer BfFacStreamer;

BfFacStreamer *bfFacStreamerNew();
void bfFacStreamerInit(BfFacStreamer *facStreamer, BfTree const *rowTree,
                       BfTree const *colTree, BfStreamingMode mode,
                       BfSize rowTreeInitDepth, BfSize colTreeInitDepth,
                       BfReal tol);
void bfFacStreamerFeed(BfFacStreamer *facStreamer, BfMat const *mat);
bool bfFacStreamerDone(BfFacStreamer const *facStreamer);