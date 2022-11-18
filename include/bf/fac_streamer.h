#pragma once

#include "mat.h"
#include "tree.h"

typedef struct BfFacStreamer BfFacStreamer;

BfFacStreamer *bfFacStreamerNew();
void bfFacStreamerInit(BfFacStreamer *facStreamer,
                       BfTree const *rowTree, BfTree const *colTree,
                       BfSize rowTreeInitDepth, BfSize colTreeInitDepth,
                       BfReal tol);
void bfFacStreamerFeed(BfFacStreamer *facStreamer, BfMat const *mat);
bool bfFacStreamerDone(BfFacStreamer const *facStreamer);
BfMat *bfFacStreamerGetFac(BfFacStreamer const *facStreamer);
