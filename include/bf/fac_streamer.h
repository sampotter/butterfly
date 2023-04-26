#pragma once

#include "fac.h"
#include "mat.h"
#include "tree.h"

typedef struct BfFacStreamer BfFacStreamer;

BfFacStreamer *bfFacStreamerNew();
void bfFacStreamerInit(BfFacStreamer *facStreamer, BfFacSpec const *facSpec);
void bfFacStreamerDeinit(BfFacStreamer *facStreamer);
BfSize bfFacStreamerGetNumRows(BfFacStreamer const *facStreamer);
void bfFacStreamerFeed(BfFacStreamer *facStreamer, BfMat const *Phi);
bool bfFacStreamerIsDone(BfFacStreamer const *facStreamer);
BfFac *bfFacStreamerGetFac(BfFacStreamer const *facStreamer);
BfTreeNode *bfFacStreamerGetCurrentColumnNode(BfFacStreamer const *facStreamer);
BfPerm const *bfFacStreamerGetRowTreeReversePerm(BfFacStreamer const *facStreamer);
