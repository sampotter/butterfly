#pragma once

#include "fac.h"
#include "mat.h"
#include "tree.h"

typedef struct BfFacStreamer BfFacStreamer;

BfFacStreamer *bfFacStreamerNew(void);
void bfFacStreamerInit(BfFacStreamer *facStreamer, BfFacSpec const *facSpec);
void bfFacStreamerDeinit(BfFacStreamer *facStreamer);
void bfFacStreamerDealloc(BfFacStreamer **facStreamer);
void bfFacStreamerDelete(BfFacStreamer **facStreamer);
BfSize bfFacStreamerGetNumRows(BfFacStreamer const *facStreamer);
void bfFacStreamerFeed(BfFacStreamer *facStreamer, BfMat *Phi);
bool bfFacStreamerIsDone(BfFacStreamer const *facStreamer);
BfFac *bfFacStreamerGetFac(BfFacStreamer const *facStreamer);
BfFacSpan *bfFacStreamerGetFacSpan(BfFacStreamer const *facStreamer);
BfTreeNode *bfFacStreamerGetCurrentColumnNode(BfFacStreamer const *facStreamer);
BfPerm const *bfFacStreamerGetRowTreeReversePerm(BfFacStreamer const *facStreamer);
