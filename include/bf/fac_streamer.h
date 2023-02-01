#pragma once

#include "mat.h"
#include "tree.h"

typedef struct BfFacStreamer BfFacStreamer;

BfFacStreamer *bfFacStreamerNew();
void bfFacStreamerInit(BfFacStreamer *facStreamer, BfTree *rowTree, BfTree *colTree, BfSize rowTreeInitDepth, BfSize colTreeInitDepth, BfReal tol);
void bfFacStreamerDeinit(BfFacStreamer *facStreamer);
void bfFacStreamerFeed(BfFacStreamer *facStreamer, BfMat const *Phi, BfVecReal const *Lam);
bool bfFacStreamerIsDone(BfFacStreamer const *facStreamer);
BfMat *bfFacStreamerGetFac(BfFacStreamer const *facStreamer);
BfTreeNode *bfFacStreamerGetCurrentColumnNode(BfFacStreamer const *facStreamer);
