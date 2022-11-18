#pragma once

typedef struct BfIntervalTree BfIntervalTree;

/* Upcasting: */
BfTree *bfIntervalTreeToTree(BfIntervalTree *intervalTree);

BfIntervalTree *bfIntervalTreeNew();
void bfIntervalTreeInit(BfIntervalTree *intervalTree, BfReal a, BfReal b,
                        BfSize initDepth, BfSize k);
