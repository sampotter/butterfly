#pragma once

#include "interval.h"

/** DisjointIntervalList: */

BfDisjointIntervalList *bfDisjointIntervalListNewEmpty();
void bfDisjointIntervalListInitEmpty(BfDisjointIntervalList *list);
void bfDisjointIntervalListDeinit(BfDisjointIntervalList *list);
void bfDisjointIntervalListDealloc(BfDisjointIntervalList **list);
void bfDisjointIntervalListDeinitAndDealloc(BfDisjointIntervalList **list);
BfSize bfDisjointIntervalListGetNumIntervals(BfDisjointIntervalList const *list);
BfInterval const *bfDisjointIntervalListGetPtrConst(BfDisjointIntervalList const *list, BfSize i);
void bfDisjointIntervalListAdd(BfDisjointIntervalList *list, BfInterval const *interval);
void bfDisjointIntervalListRemove(BfDisjointIntervalList *list, BfInterval const *interval);
bool bfDisjointIntervalListIsEmpty(BfDisjointIntervalList const *list);
BfInterval const *bfDisjointIntervalListGetFirstPtrConst(BfDisjointIntervalList const *list);
