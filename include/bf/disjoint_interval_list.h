#pragma once

#include "interval.h"

/** DisjointIntervalList: */

BfDisjointIntervalList *bfDisjointIntervalListNew();
void bfDisjointIntervalListInitEmpty(BfDisjointIntervalList *list);
void bfDisjointIntervalListDeinit(BfDisjointIntervalList *list);
void bfDisjointIntervalListDealloc(BfDisjointIntervalList **list);
void bfDisjointIntervalListDeinitAndDealloc(BfDisjointIntervalList **list);
BfSize bfDisjointIntervalListGetNumIntervals(BfDisjointIntervalList const *list);
BfInterval const *bfDisjointIntervalListGetPtrConst(BfDisjointIntervalList const *list, BfSize i);
void bfDisjointIntervalListAdd(BfDisjointIntervalList *list, BfInterval interval);
void bfDisjointIntervalListRemove(BfDisjointIntervalList *list, BfInterval interval);
