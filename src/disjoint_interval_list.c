#include <bf/disjoint_interval_list.h>

#include <bf/array.h>
#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/interval.h>
#include <bf/mem.h>

#include <math.h>

struct BfDisjointIntervalList {
  BfArray *intervals;
};

BfDisjointIntervalList *bfDisjointIntervalListNewEmpty() {
  BF_ERROR_BEGIN();

  BfDisjointIntervalList *list = bfMemAlloc(1, sizeof(BfDisjointIntervalList));
  HANDLE_ERROR();

  bfDisjointIntervalListInitEmpty(list);

  BF_ERROR_END() {
    BF_DIE();
  }

  return list;
}

void bfDisjointIntervalListInitEmpty(BfDisjointIntervalList *list) {
  BF_ERROR_BEGIN();

  list->intervals = bfArrayNewEmpty(sizeof(BfInterval));
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfDisjointIntervalListDeinit(BfDisjointIntervalList *list) {
  bfArrayDeinitAndDealloc(&list->intervals);
}

void bfDisjointIntervalListDealloc(BfDisjointIntervalList **list) {
  bfMemFree(*list);
  *list = NULL;
}

void bfDisjointIntervalListDeinitAndDealloc(BfDisjointIntervalList **list) {
  bfDisjointIntervalListDeinit(*list);
  bfDisjointIntervalListDealloc(list);
}

BfSize bfDisjointIntervalListGetNumIntervals(BfDisjointIntervalList const *list) {
  return bfArrayGetSize(list->intervals);
}

BfInterval const *bfDisjointIntervalListGetPtrConst(BfDisjointIntervalList const *list, BfSize i) {
  return bfArrayGetPtr(list->intervals, i);
}

static BfInterval merge(BfInterval const *I1, BfInterval const *I2) {
  BF_ASSERT(bfIntervalOverlaps(I1, I2));

  BfInterval I_merged;

  I_merged.closed[0] = I_merged.closed[1] = false;

  I_merged.endpoint[0] = fmin(I1->endpoint[0], I2->endpoint[0]);
  I_merged.endpoint[1] = fmax(I1->endpoint[1], I2->endpoint[1]);

  I_merged.closed[0] = I_merged.endpoint[0] == I1->endpoint[0] && !I1->endpoint[0];
  I_merged.closed[0] = I_merged.endpoint[0] == I2->endpoint[0] && !I2->endpoint[0];

  I_merged.closed[1] = I_merged.endpoint[1] == I1->endpoint[1] && !I1->endpoint[1];
  I_merged.closed[1] = I_merged.endpoint[1] == I2->endpoint[1] && !I2->endpoint[1];

  return I_merged;
}

void bfDisjointIntervalListAdd(BfDisjointIntervalList *list, BfInterval const *interval) {
  BF_ERROR_BEGIN();

  BfInterval mergeInterval = *interval;

  BfSize i = bfArrayGetSize(list->intervals);
  BfInterval *otherInterval = NULL;
  while (i > 0) {
    otherInterval = bfArrayGetPtr(list->intervals, --i);

    if (bfIntervalLeftOf(&mergeInterval, otherInterval))
      break;

    if (bfIntervalOverlaps(&mergeInterval, otherInterval)) {
      bfArrayRemove(list->intervals, i);
      mergeInterval = merge(&mergeInterval, otherInterval);
    }
  };

  bfArrayInsert(list->intervals, i, &mergeInterval);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfDisjointIntervalListRemove(BfDisjointIntervalList *list, BfInterval const *interval) {
  BF_ERROR_BEGIN();

  BfSize i = 0;
  while (i < bfArrayGetSize(list->intervals)) {
    BfInterval const *otherInterval = bfArrayGetPtr(list->intervals, i);

    if (bfIntervalEquals(interval, otherInterval)) {
      bfArrayRemove(list->intervals, i++);
      continue;
    }

    BfInterval newIntervals[2];
    BfSize numNewIntervals = bfIntervalDifference(otherInterval, interval, newIntervals);
    if (numNewIntervals > 0) {
      bfArrayRemove(list->intervals, i);
      for (BfSize j = 0; j < numNewIntervals; ++j) {
        bfArrayInsert(list->intervals, i++, &newIntervals[j]);
        HANDLE_ERROR();
      }
    } else {
      ++i;
    }
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

bool bfDisjointIntervalListIsEmpty(BfDisjointIntervalList const *list) {
  return bfArrayIsEmpty(list->intervals);
}

BfInterval const *bfDisjointIntervalListGetFirstPtrConst(BfDisjointIntervalList const *list) {
  return bfArrayGetFirstPtrConst(list->intervals);
}
