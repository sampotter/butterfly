#include "cmocka_setup.h"

#include <bf/disjoint_interval_list.h>

static void testDisjointIntervalListAddAndDelete(void **state) {
  (void)state;

  BfDisjointIntervalList *list = bfDisjointIntervalListNewEmpty();

  bfDisjointIntervalListAdd(list, &(BfInterval) {.endpoint = {-0.5, 10.5}, .closed = {true, true}});

  for (BfSize i = 1; i <= 9; i += 2)
    bfDisjointIntervalListAdd(list, &(BfInterval) {.endpoint = {i, i + 1}, .closed = {true, true}});

  assert_int_equal(bfDisjointIntervalListGetNumIntervals(list), 1);

  bfDisjointIntervalListRemove(list, &(BfInterval) {.endpoint = {-1, 0}, .closed = {false, false}});
  bfDisjointIntervalListRemove(list, &(BfInterval) {.endpoint = {10, 11}, .closed = {false, false}});

  assert_int_equal(bfDisjointIntervalListGetNumIntervals(list), 1);
  {
    BfInterval const *interval = bfDisjointIntervalListGetPtrConst(list, 0);
    assert_true(interval->endpoint[0] == 0);
    assert_true(interval->endpoint[1] == 10);
    assert_true(interval->closed[0]);
    assert_true(interval->closed[1]);
  }

  for (BfSize i = 1; i <= 9; i += 2)
    bfDisjointIntervalListRemove(list, &(BfInterval) {.endpoint = {i, i + 1}, .closed = {true, true}});

  assert_int_equal(bfDisjointIntervalListGetNumIntervals(list), 5);
  {
    BfInterval const *interval = bfDisjointIntervalListGetPtrConst(list, 0);
    assert_true(interval->endpoint[0] == 0);
    assert_true(interval->endpoint[1] == 1);
    assert_true(interval->closed[0]);
    assert_true(!interval->closed[1]);
  }
  {
    BfInterval const *interval = bfDisjointIntervalListGetPtrConst(list, 1);
    assert_true(interval->endpoint[0] == 2);
    assert_true(interval->endpoint[1] == 3);
    assert_true(!interval->closed[0]);
    assert_true(!interval->closed[1]);
  }
  {
    BfInterval const *interval = bfDisjointIntervalListGetPtrConst(list, 2);
    assert_true(interval->endpoint[0] == 4);
    assert_true(interval->endpoint[1] == 5);
    assert_true(!interval->closed[0]);
    assert_true(!interval->closed[1]);
  }
  {
    BfInterval const *interval = bfDisjointIntervalListGetPtrConst(list, 3);
    assert_true(interval->endpoint[0] == 6);
    assert_true(interval->endpoint[1] == 7);
    assert_true(!interval->closed[0]);
    assert_true(!interval->closed[1]);
  }
  {
    BfInterval const *interval = bfDisjointIntervalListGetPtrConst(list, 4);
    assert_true(interval->endpoint[0] == 8);
    assert_true(interval->endpoint[1] == 9);
    assert_true(!interval->closed[0]);
    assert_true(!interval->closed[1]);
  }

  bfDisjointIntervalListDeinitAndDealloc(&list);
}

int main(void) {
  struct CMUnitTest const tests[] = {
    cmocka_unit_test(testDisjointIntervalListAddAndDelete)
  };

  return cmocka_run_group_tests(tests, NULL, NULL);
}
