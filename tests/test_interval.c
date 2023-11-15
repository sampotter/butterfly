#include "cmocka_setup.h"

#include <bf/interval.h>

static void test_interval_difference(void **state) {
  (void)state;

  /* Expected behavior: [-1, 2] \ (0, 1) = [-1, 0] U [1, 2]. */

  BfInterval I1 = {.endpoint = {-1, 2}, .closed = {true, true}};
  BfInterval I2 = {.endpoint = {0, 1}, .closed = {false, false}};

  BfInterval diff[2];
  BfSize n = bfIntervalDifference(&I1, &I2, diff);

  assert_int_equal(n, 2);
  assert_true(diff[0].endpoint[0] == -1);
  assert_true(diff[0].endpoint[1] == 0);
  assert_true(diff[0].closed[0]);
  assert_true(diff[0].closed[1]);
  assert_true(diff[1].endpoint[0] == 1);
  assert_true(diff[1].endpoint[1] == 2);
  assert_true(diff[1].closed[0]);
  assert_true(diff[1].closed[1]);
}

int main(void) {
  struct CMUnitTest const tests[] = {
    cmocka_unit_test(test_interval_difference)
  };

  return cmocka_run_group_tests(tests, NULL, NULL);
}
