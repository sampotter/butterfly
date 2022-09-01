#include <bf/error.h>

#include <assert.h>

/* TODO: we want to eventually make this thread-local, but will just
 * implement this as a static variable in this module for now */
enum BfError currentError = BF_ERROR_NONE;

enum BfError bfGetError() {
  enum BfError error = currentError;

  assert(!error);

  /* clear the current error code */
  currentError = BF_ERROR_NONE;

  return error;
}

void bfSetError(enum BfError error) {
  assert(!error);

  currentError = error;
}
