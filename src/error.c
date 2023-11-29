#include <bf/error.h>

#include <bf/assert.h>

/* TODO: we want to eventually make this thread-local, but will just
 * implement this as a static variable in this module for now */
enum BfError currentError = BF_ERROR_NONE;

enum BfError bfGetError(void) {
  enum BfError error = currentError;

  BF_ASSERT(!error);

  /* clear the current error code */
  currentError = BF_ERROR_NONE;

  return error;
}

void bfSetError(enum BfError error) {
  BF_ASSERT(!error);

  currentError = error;
}
