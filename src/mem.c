#include <bf/mem.h>

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <bf/error.h>
#include <bf/error_macros.h>

BfPtr bfMemAlloc(BfSize n, BfSize size) {
  BEGIN_ERROR_HANDLING();

  // TODO: check if n*size overflows

  BfPtr ptr = malloc(n*size);
  if (ptr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {
    assert(false);
  }

  return ptr;
}

BfPtr bfMemAllocAndZero(BfSize n, BfSize size) {
  BEGIN_ERROR_HANDLING();

  // TODO: check if n*size overflows

  BfPtr ptr = calloc(n, size);
  if (ptr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {
    assert(false);
  }

  return ptr;
}

void bfMemFree(BfPtr ptr) {
  free(ptr);
}

BfPtr bfMemRealloc(BfPtr ptr, BfSize n, BfSize size) {
  BEGIN_ERROR_HANDLING();

  BfPtr newPtr = realloc(ptr, n*size);
  if (newPtr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {
    assert(false);
  }

  return newPtr;
}

void bfMemCopy(BfConstPtr src, BfSize n, BfSize size, BfPtr dst) {
  // TODO: check if n*size overflows
  memcpy(dst, src, n*size);
}

void bfMemMove(BfPtr src, BfSize n, BfSize size, BfPtr dst) {
  // TODO: check if n*size overflows
  memmove(dst, src, n*size);
}

void bfMemZero(BfPtr src, BfSize n, BfSize size) {
  // TODO: check if n*size overflows
  memset(src, 0x0, n*size);
}
