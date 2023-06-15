#include <bf/mem.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>

BfPtr bfMemAlloc(BfSize n, BfSize size) {
  BF_ERROR_BEGIN();

  // TODO: check if n*size overflows

  // Some things to try:
  // - [ ] assert that size is a power of two
  // - [ ] do aligned allocation
  // - [ ] collect allocation statistics
  // - [ ] better debugging...

  BfPtr ptr = malloc(n*size);
  if (ptr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BF_ERROR_END() {
    BF_DIE();
  }

  return ptr;
}

BfPtr bfMemAllocAndZero(BfSize n, BfSize size) {
  BF_ERROR_BEGIN();

  BfPtr ptr = bfMemAlloc(n, size);
  if (ptr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  bfMemZero(ptr, n, size);

  BF_ERROR_END() {
    BF_DIE();
  }

  return ptr;
}

BfPtr bfMemAllocCopy(BfConstPtr src, BfSize n, BfSize size) {
  BF_ERROR_BEGIN();

  BfPtr copy = bfMemAlloc(n, size);
  HANDLE_ERROR();

  bfMemCopy(src, n, size, copy);

  BF_ERROR_END() {
    BF_DIE();
  }

  return copy;
}

void bfMemFree(BfPtr ptr) {
  free(ptr);
}

BfPtr bfMemRealloc(BfPtr ptr, BfSize n, BfSize size) {
  BF_ERROR_BEGIN();

  BfPtr newPtr = realloc(ptr, n*size);
  if (newPtr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BF_ERROR_END() {
    BF_DIE();
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
