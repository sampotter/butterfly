#include <bf/array.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>

struct BfArray {
  BfSize eltSize;
  BfSize capacity;
  BfSize size;
  BfPtr *data;
  bool isView;
};

static void invalidate(BfArray *array) {
  array->eltSize = BF_SIZE_BAD_VALUE;
  array->capacity = BF_SIZE_BAD_VALUE;
  array->size = BF_SIZE_BAD_VALUE;
  array->data = NULL;
  array->isView = false;
}

BfArray *bfArrayNewUninitialized(void) {
  BF_ERROR_BEGIN();

  BfArray *array = bfMemAlloc(1, sizeof(BfArray));
  HANDLE_ERROR();

  invalidate(array);

  BF_ERROR_END() {
    BF_DIE();
  }

  return array;
}

BfArray *bfArrayNewEmpty(BfSize eltSize) {
  BF_ERROR_BEGIN();

  BfArray *array = bfArrayNewUninitialized();
  HANDLE_ERROR();

  bfArrayInitEmpty(array, eltSize);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return array;
}

void bfArrayInitEmpty(BfArray *array, BfSize eltSize) {
  BF_ERROR_BEGIN();

  array->eltSize = eltSize;
  array->capacity = BF_ARRAY_DEFAULT_CAPACITY;
  array->size = 0;

  array->data = bfMemAlloc(array->capacity, eltSize);
  HANDLE_ERROR();

  array->isView = false;

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfArrayDeinit(BfArray *array) {
  bfMemFree(array->data);

  invalidate(array);
}

void bfArrayDealloc(BfArray **array) {
  bfMemFree(*array);
  *array = NULL;
}

void bfArrayDeinitAndDealloc(BfArray **array) {
  bfArrayDeinit(*array);
  bfArrayDealloc(array);
}
