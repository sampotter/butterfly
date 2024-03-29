#pragma once

enum BfError {
  BF_ERROR_NONE,
  BF_ERROR_INVALID_ARGUMENTS,
  BF_ERROR_RUNTIME_ERROR,
  BF_ERROR_NOT_IMPLEMENTED,
  BF_ERROR_MEMORY_ERROR,
  BF_ERROR_OUT_OF_RANGE,
  BF_ERROR_FILE_ERROR,
  BF_ERROR_TYPE_ERROR,
  BF_ERROR_INCOMPATIBLE_SHAPES,
#ifdef BF_EMBREE
  BF_ERROR_EMBREE
#endif
};


enum BfError bfGetError(void);
void bfSetError(enum BfError error);
