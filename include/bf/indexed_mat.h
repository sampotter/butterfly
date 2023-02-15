#pragma once

#include "types.h"

struct BfIndexedMat {
  /* The leading row index for `mat`. */
  BfSize i0;

  /* The leading column index for `mat`. */
  BfSize j0;

  /* The indexed matrix. */
  BfMat *mat;
};
