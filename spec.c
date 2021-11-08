#include "spec.h"

FlyError flySpecValidate(FlySpec const *spec) {
  if (spec->dim != FLY_DIM_2D) {
    return FLY_ERROR_BAD_ARGUMENTS;
  }
}
