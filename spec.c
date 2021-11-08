#include "spec.h"

enum BfError bfSpecValidate(BfSpec const *spec) {
  if (spec->domain != BF_DOMAIN_R2) {
    return BF_ERROR_INVALID_ARGUMENTS;
  }
}
