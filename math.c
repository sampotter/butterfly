#include "math.h"

#include <gsl/gsl_sf_bessel.h>

BfReal bf_j0(BfReal x) {
  return gsl_sf_bessel_J0(x);
}

BfReal bf_y0(BfReal x) {
  return gsl_sf_bessel_Y0(x);
}
