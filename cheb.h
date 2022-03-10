#pragma once

#include "def.h"

/* Chebyshev series defined on the interval [a, b]. */
typedef struct {
  /* Chebyshev series coefficients: */
  BfReal *c;

  /* Order of Chebyshev series: */
  BfSize order;

  /* TODO: define a BfChebStd where [a, b] is assumed to equal [-1, 1] */
  BfReal a_plus_b;
  BfReal b_minus_a;
} BfCheb;

BfReal bfChebEval(BfCheb const *cheb, BfReal x);
