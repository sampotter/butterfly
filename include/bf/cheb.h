#pragma once

#include "def.h"

/* Chebyshev series defined on the interval [-1, 1]. */
typedef struct {
  BfReal *c;
  BfSize order;
} BfChebStd;

BfReal bfChebStdEval(BfChebStd const *cheb, BfReal x);

/* Chebyshev series defined on the interval [a, b]. */
typedef struct {
  /* Chebyshev series coefficients: */
  BfReal *c;

  /* Order of Chebyshev series: */
  BfSize order;

  BfReal a_plus_b;
  BfReal b_minus_a;
} BfCheb;

BfReal bfChebEval(BfCheb const *cheb, BfReal x);
