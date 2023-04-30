#pragma once

#include "def.h"

void bfGetChebPts(BfSize n, BfReal *x);

/* Chebyshev series defined on the interval [-1, 1]. */
typedef struct {
  /* The order of the Chebyshev polynomial (equal to deg(p) + 1). */
  BfSize order;

  /* The Chebyshev coefficients (a length `order` array). */
  BfReal *c;
} BfChebStd;

BfChebStd *bfChebStdNewWithDegree(BfSize d);
void bfChebStdDeinit(BfChebStd *cheb);
void bfChebStdDealloc(BfChebStd **cheb);
void bfChebStdDelete(BfChebStd **cheb);
BfReal bfChebStdGetErrorEstimate(BfChebStd const *cheb);
void bfChebStdInterp(BfChebStd *cheb, BfReal (*f)(BfReal), BfReal a, BfReal b, BfReal const *x);
BfReal bfChebStdEval(BfChebStd const *cheb, BfReal x);

/* Chebyshev series defined on the interval [a, b]. */
typedef struct {
  /* Chebyshev series coefficients: */
  BfReal *c;

  /* Order of Chebyshev series: */
  BfSize order;

  BfReal a;
  BfReal b;
} BfCheb;

void bfChebInitWithDegree(BfCheb *cheb, BfSize d, BfSize a, BfSize b);
void bfChebDeinit(BfCheb *cheb);
void bfChebInterp(BfCheb *cheb, BfReal (*f)(BfReal), BfReal const *x);
BfReal bfChebEval(BfCheb const *cheb, BfReal x);
