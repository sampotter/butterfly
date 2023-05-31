#include <bf/cheb.h>

#include <math.h>

#include <bf/assert.h>
#include <bf/const.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>

/** Utilities: */

void bfGetChebPts(BfSize n, BfReal *x) {
  for (int i = 0; i < (int)n; ++i)
    x[i] = sin((BF_PI_OVER_TWO*(2*i + 1 - (int)n))/n);
}

/** Implementation: ChebStd */

BfChebStd *bfChebStdNewWithDegree(BfSize d) {
  BfChebStd *cheb = bfMemAlloc(1, sizeof(BfChebStd));
  cheb->order = d + 1;
  cheb->c = bfMemAlloc(cheb->order, sizeof(BfReal));
  return cheb;
}

void bfChebStdDeinit(BfChebStd *cheb) {
  cheb->order = BF_SIZE_BAD_VALUE;
  bfMemFree(cheb->c);
  cheb->c = NULL;
}

void bfChebStdDealloc(BfChebStd **cheb) {
  bfMemFree(*cheb);
  *cheb = NULL;
}

void bfChebStdDelete(BfChebStd **cheb) {
  bfChebStdDeinit(*cheb);
  bfChebStdDealloc(cheb);
}

BfReal bfChebStdGetErrorEstimate(BfChebStd const *cheb) {
  return fabs(cheb->c[cheb->order - 1]);
}

/* This function simultaneously constructs the transpose of the
 * pseudo-Vandermonde for the Chebyshev polynomials and projects onto
 * them, returning the Chebysehv coefficients for f transformed to the
 * domain [-1, 1]. */
void bfChebStdInterp(BfChebStd *cheb, BfReal (*f)(BfReal), BfReal a, BfReal b, BfReal const *x) {
  BF_ERROR_BEGIN();

  BfSize n = cheb->order;
  BfReal *c = cheb->c;

  bool shouldFreeX = false;
  if (x == NULL) {
    x = bfMemAlloc(n, sizeof(BfReal));
    HANDLE_ERROR();

    bfGetChebPts(n, (BfReal *)x);

    shouldFreeX = true;
  }

  BfReal *y = bfMemAlloc(n, sizeof(BfReal));
  for (BfSize j = 0; j < n; ++j)
    y[j] = f((b - a)*(x[j] + 1)/2 + a);

  BfReal *v0 = bfMemAlloc(n, sizeof(BfReal));
  for (BfSize j = 0; j < n; ++j)
    v0[j] = 1;

  BfReal *v1 = bfMemAlloc(n, sizeof(BfReal));
  for (BfSize j = 0; j < n; ++j)
    v1[j] = x[j];

  BfReal *v = bfMemAlloc(n, sizeof(BfReal));

  c[0] = 0;
  for (BfSize j = 0; j < n; ++j)
    c[0] += y[j];

  c[1] = 0;
  for (BfSize j = 0; j < n; ++j)
    c[1] += x[j]*y[j];

  for (BfSize i = 2; i < n; ++i) {
    c[i] = 0;
    for (BfSize j = 0; j < n; ++j) {
      v[j] = 2*v1[j]*x[j] - v0[j];
      c[i] += v[j]*y[j];
      v0[j] = v1[j];
      v1[j] = v[j];
    }
  }

  c[0] /= n;
  for (BfSize i = 1; i < n; ++i) c[i] *= 2.0/n;

  BF_ERROR_END() {
    BF_DIE();
  }

  if (shouldFreeX)
    bfMemFree((BfReal *)x);

  bfMemFree(v);
  bfMemFree(v1);
  bfMemFree(v0);
  bfMemFree(y);
}

BfReal bfChebStdEval(BfChebStd const *cheb, BfReal x) {
  BfSize n = cheb->order;
  BfReal const *c = cheb->c;

  BfReal c0, c1, tmp;
  if (n == 1) {
    c0 = c[0];
    c1 = 0;
  } else if (n == 2) {
    c0 = c[0];
    c1 = c[1];
  } else {
    c0 = c[n - 2];
    c1 = c[n - 1];
    for (BfSize i = 3; i <= n; ++i) {
      tmp = c0;
      c0 = c[n - i] - c1;
      c1 = tmp + 2*c1*x;
    }
  }
  return c0 + c1*x;
}

void bfChebInitWithDegree(BfCheb *cheb, BfSize d, BfSize a, BfSize b) {
  BF_ERROR_BEGIN();

  cheb->order = d + 1;

  cheb->c = bfMemAlloc(cheb->order, sizeof(BfReal));
  HANDLE_ERROR();

  cheb->a = a;
  cheb->b = b;

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfChebDeinit(BfCheb *cheb) {
  cheb->order = BF_SIZE_BAD_VALUE;

  bfMemFree(cheb->c);
  cheb->c = NULL;

  cheb->a = BF_NAN;
  cheb->b = BF_NAN;
}

void bfChebInterp(BfCheb *cheb, BfReal (*f)(BfReal), BfReal const *x) {
  BF_ERROR_BEGIN();

  BfChebStd chebStd = {.order = cheb->order, .c = cheb->c};

  bfChebStdInterp(&chebStd, f, cheb->a, cheb->b, x);
  HANDLE_ERROR();

  BF_ERROR_END() {}
}

BfReal bfChebEval(BfCheb const *cheb, BfReal x) {
  BfChebStd chebStd = {.order = cheb->order, .c = cheb->c};
  BfReal xStd = (2*x - cheb->a - cheb->b)/(cheb->b - cheb->a);
  return bfChebStdEval(&chebStd, xStd);
}
