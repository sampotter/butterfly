#include "cheb.h"

BfReal bfChebStdEval(BfChebStd const *cheb, BfReal x) {
  BfReal d  = 0.0;
  BfReal dd = 0.0;

  for (int j = cheb->order; j >= 1; j--) {
    BfReal temp = d;
    d = 2*x*d - dd + cheb->c[j];
    dd = temp;
  }

  d = x*d - dd + 0.5 * cheb->c[0];

  return d;
}

BfReal bfChebEval(BfCheb const *cheb, BfReal x) {
  BfReal d  = 0.0;
  BfReal dd = 0.0;

  BfReal y  = (2.0*x - cheb->a_plus_b)/cheb->b_minus_a;
  BfReal y2 = 2.0 * y;

  for (int j = cheb->order; j >= 1; j--) {
    BfReal temp = d;
    d = y2*d - dd + cheb->c[j];
    dd = temp;
  }

  d = y*d - dd + 0.5 * cheb->c[0];

  return d;
}
