#define _DEFAULT_SOURCE 1

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <bf/assert.h>
#include <bf/bessel.h>
#include <bf/cheb.h>
#include <bf/const.h>
#include <bf/eval_tree.h>
#include <bf/fac_streamer.h>
#include <bf/mem.h>
#include <bf/rand.h>
#include <bf/util.h>

#include "argtable3.h"

int const MAX_NUM_ARG_ERRORS = 20;

typedef struct {
  BfReal r0, r1;
  BfSize numPoints;
  BfSize degree;
  BfSize valence;
  BfReal tol;
  char const *pointsType;
} Opts;

static int parseArgs(int argc, char *argv[], Opts *opts) {
  int code = 0, nerrors = 0;
  char const *progname = argv[0];

  struct arg_lit *help;
  struct arg_lit *verbose;
  struct arg_dbl *r0, *r1;
  struct arg_int *numPoints;
  struct arg_str *pointsType;
  struct arg_int *degree;
  struct arg_int *valence;
  struct arg_dbl *tol;
  struct arg_end *end;

  void *argtable[] = {
    help = arg_litn(NULL, "help", 0, 1, "display this help and exit"),
    verbose = arg_litn("v", "verbose", 0, 1, "verbose output"),
    r0 = arg_dbln(NULL, "r0", "<r0>", 0, 1, "the start of the approximation interval"),
    r1 = arg_dbln(NULL, "r1", "<r1>", 0, 1, "the end of the approximation interval"),
    numPoints = arg_intn("n", "num_points", "<n>", 0, 1, "number of random eval points"),
    pointsType = arg_strn(NULL, "points_type", "<points_type>", 0, 1, "type of point distribution (can be \"random\" or \"uniform\")"),
    degree = arg_intn("d", "degree", "<d>", 0, 1, "degree of Chebyshev poly used for tree evaluator"),
    valence = arg_intn("k", "valence", "<k>", 0, 1, "the valence of the tree evaluator"),
    tol = arg_dbln(NULL, "tol", "<tol>", 0, 1, "tolerance used to build tree evaluator"),
    end = arg_end(MAX_NUM_ARG_ERRORS),
  };
  nerrors = arg_parse(argc, argv, argtable);

  /* special case: '--help' takes precedence over error reporting */
  if (help->count > 0) {
    printf("Usage: %s", progname);
    arg_print_syntax(stdout, argtable, "\n");
    printf("Demonstrate command-line parsing in argtable3.\n\n");
    arg_print_glossary(stdout, argtable, "  %-25s %s\n");
    code = 0;
    goto exit;
  }

  /* If the parser returned any errors then display them and exit */
  if (nerrors > 0) {
    /* Display the error details contained in the arg_end struct.*/
    arg_print_errors(stdout, end, progname);
    printf("Try '%s --help' for more information.\n", progname);
    code = 1;
    goto exit;
  }

  opts->r0 = *r0->dval;
  opts->r1 = *r1->dval;
  opts->numPoints = *numPoints->ival;
  opts->pointsType = *pointsType->sval;
  opts->degree = *degree->ival;
  opts->valence = *valence->ival;
  opts->tol = *tol->dval;

exit:
  arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

  return code;
}

int compar(void const *x, void const *y) {
  return (int)(*(BfReal const *)x - *(BfReal const *)y);
}

int main(int argc, char *argv[]) {
  Opts opts;
  int code = parseArgs(argc, argv, &opts);
  if (code) exit(EXIT_FAILURE);

  printf("approximation interval: [%g, %g]\n", opts.r0, opts.r1);
  printf("number of evaluation points: %lu\n", opts.numPoints);

  bfSeed(0);

  BfReal *X = bfMemAlloc(opts.numPoints, sizeof(BfReal));
  if (!strcmp(opts.pointsType, "random")) {
    bfRealUniform(opts.numPoints, X);
    for (BfSize i = 0; i < opts.numPoints; ++i)
      X[i] = (opts.r1 - opts.r0)*X[i] + opts.r0;
    qsort(X, opts.numPoints, sizeof(BfReal), compar);
  } else if (!strcmp(opts.pointsType, "uniform")) {
    for (BfSize i = 0; i < opts.numPoints; ++i) {
      BfReal t = i;
      t /= opts.numPoints - 1;
      X[i] = (1 - t)*opts.r0 + t*opts.r1;
      X[i] = fmax(opts.r0, fmin(X[i], opts.r1));
    }
  } else {
    BF_ASSERT(false);
  }

  /* check sqrt timings */

  BfReal *sqrtX = bfMemAlloc(opts.numPoints, sizeof(BfReal));

  bfToc();
  for (BfSize i = 0; i < opts.numPoints; ++i)
    sqrtX[i] = sqrt(X[i]);
  BfReal time_sqrt = bfToc();

  printf("C std lib sqrt eval time: %g s\n", time_sqrt);
  printf("C std lib sqrt eval rate: %g pps\n", opts.numPoints/time_sqrt);

  /* check current C std library timings for j0 */

  BfReal *J0_std = bfMemAlloc(opts.numPoints, sizeof(BfReal));

  bfToc();
  for (BfSize i = 0; i < opts.numPoints; ++i)
    J0_std[i] = j0(X[i]);
  BfReal time_std = bfToc();

  printf("C std lib J0 eval time: %g s\n", time_std);
  printf("C std lib J0 eval rate: %g pps\n", opts.numPoints/time_std);

  /* check current GSL timings */

  BfReal *J0_gsl = bfMemAlloc(opts.numPoints, sizeof(BfReal));

  bfToc();
  for (BfSize i = 0; i < opts.numPoints; ++i)
    J0_gsl[i] = bf_j0(X[i]);
  BfReal time_gsl = bfToc();

  printf("GSL J0 eval time: %g s\n", time_gsl);
  printf("GSL J0 eval rate: %g pps\n", opts.numPoints/time_gsl);

  /* degree 10 clenshaw to check upper bound on rate */

  BfReal *tmp = bfMemAlloc(opts.numPoints, sizeof(BfReal));

  BfChebStd randCheb = {
    .c = bfMemAlloc((opts.degree + 1), sizeof(BfReal)),
    .order = opts.degree + 1
  };
  bfRealRandn(randCheb.order, randCheb.c);

  bfToc();
  for (BfSize i = 0; i < opts.numPoints; ++i) {
    BfReal t = 2*(X[i] - opts.r0)/(opts.r1 - opts.r0) - 1;
    tmp[i] = bfChebStdEval(&randCheb, t);
  }
  BfReal time_clenshaw = bfToc();

  printf("degree %lu Clenshaw eval time: %g s\n", opts.degree, time_clenshaw);
  printf("degree %lu Clenshaw eval rate: %g pps\n", opts.degree, opts.numPoints/time_clenshaw);

  /* tree evaluator */

  bfToc();
  BfEvalTree *evalTree = bfEvalTreeNew();
  BfEvalTreeSpec spec = {
    .f = j0,
    .a = opts.r0,
    .b = opts.r1,
    .d = opts.degree,
    .k = opts.valence,
    .tol = opts.tol
  };
  bfEvalTreeInit(evalTree, &spec);
  BfReal time_makeEvalTree = bfToc();

  BfReal *J0_evalTree = bfMemAlloc(opts.numPoints, sizeof(BfReal));

  bfToc();
  for (BfSize i = 0; i < opts.numPoints; ++i)
    J0_evalTree[i] = bfEvalTreeGetValue(evalTree, X[i]);
  BfReal time_evalTree = bfToc();

  printf("time to build eval tree: %g s\n", time_makeEvalTree);
  printf("EvalTree eval time: %g s\n", time_evalTree);
  printf("EvalTree eval rate: %g pps\n", opts.numPoints/time_evalTree);

  BfReal maxAbsError = 0;
  for (BfSize i = 0; i < opts.numPoints; ++i) {
    BfReal absError = fabs(J0_gsl[i] - J0_evalTree[i]);
    maxAbsError = fmax(absError, maxAbsError);
  }
  printf("max abs error between GSL and eval tree: %g\n", maxAbsError);

  /* clean up */

  bfMemFree(X);
  bfMemFree(J0_gsl);
  bfMemFree(tmp);
  bfMemFree(J0_evalTree);
}
