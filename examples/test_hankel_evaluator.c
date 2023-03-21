#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <bf/bessel.h>
#include <bf/const.h>
#include <bf/fac_streamer.h>
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

BfReal clenshaw_spoof(BfReal const *c, BfReal x) {
  BfReal x2 = 2*x;
  BfReal c0 = c[0];
  BfReal c1 = c[1];
  for (int i = 2; i < 10; ++i) {
    BfReal tmp = c1;
    c1 = c[i] - c0;
    c0 = tmp + c0*x2;
  }
  return c1 + c0*x;
}

////////////////////////////////////////////////////////////////////////////

static BfReal clenshaw(BfReal const *c, BfSize n, BfReal x) {
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

static void getChebyshevPoints(BfSize n, BfReal *x) {
  for (int i = 0; i < (int)n; ++i)
    x[i] = sin((BF_PI_OVER_TWO*(2*i + 1 - (int)n))/n);

//   for (BfSize i = 0; i < n; ++i)
//     x[i] = BF_NAN;
//   for (BfSize i = 0; i < n/2; ++i) {
//     BfReal t = i + 1;
//     t /= n + 1;
//     x[i] = -cos(BF_PI*t);
//     x[n - i - 1] = -x[i];
//   }
//   if (n % 2 == 1)
//     x[n/2] = 0;

//   BfReal w = BF_PI/(n - 1);
//   BfSize i = 0;
//   for (; i < (n + 1)/2; ++i)
//     x[i] = cos(w*i);
//   for (; i < n; ++i)
//     x[i] = -x[n - i - 1];
}

/* This function simultaneously constructs the transpose of the
 * pseudo-Vandermonde for the Chebyshev polynomials and projects onto
 * them, returning the Chebysehv coefficients for f transformed to the
 * domain [-1, 1]. */
static void getChebyshevCoefs(BfReal (*f)(BfReal), BfReal a, BfReal b, BfSize n, BfReal const *x, BfReal *c) {
  BfReal *y = malloc(n*sizeof(BfReal));
  for (BfSize j = 0; j < n; ++j)
    y[j] = f((b - a)*(x[j] + 1)/2 + a);

  BfReal *v0 = malloc(n*sizeof(BfReal));
  for (BfSize j = 0; j < n; ++j)
    v0[j] = 1;

  BfReal *v1 = malloc(n*sizeof(BfReal));
  for (BfSize j = 0; j < n; ++j)
    v1[j] = x[j];

  BfReal *v = malloc(n*sizeof(BfReal));

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

  free(v);
  free(v1);
  free(v0);
  free(y);
}

/* rough draft of tree interpolator follows */

typedef struct EvalTree EvalTree;
typedef struct EvalTreeNode EvalTreeNode;

// NOTE: a better way to do these is to allocate a pool of coefficient
// vectors and EvalTreeNodes in EvalTree... can use bitfields to
// combine "isLeaf" and the index into one. The idea here is that
// since every internal node has the same number of children, we only
// need the index to the starting position.

struct EvalTreeNode {
  BfReal a, b;
  bool isLeaf;
  BfReal *c;
  EvalTreeNode *children;
};

struct EvalTree {
  BfReal (*f)(BfReal);
  BfReal a, b;
  BfSize d;
  BfSize k;
  BfReal tol;
  BfReal *x; // chebyshev points

  EvalTreeNode *root;
};

void evalTreeNodeInitRec(EvalTreeNode *node, EvalTree const *tree) {
  BfSize n = tree->d + 1;

  node->c = malloc(n*sizeof(BfReal));
  getChebyshevCoefs(tree->f, node->a, node->b, n, tree->x, node->c);

  node->isLeaf = fabs(node->c[tree->d]) <= tree->tol;

  if (node->isLeaf) {
    node->children = NULL;
    return;
  }

  free(node->c);
  node->c = NULL;

  node->children = malloc(tree->k*sizeof(EvalTreeNode));

  BfReal delta = node->b - node->a;
  delta /= tree->k;
  for (BfSize i = 0; i < tree->k; ++i) {
    node->children[i].a = i == 0 ? node->a : node->a + i*delta;
    node->children[i].b = i == tree->k - 1 ? node->b : node->a + (i + 1)*delta;
    evalTreeNodeInitRec(&node->children[i], tree);
  }

  if (node->isLeaf) {
    assert(node->c != NULL);
    assert(node->children == NULL);
  } else {
    assert(node->c == NULL);
    assert(node->children != NULL);
  }
}

BfReal evalTreeNodeGetValueRec(EvalTreeNode const *node, EvalTree const *tree, BfReal x) {
  if (node->isLeaf) {
    assert(node->a <= x && x <= node->b);
    return clenshaw(node->c, tree->d + 1, 2*(x - node->a)/(node->b - node->a) - 1);
  } else {
    for (BfSize i = 0; i < tree->k; ++i) {
      EvalTreeNode const *child = &node->children[i];
      if (child->a <= x && x <= child->b)
        return evalTreeNodeGetValueRec(child, tree, x);
    }
  }
  assert(false);
}

static EvalTree *
makeEvalTree(BfReal (*f)(BfReal), BfReal a, BfReal b, BfSize d, BfSize k, BfReal tol) {
  EvalTree *tree = malloc(sizeof(EvalTree));
  tree->f = f;
  tree->a = a;
  tree->b = b;
  tree->d = d;
  tree->k = k;

  tree->tol = tol;

  tree->x = malloc((tree->d + 1)*sizeof(BfReal));
  getChebyshevPoints(tree->d + 1, tree->x);

  tree->root = malloc(sizeof(EvalTreeNode));
  tree->root->a = a;
  tree->root->b = b;
  evalTreeNodeInitRec(tree->root, tree);

  return tree;
}

BfReal evalTreeGetValue(EvalTree const *tree, BfReal x) {
  return evalTreeNodeGetValueRec(tree->root, tree, x);
}

////////////////////////////////////////////////////////////////////////////

int compar(void const *x, void const *y) {
  return (int)(*(BfReal const *)x - *(BfReal const *)y);
}

int main(int argc, char *argv[]) {
  Opts opts;
  int code = parseArgs(argc, argv, &opts);
  if (code) exit(EXIT_FAILURE);

  printf("approximation interval: [%g, %g]\n", opts.r0, opts.r1);
  printf("number of evaluation points: %lu\n", opts.numPoints);

  BfReal *X = malloc(opts.numPoints*sizeof(BfReal));
  if (!strcmp(opts.pointsType, "random")) {
    bfSeed(0);
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
    assert(false);
  }

  /* check current GSL timings */

  BfReal *J0_gsl = malloc(opts.numPoints*sizeof(BfReal));

  bfToc();
  for (BfSize i = 0; i < opts.numPoints; ++i)
    J0_gsl[i] = bf_j0(X[i]);
  BfReal time_gsl = bfToc();

  printf("GSL J0 eval time: %g s\n", time_gsl);
  printf("GSL J0 eval rate: %g pps\n", opts.numPoints/time_gsl);

  /* degree 10 clenshaw to check upper bound on rate */

  BfReal *tmp = malloc(opts.numPoints*sizeof(BfReal));

  BfReal const c[10] = {
    -1.26700706,  1.07870778, -0.32287527,  0.60970485, -1.33727439, -1.87921393,
    1.05588976, -0.91636183,  1.15975378,  0.82640578
  };

  bfToc();
  for (BfSize i = 0; i < opts.numPoints; ++i)
    tmp[i] = clenshaw_spoof(c, X[i]);
  BfReal time_clenshaw10 = bfToc();

  printf("degree 10 Clenshaw eval time: %g s\n", time_clenshaw10);
  printf("degree 10 Clenshaw eval rate: %g pps\n", opts.numPoints/time_clenshaw10);

  /* tree evaluator */

  bfToc();
  EvalTree *evalTree = makeEvalTree(bf_j0, opts.r0, opts.r1, opts.degree, opts.valence, opts.tol);
  BfReal time_makeEvalTree = bfToc();

  BfReal *J0_evalTree = malloc(opts.numPoints*sizeof(BfReal));

  bfToc();
  for (BfSize i = 0; i < opts.numPoints; ++i)
    J0_evalTree[i] = evalTreeGetValue(evalTree, X[i]);
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

  free(X);
  free(J0_gsl);
  free(tmp);
  free(J0_evalTree);
}
