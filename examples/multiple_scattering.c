#include <bf/bbox.h>
#include <bf/const.h>
#include <bf/ellipse.h>
#include <bf/helm2.h>
#include <bf/linalg.h>
#include <bf/mat_diag_real.h>
#include <bf/poisson_disk_sampling.h>
#include <bf/quadtree.h>
#include <bf/rand.h>
#include <bf/util.h>
#include <bf/vectors.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "argtable3.h"

int const MAX_NUM_ARG_ERRORS = 20;

typedef struct {
  BfReal minDist;
  BfReal axisLow;
  BfReal axisHigh;
  BfReal h;
  BfReal wavenumber;
} Opts;

struct K_helm2_wkspc {
  BfPoints2 const *points;
  BfVectors2 const *normals;
  BfReal K;
  BfReal nu;
};

BfComplex K_helm2(BfSize i, BfSize j, void *aux) {
  struct K_helm2_wkspc *wkspc = aux;
  BfReal const *xsrc = &wkspc->points->data[i][0];
  BfReal const *xtgt = &wkspc->points->data[j][0];
  BfReal const *nsrc = &wkspc->normals->data[i][0];
  BfReal G = bfHelm2GetKernelValue(
    xsrc, xtgt, NULL, NULL, wkspc->K, BF_LAYER_POTENTIAL_SINGLE);
  BfReal dGdN = bfHelm2GetKernelValue(
    xsrc, xtgt, nsrc, NULL, wkspc->K, BF_LAYER_POTENTIAL_PV_DOUBLE);
  return G + 1i*wkspc->nu*dGdN;
}

static bool parseArgs(int argc, char *argv[], Opts *opts) {
  bool success = true;

  struct arg_lit *help;
  struct arg_dbl *minDist;
  struct arg_dbl *axisLow;
  struct arg_dbl *axisHigh;
  struct arg_dbl *h;
  struct arg_dbl *wavenumber;
  struct arg_end *end;

  void *argtable[] = {
    help = arg_litn(NULL, "help", 0, 1, "Display help and exit"),
    minDist = arg_dbln("r", "minDist", "<r>", 0, 1, "Minimum distance between ellipse centers"),
    axisLow = arg_dbln("a", "axisLow", "<a>", 0, 1, "Lower bound for uniformly sampled ellipse axes"),
    axisHigh = arg_dbln("b", "axisHigh", "<b>", 0, 1, "Upper bound for uniformly sampled ellipse axes"),
    h = arg_dbln("h", NULL, "<h>", 0, 1, "The mesh fineness"),
    wavenumber = arg_dbln("k", "wavenumber", "<k>", 0, 1, "The wavenumber for the problem"),
    end = arg_end(MAX_NUM_ARG_ERRORS)
  };

  BfSize numErrors = arg_parse(argc, argv, argtable);

  if (help->count > 0) {
    printf("Usage: %s", argv[0]);
    arg_print_syntax(stdout, argtable, "\n");
    printf("Run a test problem demonstrating multiple scattering\n");
    arg_print_glossary(stdout, argtable, "  %-25s %s\n");
    goto cleanup;
  }

  if (numErrors > 0) {
    arg_print_errors(stdout, end, argv[0]);
    printf("Try '%s --help' for more information.\n", argv[0]);
    success = false;
    goto cleanup;
  }

  opts->minDist = *minDist->dval;
  opts->axisLow = *axisLow->dval;
  opts->axisHigh = *axisHigh->dval;
  opts->h = *h->dval;
  opts->wavenumber = *wavenumber->dval;

cleanup:
  arg_freetable(argtable, sizeof(argtable)/sizeof(argtable[0]));

  return success;
}

int main(int argc, char *argv[]) {
  Opts opts;
  bool success = parseArgs(argc, argv, &opts);
  if (!success) exit(EXIT_SUCCESS);

  bfSeed(0);

  puts("solving multiple scattering problem:");

  BfReal k = opts.wavenumber;

  BfReal nu = k;

  BfVector2 d;
  bfSampleRandomUnitVector2(d);

  BfReal tol = 1e-12;

  BfSize nxEval = 257;
  BfSize nyEval = 257;

  /** Set up problem geometry: randomly sample lots of little
   ** well-separated ellipses with different semimajor/minor axes and
   ** orientations: */

  bfToc();

  /* Set up bounding box of sampling domain---give a little margin to
   * keep ellipses fully inside */
  BfReal R = 1 - opts.minDist;
  BfBbox2 ellipseBbox = {.min = {-R, -R}, .max = {R, R}};
  BfBbox2 bbox = {.min = {-1, -1}, .max = {1, 1}};

  /* Sample ellipse centers */
  BfPoints2 *ellipseCenters = bfPoints2SamplePoissonDisk(&ellipseBbox, opts.minDist, 30);
  bfSavePoints2(ellipseCenters, "ellipseCenters.bin");

  BfSize numEllipses = bfPoints2GetSize(ellipseCenters);

  /* Randomly sample ellipses */
  BfEllipse *ellipse = malloc(numEllipses*sizeof(BfEllipse));
  for (BfSize i = 0; i < numEllipses; ++i) {
    /* Sample ellipse semi-major and semi-minor axes: */
    BfReal a = (opts.axisHigh - opts.axisLow)*bfRealUniform1() + opts.axisLow;
    BfReal b = (opts.axisHigh - opts.axisLow)*bfRealUniform1() + opts.axisLow;
    ellipse[i].semiMajorAxis = fmax(a, b);
    ellipse[i].semiMinorAxis = fmin(a, b);

    /* Give ellipse random orientation */
    ellipse[i].theta = BF_TWO_PI*bfRealUniform1();

    /* Set ellipse center */
    memcpy(ellipse[i].center, ellipseCenters->data[i], sizeof(BfPoint2));
  }

  printf("- created %lu random ellipses [%0.2fs]\n", numEllipses, bfToc());

  /** Sample discretization points from the ellipses using inverse
   ** curvature weighting: */

  bfToc();

  BfPoints2 *X = bfPoints2NewEmpty();
  BfVectors2 *N = bfVectors2NewEmpty();
  for (BfSize i = 0; i < numEllipses; ++i) {
    BfReal p = bfEllipseGetPerimeter(&ellipse[i]);
    BfSize n = floor(p/opts.h) + 1;
    bfEllipseSampleEquispaced(&ellipse[i], n, X, NULL, N);
    // bfEllipseSampleWithInverseCurvatureSpacing(&ellipse[i], n, X, NULL, N);
  }
  bfSavePoints2(X, "X.bin");
  bfSaveVectors2(N, "N.bin");

  BfSize n = bfPoints2GetSize(X);
  printf("- sampled %lu discretization points [%0.2fs]\n", n, bfToc());

  /** Build a quadtree on the discretization points: */

  bfToc();

  BfQuadtree quadtree;
  bfQuadtreeInit(&quadtree, X, N);

  printf("- set up quadtree [%0.2fs]\n", bfToc());

  BfTree *tree = bfQuadtreeToTree(&quadtree);
  BfPerm const *perm = bfTreeGetPermConst(tree);
  BfPerm revPerm = bfPermGetReversePerm(perm);

  /** Set up the dense kernel matrix using the combined field integral
   * equation for the exterior Dirichlet problem: */

  bfToc();

  BfMat *oneHalfEye = bfMatDiagRealToMat(bfMatDiagRealNewConstant(n, n, 0.5));
  BfMat *K = bfGetHelm2KernelMatrix(X, X, NULL, NULL, k, BF_LAYER_POTENTIAL_SINGLE);
  BfMat *D = bfGetHelm2KernelMatrix(X, X, N, NULL, k, BF_LAYER_POTENTIAL_PV_DOUBLE);
  bfMatAddInplace(K, oneHalfEye);
  bfMatScale(D, 1j*nu);
  bfMatAddInplace(K, D);

//   /* Apply KR correction */
//   struct K_helm2_wkspc K_wkspc = {.points = X, .normals = N, .K = k};
//   bf_apply_KR_correction(K, 6, K_helm2, (void *)&K_wkspc);

  printf("- set up dense kernel matrix [%0.2fs]\n", bfToc());

  /** Set up the RHS for the scattering problem: */

  bfToc();

  BfMat *rhs = NULL;
  { BfMatDenseComplex *_ = bfMatDenseComplexNew();
    bfMatDenseComplexInit(_, n, 1);
    for (BfSize i = 0; i < n; ++i) {
      BfPoint2 x;
      bfPoints2Get(X, i, x);
      *(_->data + i*_->rowStride) = cexp(1i*k*(d[0]*x[0] + d[1]*x[1]));
    }
    rhs = bfMatDenseComplexToMat(_); }

  printf("- set up RHS for scattering problem [%0.2fs]\n", bfToc());

  /** Solve the dense system using GMRES: */

  bfToc();

  BfSize numIterDense;
  BfMat *sigmaDense = bfSolveGMRES(K, rhs, NULL, tol, n, &numIterDense);

  printf("- solve dense system using GMRES in %lu iterations [%0.2fs]\n", numIterDense, bfToc());

  /** Set up evaluation grid: */

  bfToc();

  BfPoints2 *XEval = bfPoints2NewGrid(&bbox, nxEval, nyEval);

  printf("- set up %lu x %lu evaluation grid [%0.2fs]\n", nxEval, nyEval, bfToc());

  bfSavePoints2(XEval, "XEval.bin");

  /** Compute uIn on evaluation grid: */

  bfToc();

  BfMat *uIn = NULL;
  { BfMatDenseComplex *_ = bfMatDenseComplexNew();
    bfMatDenseComplexInit(_, nxEval*nyEval, 1);
    for (BfSize i = 0; i < nxEval*nyEval; ++i) {
      BfPoint2 x;
      bfPoints2Get(XEval, i, x);
      *(_->data + i*_->rowStride) = cexp(1i*k*(d[0]*x[0] + d[1]*x[1]));
    }
    uIn = bfMatDenseComplexToMat(_); }

  printf("- computed uIn on evaluation grid [%0.2fs]\n", bfToc());

  bfMatSave(uIn, "uIn.bin");

  /** Evaluate solution to dense system: */

  bfToc();

  BfMat *KEval = bfGetHelm2KernelMatrix(X, XEval, NULL, NULL, k, BF_LAYER_POTENTIAL_SINGLE);
  BfMat *DEval = bfGetHelm2KernelMatrix(X, XEval, N, NULL, k, BF_LAYER_POTENTIAL_PV_DOUBLE);
  bfMatScale(DEval, 1j*nu);
  bfMatAddInplace(KEval, DEval);

  BfMat *uScatDense = bfMatMul(KEval, sigmaDense);

  printf("- computed scattered field for dense problem [%0.2fs]\n", bfToc());

  bfMatSave(uScatDense, "uScatDense.bin");

  (void)revPerm;
  (void)sigmaDense;
}
