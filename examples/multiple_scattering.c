#include <bf/bbox.h>
#include <bf/const.h>
#include <bf/ellipse.h>
#include <bf/helm2.h>
#include <bf/linalg.h>
#include <bf/mat_diag_real.h>
#include <bf/poisson_disk_sampling.h>
#include <bf/quadtree.h>
#include <bf/rand.h>
#include <bf/real_array.h>
#include <bf/size_array.h>
#include <bf/util.h>
#include <bf/vectors.h>

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "argtable3.h"

#include <fmm2d/helmholtz.h>

int const MAX_NUM_ARG_ERRORS = 20;

typedef struct {
  /* Problem parameters: */
  BfReal wavenumber;
  BfReal minDist;
  BfReal axisLow;
  BfReal axisHigh;

  /* Discretization parameters: */
  BfReal h;

  /* Postprocessing parameters: */
  BfBbox2 bboxEval;
  BfSize nxEval;
  BfSize nyEval;
} Opts;

struct K_helm2_wkspc {
  BfPoints2 const *points;
  BfVectors2 const *normals;
  BfReal k;
};

BfComplex K_helm2(BfSize i, BfSize j, void *aux) {
  struct K_helm2_wkspc *wkspc = aux;

  BfReal const *xsrc = &wkspc->points->data[i][0];
  BfReal const *xtgt = &wkspc->points->data[j][0];
  BfReal const *nsrc = &wkspc->normals->data[i][0];
  BfReal k = wkspc->k;

  BfComplex G = bfHelm2GetKernelValue(
    xsrc, xtgt, NULL, NULL, k, BF_LAYER_POTENTIAL_SINGLE);
  BfComplex dGdN = bfHelm2GetKernelValue(
    xsrc, xtgt, nsrc, NULL, k, BF_LAYER_POTENTIAL_PV_DOUBLE);
  BfComplex K = dGdN - 1i*wkspc->k*G;

  return K;
}

static bool parseArgs(int argc, char *argv[], Opts *opts) {
  bool success = true;

  /** All CLI parameters: */

  struct arg_lit *help;

  /* Problem parameters: */
  struct arg_dbl *wavenumber;
  struct arg_dbl *minDist;
  struct arg_dbl *axisLow;
  struct arg_dbl *axisHigh;

  /* Discretization parameters: */
  struct arg_dbl *h;

  /* Postprocessing parameters: */
  struct arg_dbl *xminEval;
  struct arg_dbl *xmaxEval;
  struct arg_dbl *yminEval;
  struct arg_dbl *ymaxEval;
  struct arg_int *nxEval;
  struct arg_int *nyEval;

  struct arg_end *end;

  /** Set up argtable and parse: */

  void *argtable[] = {
    help = arg_litn(NULL, "help", 0, 1, "Display help and exit"),

    /* Problem parameters: */
    wavenumber = arg_dbln("k", "wavenumber", "<k>", 0, 1, "The wavenumber for the problem"),
    minDist = arg_dbln("r", "minDist", "<r>", 0, 1, "Minimum distance between ellipse centers"),
    axisLow = arg_dbln("a", "axisLow", "<a>", 0, 1, "Lower bound for uniformly sampled ellipse axes"),
    axisHigh = arg_dbln("b", "axisHigh", "<b>", 0, 1, "Upper bound for uniformly sampled ellipse axes"),

    /* Discretization parameters: */
    h = arg_dbln("h", NULL, "<h>", 0, 1, "The mesh fineness"),

    /* Postprocessing parameters: */
    xminEval = arg_dbln(NULL, "xmin", "<xmin>", 0, 1, "The minimum x coord. of the eval. box"),
    xmaxEval = arg_dbln(NULL, "xmax", "<xmax>", 0, 1, "The maximum x coord. of the eval. box"),
    yminEval = arg_dbln(NULL, "ymin", "<ymin>", 0, 1, "The minimum y coord. of the eval. box"),
    ymaxEval = arg_dbln(NULL, "ymax", "<ymax>", 0, 1, "The maximum y coord. of the eval. box"),
    nxEval = arg_intn(NULL, "nx", "<nx>", 0, 1, "The number of eval. box nodes in x direction"),
    nyEval = arg_intn(NULL, "ny", "<ny>", 0, 1, "The number of eval. box nodes in y direction"),

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

  /** Extract parameters from parsed CLI options: */

  /* Problem parameters: */
  opts->wavenumber = *wavenumber->dval;
  opts->minDist = *minDist->dval;
  opts->axisLow = *axisLow->dval;
  opts->axisHigh = *axisHigh->dval;

  /* Discretization parameters: */
  opts->h = *h->dval;

  /* Postprocessing parameters: */
  opts->bboxEval.min[0] = *xminEval->dval;
  opts->bboxEval.max[0] = *xmaxEval->dval;
  opts->bboxEval.min[1] = *yminEval->dval;
  opts->bboxEval.max[1] = *ymaxEval->dval;
  opts->nxEval = *nxEval->ival;
  opts->nyEval = *nyEval->ival;

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

  BfVector2 d;
  bfSampleRandomUnitVector2(d);

  printf("- uIn = exp(i*k*d*r), where:\n");
  printf("  * k = %g\n", k);
  printf("  * d = (%g, %g)\n", d[0], d[1]);

  BfReal tol = 1e-12;

  /** Set up problem geometry: randomly sample lots of little
   ** well-separated ellipses with different semimajor/minor axes and
   ** orientations: */

  bfToc();

  /* Set up bounding box of sampling domain---give a little margin to
   * keep ellipses fully inside */
  BfReal R = 1 - opts.minDist;
  BfBbox2 ellipseBbox = {.min = {-R, -R}, .max = {R, R}};

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

  /* Keep track of the offset to each ellipse---we'll use these
   * indices to figure out where to apply the KR correction and do
   * preconditioning *blockwise*. */
  BfSizeArray *ellipseOffsets = bfSizeArrayNewWithDefaultCapacity();
  bfSizeArrayAppend(ellipseOffsets, 0);

  BfPoints2 *X = bfPoints2NewEmpty();
  BfVectors2 *N = bfVectors2NewEmpty();
  BfRealArray *W = bfRealArrayNewWithDefaultCapacity();
  for (BfSize i = 0; i < numEllipses; ++i) {
    BfReal p = bfEllipseGetPerimeter(&ellipse[i]);
    BfSize n = floor(p/opts.h) + 1;
    bfEllipseSampleLinspaced(&ellipse[i], n, X, NULL, N, W);
    // bfEllipseSampleEquispaced(&ellipse[i], n, X, NULL, N);
    // bfEllipseSampleWithInverseCurvatureSpacing(&ellipse[i], n, X, NULL, N);

    /* Add the current offset: */
    bfSizeArrayAppend(ellipseOffsets, bfPoints2GetSize(X));
  }
  bfSavePoints2(X, "X.bin");
  bfSaveVectors2(N, "N.bin");

  BfSize n = bfPoints2GetSize(X);
  printf("- sampled %lu discretization points [%0.2fs]\n", n, bfToc());

  /** Build a quadtree on the discretization points: */

  /** Set up the dense kernel matrix using the combined field integral
   * equation for the exterior Dirichlet problem: */

  bfToc();

  BfMat *oneHalfEye = bfMatDiagRealToMat(bfMatDiagRealNewConstant(n, n, 0.5));
  BfMat *K = bfGetHelm2KernelMatrix(X, X, NULL, NULL, k, BF_LAYER_POTENTIAL_SINGLE);
  BfMat *D = bfGetHelm2KernelMatrix(X, X, N, NULL, k, BF_LAYER_POTENTIAL_PV_DOUBLE);
  bfMatScale(K, -1j*k);
  bfMatAddInplace(K, D);

  /* Apply KR correction blockwise: */
  for (BfSize i = 0; i < numEllipses; ++i) {
    BfSize i0 = bfSizeArrayGet(ellipseOffsets, i);
    BfSize i1 = bfSizeArrayGet(ellipseOffsets, i + 1);

    struct K_helm2_wkspc wkspc;
    wkspc.points = bfPoints2GetRangeView(X, i0, i1);
    wkspc.normals = bfVectors2GetRangeView(N, i0, i1);
    wkspc.k = k;

    BfMat *KBlock = bfMatGetBlockView(K, i0, i1, i0, i1);
    bf_apply_KR_correction(KBlock, 6, K_helm2, (void *)&wkspc);
  }

  BfVec *WVec = bfRealArrayGetVecView(W);
  bfMatScaleCols(K, WVec);

  bfMatAddInplace(K, oneHalfEye);

  printf("- set up dense kernel matrix [%0.2fs]\n", bfToc());

  bfMatSave(K, "K.bin");

  /** Set up the RHS for the scattering problem: */

  bfToc();

  BfMat *rhs = NULL;
  { BfMatDenseComplex *_ = bfMatDenseComplexNew();
    bfMatDenseComplexInit(_, n, 1);
    for (BfSize i = 0; i < n; ++i) {
      BfPoint2 x;
      bfPoints2Get(X, i, x);
      *(_->data + i*_->rowStride) = -cexp(1i*k*(d[0]*x[0] + d[1]*x[1]));
    }
    rhs = bfMatDenseComplexToMat(_); }

  printf("- set up RHS for scattering problem [%0.2fs]\n", bfToc());

  /** Set up block LU preconditioner: */

  bfToc();

  assert(false); // TODO: implement using stuff in <bf/lu.h>

  printf("- set up block LU preconditioner [%0.2fs]\n", bfToc());

  /** Solve the dense system using GMRES: */

  bfToc();

  BfSize numIterDense;
  BfMat *sigmaDense = bfSolveGMRES(K, rhs, NULL, tol, n, &numIterDense, NULL);

  printf("- solve dense system using GMRES in %lu iterations [%0.2fs]\n", numIterDense, bfToc());

  /** Check BCs: */

  BfMat *KCheck = bfGetHelm2KernelMatrix(X, X, NULL, NULL, k, BF_LAYER_POTENTIAL_SINGLE);
  BfMat *DCheck = bfGetHelm2KernelMatrix(X, X, N, NULL, k, BF_LAYER_POTENTIAL_PV_DOUBLE);
  bfMatScale(KCheck, -1j*k);
  bfMatAddInplace(KCheck, DCheck);
  bfMatScaleCols(KCheck, WVec);

  BfMat *uDenseCheck = NULL;
  { BfMatDenseComplex *_ = bfMatDenseComplexNew();
    bfMatDenseComplexInit(_, n, 1);
    for (BfSize i = 0; i < n; ++i) {
      BfPoint2 x;
      bfPoints2Get(X, i, x);
      *(_->data + i*_->rowStride) = cexp(1i*k*(d[0]*x[0] + d[1]*x[1]));
    }
    uDenseCheck = bfMatDenseComplexToMat(_); }

  BfMat *uScatDenseCheck = bfMatMul(KCheck, sigmaDense);

  bfMatAddInplace(uDenseCheck, uScatDenseCheck);

  bfMatSave(uDenseCheck, "uDenseCheck.bin");

  /** Set up evaluation grid: */

  bfToc();

  BfPoints2 *XEval = bfPoints2NewGrid(&opts.bboxEval, opts.nxEval, opts.nyEval);

  printf("- set up %lu x %lu evaluation grid [%0.2fs]\n", opts.nxEval, opts.nyEval, bfToc());

  bfSavePoints2(XEval, "XEval.bin");

  /** Compute uIn on evaluation grid: */

  bfToc();

  BfMat *uIn = NULL;
  { BfSize nEval = opts.nxEval*opts.nyEval;
    BfMatDenseComplex *_ = bfMatDenseComplexNew();
    bfMatDenseComplexInit(_, nEval, 1);
    for (BfSize i = 0; i < nEval; ++i) {
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
  bfMatScale(KEval, -1j*k);
  bfMatAddInplace(KEval, DEval);
  bfMatScaleCols(KEval, WVec);

  BfMat *uScatDense = bfMatMul(KEval, sigmaDense);

  printf("- computed scattered field for dense problem [%0.2fs]\n", bfToc());

  bfMatSave(uScatDense, "uScatDense.bin");

  /** Set up quadtree for building multilevel BF: */

  bfToc();

  BfQuadtree quadtree;
  bfQuadtreeInit(&quadtree, X, N);

  printf("- set up quadtree [%0.2fs]\n", bfToc());

  BfTree *tree = bfQuadtreeToTree(&quadtree);
  BfPerm const *perm = bfTreeGetPermConst(tree);
  BfPerm revPerm = bfPermGetReversePerm(perm);

  (void)revPerm;

  /** Solve using HF-FMM (fmm2d): */

  // TODO: starting with a quick test...
  {
    BfReal eps = 1e-10;
    BfComplex zk = k;
    BfSize ns = bfSizeArrayGet(ellipseOffsets, 1) - bfSizeArrayGet(ellipseOffsets, 0);
    BfReal const *sources = (BfReal const *)&X->data[0];
    BfComplex const *charge = &bfMatToMatDenseComplex(rhs)->data[0];
    BfComplex *pot = malloc(ns*sizeof(BfComplex));
    int64_t ier;
    hfmm2d_s_c_p(eps, zk, ns, sources, charge, pot, &ier);
  }
}
