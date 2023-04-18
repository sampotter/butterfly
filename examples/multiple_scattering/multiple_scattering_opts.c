#include "multiple_scattering_opts.h"

#include <bf/const.h>
#include <bf/rand.h>

#include <assert.h>

#include "argtable3.h"

int const MAX_NUM_ARG_ERRORS = 20;

bool parseArgs(Opts *opts, int argc, char *argv[]) {
  bool success = true;

  /** All CLI parameters: */

  struct arg_lit *help;

  struct arg_int *seed;

  /* Problem parameters: */
  struct arg_dbl *wavenumber;
  struct arg_dbl *minDist;
  struct arg_dbl *axisLow;
  struct arg_dbl *axisHigh;
  struct arg_dbl *alphaRe;
  struct arg_dbl *alphaIm;
  struct arg_dbl *betaRe;
  struct arg_dbl *betaIm;
  struct arg_dbl *dx;
  struct arg_dbl *dy;

  /* Discretization parameters: */
  struct arg_dbl *h;
  struct arg_dbl *tol;
  struct arg_int *orderKR;

  /* Postprocessing parameters: */
  struct arg_dbl *xminEval;
  struct arg_dbl *xmaxEval;
  struct arg_dbl *yminEval;
  struct arg_dbl *ymaxEval;
  struct arg_int *nxEval;
  struct arg_int *nyEval;

  /* Which problems to solve: */
  struct arg_lit *doLu;
  struct arg_lit *doDenseGmres;
  struct arg_lit *doDensePreconditionedGmres;
  struct arg_lit *doButterflyGmres;
  struct arg_lit *doButterflyPreconditionedGmres;
  struct arg_lit *doFmmGmres;
  struct arg_lit *doFmmPreconditionedGmres;

  /* Miscellaneous options: */
  struct arg_lit *extractDenseButterfliedK;
  struct arg_lit *doPostprocessing;

  struct arg_end *end;

  /** Set up argtable and parse: */

  void *argtable[] = {
    help = arg_litn(NULL, "help", 0, 1, "Display help and exit"),

    seed = arg_intn(NULL, "seed", NULL, 0, 1, "Seed for random number generator (default: 0)"),

    /* Problem parameters: */
    wavenumber = arg_dbln("k", "wavenumber", "<k>", 0, 1, "wavenumber for the problem"),
    minDist = arg_dbln("r", "minDist", "<r>", 0, 1, "minimum distance between ellipse centers"),
    axisLow = arg_dbln("a", "axisLow", "<a>", 0, 1, "lower bound for uniformly sampled ellipse axes"),
    axisHigh = arg_dbln("b", "axisHigh", "<b>", 0, 1, "upper bound for uniformly sampled ellipse axes"),
    alphaRe = arg_dbln(NULL, "alphaRe", NULL, 0, 1, "real part of alpha (default: alpha = -1j*k)"),
    alphaIm = arg_dbln(NULL, "alphaIm", NULL, 0, 1, "imaginary part of alpha (default: alpha = -1j*k)"),
    betaRe = arg_dbln(NULL, "betaRe", NULL, 0, 1, "real part of beta (default: beta = -1)"),
    betaIm = arg_dbln(NULL, "betaIm", NULL, 0, 1, "imaginary part of beta (default: beta = -1)"),
    dx = arg_dbln(NULL, "dx", NULL, 0, 1, "x component of the wave direction vector (default: random)"),
    dy = arg_dbln(NULL, "dy", NULL, 0, 1, "y component of the wave direction vector (default: random)"),

    /* Discretization parameters: */
    h = arg_dbln("h", NULL, "<h>", 0, 1, "mesh fineness"),
    tol = arg_dbln(NULL, "tol", NULL, 0, 1, "relative tol used for low-rank approx (default: 1e-13)"),
    orderKR = arg_intn(NULL, "orderKR", NULL, 0, 1, "order of Kapur-Rokhlin quadrature (default: 6)"),

    /* Postprocessing parameters: */
    xminEval = arg_dbln(NULL, "xmin", "<xmin>", 0, 1, "The minimum x coord. of the eval. box"),
    xmaxEval = arg_dbln(NULL, "xmax", "<xmax>", 0, 1, "The maximum x coord. of the eval. box"),
    yminEval = arg_dbln(NULL, "ymin", "<ymin>", 0, 1, "The minimum y coord. of the eval. box"),
    ymaxEval = arg_dbln(NULL, "ymax", "<ymax>", 0, 1, "The maximum y coord. of the eval. box"),
    nxEval = arg_intn(NULL, "nx", "<nx>", 0, 1, "The number of eval. box nodes in x direction"),
    nyEval = arg_intn(NULL, "ny", "<ny>", 0, 1, "The number of eval. box nodes in y direction"),

    /* Which problems to solve: */
    doLu = arg_litn(NULL, "doLu", 0, 1, ""),
    doDenseGmres = arg_litn(NULL, "doDenseGmres", 0, 1, ""),
    doDensePreconditionedGmres = arg_litn(NULL, "doDensePreconditionedGmres", 0, 1, ""),
    doButterflyGmres = arg_litn(NULL, "doButterflyGmres", 0, 1, ""),
    doButterflyPreconditionedGmres = arg_litn(NULL, "doButterflyPreconditionedGmres", 0, 1, ""),
    doFmmGmres = arg_litn(NULL, "doFmmGmres", 0, 1, ""),
    doFmmPreconditionedGmres = arg_litn(NULL, "doFmmPreconditionedGmres", 0, 1, ""),

    /* Miscellaneous options: */
    extractDenseButterfliedK = arg_litn(NULL, "extractDenseButterfliedK", 0, 1, ""),
    doPostprocessing = arg_litn(NULL, "doPostprocessing", 0, 1, ""),

    end = arg_end(MAX_NUM_ARG_ERRORS)
  };

  /** Set default values then parse arguments: */

  *seed->ival = 0;

  /* We set alpha to NAN by default to signal that we should set alpha
   * based on the wavenumber: */
  *alphaRe->dval = BF_NAN;
  *alphaIm->dval = BF_NAN;

  *betaRe->dval = -1.0;
  *betaIm->dval = 0.0;

  *dx->dval = BF_NAN;
  *dy->dval = BF_NAN;

  *tol->dval = 1e-13;
  *orderKR->ival = 6;

  BfSize numErrors = arg_parse(argc, argv, argtable);

  /** Handle help and errors: */

  if (help->count > 0) {
    printf("Usage: %s", argv[0]);
    arg_print_syntax(stdout, argtable, "\n");
    printf("Run a test problem demonstrating multiple scattering\n");
    arg_print_glossary(stdout, argtable, "  %-25s %s\n");
    success = false;
    goto cleanup;
  }

  if (numErrors > 0) {
    arg_print_errors(stdout, end, argv[0]);
    printf("Try '%s --help' for more information.\n", argv[0]);
    success = false;
    goto cleanup;
  }

  /** Extract parameters from parsed CLI options: */

  opts->seed = *seed->ival;

  /* Seed the PRNG now! */
  bfSeed(opts->seed);

  opts->wavenumber = *wavenumber->dval;
  opts->minDist = *minDist->dval;
  opts->axisLow = *axisLow->dval;
  opts->axisHigh = *axisHigh->dval;

  if (!(isnan(*alphaRe->dval) && isnan(*alphaIm->dval))
      && (!isfinite(*alphaRe->dval) || !isfinite(*alphaIm->dval))) {
    printf("Got bad value: alpha = %g + i*%g\n", *alphaRe->dval, *alphaIm->dval);
    success = false;
    goto cleanup;
  }
  opts->alpha = isnan(*alphaRe->dval) && isnan(*alphaIm->dval) ?
    -1j*opts->wavenumber :
    *alphaRe->dval + *alphaIm->dval*1j;

  if (!(isnan(*betaRe->dval) && isnan(*betaIm->dval))
      && (!isfinite(*betaRe->dval) || !isfinite(*betaIm->dval))) {
    printf("Got bad value: beta = %g + i*%g\n", *betaRe->dval, *betaIm->dval);
    success = false;
    goto cleanup;
  }
  opts->beta = isnan(*betaRe->dval) && isnan(*betaIm->dval) ?
    -1.0 :
    *betaRe->dval + *betaIm->dval*1j;

  if (!(isnan(*dx->dval) && isnan(*dy->dval))
      && (!isfinite(*dx->dval) || !isfinite(*dy->dval))) {
    printf("Got bad value: d = (%g, %g)\n", *dx->dval, *dy->dval);
    success = false;
    goto cleanup;
  }
  if (isnan(*dx->dval) && isnan(*dy->dval)) {
    bfSampleRandomUnitVector2(opts->d);
  } else {
    opts->d[0] = *dx->dval;
    opts->d[1] = *dy->dval;
  }

  opts->h = *h->dval;

  if (*tol->dval <= 0) {
    printf("Got bad value: tol = %g\n", *tol->dval);
    success = false;
    goto cleanup;
  }
  opts->tol = *tol->dval;

  if (*orderKR->ival != 2 && *orderKR->ival != 6 && *orderKR->ival != 10) {
    printf("Got bad value: orderK = %d (should be 2, 6, or 10)\n", *orderKR->ival);
    success = false;
    goto cleanup;
  }
  opts->orderKR = *orderKR->ival;

  /* Postprocessing parameters: */
  opts->bboxEval.min[0] = *xminEval->dval;
  opts->bboxEval.max[0] = *xmaxEval->dval;
  opts->bboxEval.min[1] = *yminEval->dval;
  opts->bboxEval.max[1] = *ymaxEval->dval;
  opts->nxEval = *nxEval->ival;
  opts->nyEval = *nyEval->ival;

  /* Which problems to solve: */
  opts->doLu = doLu->count;
  opts->doDenseGmres = doDenseGmres->count;
  opts->doDensePreconditionedGmres = doDensePreconditionedGmres->count;
  opts->doButterflyGmres = doButterflyGmres->count;
  opts->doButterflyPreconditionedGmres = doButterflyPreconditionedGmres->count;
  opts->doFmmGmres = doFmmGmres->count;
  opts->doFmmPreconditionedGmres = doFmmPreconditionedGmres->count;

  /* Miscellaneous options: */
  opts->extractDenseButterfliedK = extractDenseButterfliedK->count;
  opts->doPostprocessing = doPostprocessing->count;

cleanup:
  arg_freetable(argtable, sizeof(argtable)/sizeof(argtable[0]));

  return success;
}

bool shouldAssembleDenseK(Opts const *opts) {
  return opts->doLu || opts->doDenseGmres || opts->doDensePreconditionedGmres;
}

bool shouldAssembleButterfliedK(Opts const *opts) {
  return opts->doButterflyGmres || opts->doButterflyPreconditionedGmres;
}

bool shouldExtractDenseButterfliedK(Opts const *opts) {
  return opts->extractDenseButterfliedK;
}

bool shouldAssembleFmmK(Opts const *opts) {
  return opts->doFmmGmres || opts->doFmmPreconditionedGmres;
}

bool shouldAssemblePreconditioner(Opts const *opts) {
  return opts->doDensePreconditionedGmres
    || opts->doButterflyPreconditionedGmres
    || opts->doFmmPreconditionedGmres;
}

bool shouldComputeLu(Opts const *opts) {
  return opts->doLu;
}

bool shouldSolveLu(Opts const *opts) {
  return opts->doLu;
}

bool shouldSolveDenseGmres(Opts const *opts) {
  return opts->doDenseGmres;
}

bool shouldSolveDensePreconditionedGmres(Opts const *opts) {
  return opts->doDensePreconditionedGmres;
}

bool shouldSolveButterflyGmres(Opts const *opts) {
  return opts->doButterflyGmres;
}

bool shouldSolveButterflyPreconditionedGmres(Opts const *opts) {
  return opts->doButterflyPreconditionedGmres;
}

bool shouldSolveFmmGmres(Opts const *opts) {
  return opts->doFmmGmres;
}

bool shouldSolveFmmPreconditionedGmres(Opts const *opts) {
  return opts->doFmmPreconditionedGmres;
}

bool shouldCollectAndPrintStats(Opts const *opts) {
  return opts->doLu || opts->doDenseGmres || opts->doDensePreconditionedGmres
    || opts->doButterflyGmres || opts->doButterflyPreconditionedGmres
    || opts->doFmmGmres || opts->doFmmPreconditionedGmres;
}

bool shouldDoPostprocessing(Opts const *opts) {
  return opts->doPostprocessing;
}
