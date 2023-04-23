#pragma once

#include <bf/def.h>
#include <bf/bbox.h>

typedef struct {
  BfSize seed;

  /** Parameters which get forwarded to MultipleScatteringContext: */

  /* Problem parameters: */
  BfReal wavenumber;
  BfReal minDist;
  BfReal axisLow;
  BfReal axisHigh;
  BfComplex alpha;
  BfComplex beta;
  BfVector2 d;

  /* Discretization parameters: */
  BfReal h;
  BfReal tol;
  BfSize orderKR;

  /* Postprocessing parameters: */
  BfBbox2 bboxEval;
  BfSize nxEval;
  BfSize nyEval;

  /** "Job control" parameters: */

  /* Which problems to solve: */
  bool doLu;
  bool doDenseGmres;
  bool doDensePreconditionedGmres;
  bool doButterflyGmres;
  bool doButterflyPreconditionedGmres;
  bool doFmmGmres;
  bool doFmmPreconditionedGmres;

  /* Miscellaneous options: */
  bool extractDenseButterfliedK;
  bool doPostprocessing;
} Opts;

bool parseArgs(Opts *opts, int argc, char *argv[]);

bool shouldAssembleDenseK(Opts const *opts);
bool shouldAssembleButterfliedK(Opts const *opts);
bool shouldExtractDenseButterfliedK(Opts const *opts);
bool shouldAssembleFmmK(Opts const *opts);
bool shouldEstimateFmmTol(Opts const *opts);
bool shouldAssemblePreconditioner(Opts const *opts);
bool shouldComputeLu(Opts const *opts);
bool shouldSolveLu(Opts const *opts);
bool shouldSolveDenseGmres(Opts const *opts);
bool shouldSolveDensePreconditionedGmres(Opts const *opts);
bool shouldSolveButterflyGmres(Opts const *opts);
bool shouldSolveButterflyPreconditionedGmres(Opts const *opts);
bool shouldSolveFmmGmres(Opts const *opts);
bool shouldSolveFmmPreconditionedGmres(Opts const *opts);
bool shouldCollectAndPrintStats(Opts const *opts);
bool shouldDoPostprocessing(Opts const *opts);
