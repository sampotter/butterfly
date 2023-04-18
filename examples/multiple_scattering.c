#include <stdlib.h>

#include <bf/rand.h>

#include "multiple_scattering_context.h"

int main(int argc, char *argv[]) {
  Opts opts;
  bool success = parseArgs(&opts, argc, argv);
  if (!success) exit(EXIT_SUCCESS);

  MultipleScatteringContext context;
  init(&context, &opts);

  printInfo(&context);

  /* Common setup: */
  setUpGeometry(&context);
  setUpDiscretization(&context);
  buildQuadtrees(&context);
  setUpKrWorkspace(&context);
  setUpRhs(&context);

  if (shouldAssembleDenseK(&opts))
    assembleDenseK(&context);

  if (shouldAssembleButterfliedK(&opts))
    assembleButterfliedK(&context);

  if (shouldExtractDenseButterfliedK(&opts))
    extractDenseButterfliedK(&context);

  if (shouldAssemblePreconditioner(&opts))
    assemblePreconditioner(&context);

  if (shouldAssembleFmmK(&opts))
    assembleFmmK(&context);

  if (shouldComputeLu(&opts))
    computeLu(&context);

  if (shouldSolveLu(&opts))
    solveLu(&context);

  if (shouldSolveDenseGmres(&opts))
    solveDenseGmres(&context);

  if (shouldSolveDensePreconditionedGmres(&opts))
    solveDensePreconditionedGmres(&context);

  if (shouldSolveButterflyGmres(&opts))
    solveButterflyGmres(&context);

  if (shouldSolveButterflyPreconditionedGmres(&opts))
    solveButterflyPreconditionedGmres(&context);

  if (shouldSolveFmmGmres(&opts))
    solveFmmGmres(&context);

  if (shouldSolveFmmPreconditionedGmres(&opts))
    solveFmmPreconditionedGmres(&context);

  if (shouldCollectAndPrintStats(&opts))
    collectAndPrintStats(&context);

  if (shouldDoPostprocessing(&opts))
    doPostprocessing(&context);

  // deinit(&context, &opts);
}
